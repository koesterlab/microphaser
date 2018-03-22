use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;
use rust_htslib::{bam, bcf};
use rust_htslib::bam::record::Cigar;

use common::{Gene, Variant, Interval, Transcript};


pub fn bitvector_is_set(b: u64, k: usize) -> bool {
    (b & (1 << k)) != 0
}


pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        &Variant::SNV { pos, alt, .. } => {
            let b = read.seq()[read.cigar().read_pos(pos, false, false)?.expect("bug: read does not enclose variant") as usize];
            Ok(b == alt)
        },
        &Variant::Insertion {..} => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                match c {
                    &Cigar::Ins(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        },
        &Variant::Deletion {..} => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                match c {
                    &Cigar::Del(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        }
    }
}


pub struct Observation{
    read: bam::Record,
    haplotype: u64
}


pub struct ObservationMatrix {
    observations: BTreeMap<u32, Vec<Observation>>,
    variants: VecDeque<Variant>
}


impl ObservationMatrix {
    pub fn new() -> Self {
        ObservationMatrix {
            observations: BTreeMap::new(),
            variants: VecDeque::new()
        }
    }

    pub fn shrink_left(&mut self, k: usize) {
        self.variants.drain(..k);
        let mask = 2u64.pow(self.ncols()) - 1;
        for obs in self.observations.values_mut().flatten() {
            obs.haplotype = obs.haplotype & mask;
        }
    }

    pub fn extend_right(
        &mut self, mut variants: Vec<Variant>,
        new_variant_count:usize, old_variant_count:usize
    ) -> Result<(), Box<Error>> {
        let k = new_variant_count;
        if k > 0 {
            for obs in self.observations.values_mut().flatten() {
                obs.haplotype <<= k;
            }
        }
        for obs in self.observations.values_mut().flatten() {
            for (i, variant) in variants.iter().rev().enumerate() {
                if (obs.read.pos() as u32) > variant.pos() {
                continue
                }
                if supports_variant(&obs.read, &variant)? {
                    debug!("Read {} supports the variant at {}", String::from_utf8_lossy(obs.read.qname()), variant.pos());
                    obs.haplotype |= 1 << i;
                    debug!("Haplotype {}", obs.haplotype);
                }
            }
        }
        let new_variants = variants.split_off(old_variant_count);
        self.variants.extend(new_variants.into_iter());

        Ok(())
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u32) {
        let observations = self.observations.split_off(&interval_end);
        self.observations = observations;
    }

    pub fn already_in(&mut self, read: bam::Record) -> Result<bool, Box<Error>> {
        let end_pos = read.cigar().end_pos()? as u32;
        let qname = read.qname();
        if self.observations.contains_key(&end_pos) {
            for obs in self.observations.get(&end_pos).unwrap() {
                if obs.read.qname() == qname {
                    return Ok(true)
                }
            }
           return Ok(false)
        }
        Ok(false)
    }

    /// Add read, while considering given interval end. TODO: Think about reads that are not overlapping variant
    pub fn push_read(&mut self, read: bam::Record, interval_end:u32, interval_start:u32) -> Result<(), Box<Error>> {
        let end_pos = read.cigar().end_pos()? as u32;
        let start_pos = read.pos() as u32;
        if end_pos >= interval_end && start_pos <= interval_start {
            // only insert if end_pos is larger than the interval end
            self.observations.entry(end_pos).or_insert_with(|| Vec::new()).push(
                Observation { read: read, haplotype: 0 }
            );
        }
        Ok(())
    }

    pub fn ncols(&self) -> u32 {
        self.variants.len() as u32
    }

    pub fn nrows(&self) -> usize {
        self.observations.values().map(|o| o.len()).sum()
    }

    pub fn print_haplotypes<O: io::Write>(
        &self,
        gene: &Gene,
        transcript: &Transcript,
        offset: u32,
        window_len: u32,
        refseq: &[u8],
        fasta_writer: &mut fasta::Writer<O>
    ) -> Result<(), Box<Error>> {
        let variants = self.variants.iter().collect_vec();
        // count haplotypes
        let mut haplotypes: VecMap<usize> = VecMap::new();
        for obs in self.observations.values().flatten() {
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let mut seq = Vec::with_capacity(window_len as usize);
        for (haplotype, count) in haplotypes.iter() {
            // VecMap forces usize as type for keys, but our haplotypes as u64
            let haplotype = haplotype as u64;
            // build haplotype sequence
            seq.clear();
            let mut is_germline = true;
            let mut is_variant = false;
            let freq = *count as f64 / self.nrows() as f64;
            let mut i = offset;
            let mut j = 0;
            if variants.is_empty() {
                seq.extend(&refseq[(offset - gene.start()) as usize..(offset + window_len - gene.start()) as usize]);
            } else {
                while i < offset + window_len {
                    // TODO what happens if an insertion starts upstream of window and overlaps it
                    while j < variants.len() && i == variants[j].pos() {
                        if bitvector_is_set(haplotype, j) {
                            match variants[j] {
                                &Variant::SNV { alt, .. } => {
                                    seq.push(alt.to_ascii_lowercase());
                                    i += 1;
                                },
                                &Variant::Insertion { seq: ref s, .. } => {
                                    seq.extend(s.to_ascii_lowercase().into_iter());
                                    i += 1;
                                },
                                &Variant::Deletion { len, .. } => i += len
                            }
                            is_germline = variants[j].is_germline();
                            is_variant = true;
                            break;
                        } else {
                            j += 1;
                        }
                    }
                    seq.push(refseq[(i - gene.start()) as usize]);
                    i += 1
                }
            }

            fasta_writer.write(
                &format!(
                    "{}:{{\"offset\":{},\"af\":{:.2},\"variant\":{},\"somatic\":{}}}",
                    transcript.id, offset, freq, is_variant, !is_germline
                ),
                None,
                // restrict to window len (it could be that we insert too much above)
                &seq[..window_len as usize]
            )?;
        }
        Ok(())
    }
}


pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    gene: &Gene,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::RecordBuffer,
    fasta_writer: &mut fasta::Writer<O>,
    window_len: u32,
    refseq: &mut Vec<u8>
) -> Result<(), Box<Error>> {
    fasta_reader.read(&gene.chrom, gene.start() as u64, gene.end() as u64, refseq)?;
    for transcript in &gene.transcripts {
        let mut observations = ObservationMatrix::new();

        for exon in &transcript.exons {
            let mut offset = exon.start;
            while offset + window_len < exon.end {
                // advance window to next position
                let (added_vars, deleted_vars) = variant_buffer.fetch(
                    &gene.chrom.as_bytes(), offset, offset + window_len
                )?; // -1 because gtf is 1-based, bcf is 0-based
                read_buffer.fetch(
                    &gene.chrom.as_bytes(), offset, offset + window_len
                )?;

                {
                    // delete rows
                    observations.cleanup_reads(offset + window_len);

                    // delete columns
                    observations.shrink_left(deleted_vars);

                    // add new reads
                    for read in read_buffer.iter() {
                        if observations.already_in(read.clone())? == false {
                            observations.push_read(read.clone(), offset + window_len, offset)?;
                        }
                    }
                    // TODO do not re-annotate old observations
                    // 1. annotate new observations with old variants
                    // 2. annotate all observations with new variants

                    // collect variants
                    let nvars = variant_buffer.len();
                    let variants = variant_buffer.iter_mut()//.skip(
                        //nvars - added_vars
                    .map(|rec| Variant::new(rec).unwrap()).flatten().collect_vec();

                    // add columns
                    observations.extend_right(variants, added_vars, nvars - added_vars)?;

                    // print haplotypes
                    observations.print_haplotypes(
                        gene, transcript, offset, window_len, refseq, fasta_writer
                    )?;

                    offset += 3;
                }
            }
        }
    }
    Ok(())
}


pub fn phase<F: io::Read + io::Seek, G: io::Read, O: io::Write>(
    fasta_reader: &mut fasta::IndexedReader<F>,
    gtf_reader: &mut gff::Reader<G>,
    bcf_reader: bcf::Reader,
    bam_reader: bam::IndexedReader,
    fasta_writer: &mut fasta::Writer<O>,
    window_len: u32
) -> Result<(), Box<Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader);
    let mut variant_buffer = bcf::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence

    let mut gene = None;
    let mut phase_last_gene = |gene| -> Result<(), Box<Error>> {
        if let Some(ref gene) = gene {
            phase_gene(
                &gene, fasta_reader, &mut read_buffer,
                &mut variant_buffer, fasta_writer,
                window_len,
                &mut refseq
            )?;
        }
        Ok(())
    };

    for record in gtf_reader.records() {
        let record = record?;
        match record.feature_type() {
            "gene" => {
                // first, phase the last gene
                phase_last_gene(gene)?;
                if record.attributes().get("gene_biotype").expect("missing gene_biotype in GTF") == "protein_coding" {
                    // if protein coding, start new gene
                    gene = Some(Gene::new(
                        record.attributes().get("gene_id").expect("missing gene_id in GTF"),
                        record.seqname(),
                        Interval::new(*record.start() as u32, *record.end() as u32)
                    ));
                } else {
                    // ignore until protein coding gene occurs
                    gene = None
                }
            },
            "transcript" => {
                // register new transcript
                gene.as_mut()
                    .expect("no gene record before transcript in GTF").transcripts.push(
                    Transcript::new(record.attributes().get("transcript_id").expect(
                        "missing transcript_id attribute in GTF"
                    ))
                );
            },
            "exon" => {
                // register exon
                gene.as_mut().expect("no gene record before exon in GTF")
                    .transcripts.last_mut()
                    .expect("no transcript record before exon in GTF")
                    .exons.push(Interval::new(
                    *record.start() as u32,
                    *record.end() as u32
                ));
            },
            "start_codon" => {
                gene.as_mut().expect("no gene record before start_codon in GTF")
                    .transcripts.last_mut()
                    .expect("no transcript record before start codon in GTF")
                    .exons.last_mut()
                    .expect("no exon record before start codon in GTF")
                    .start = *record.start() as u32;
            },
            "stop_codon" => {
                gene.as_mut().expect("no gene record before stop_codon in GTF")
                    .transcripts.last_mut()
                    .expect("no transcript record before stop codon in GTF")
                    .exons.last_mut()
                    .expect("no exon record before stop codon in GTF")
                    .end = *record.end() as u32;
            }
            _ => continue
        }
    }
    phase_last_gene(gene)?;

    Ok(())
}
