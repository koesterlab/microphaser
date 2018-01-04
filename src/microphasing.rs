use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;
use std::str;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;
use rust_htslib::{bam, bcf};
use rust_htslib::bam::record::Cigar;
use rust_htslib::prelude::*;

use common::{Gene, Variant, Interval, Transcript};


pub fn bitvector_is_set(b: u64, k: u32) -> bool {
    (b & (1 << k)) != 0
}


pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        Variant::SNV { pos: pos, alt: alt } => {
            let b = read.seq()[read.cigar().read_pos(pos)?.unwrap()];
            Ok(b == alt)
        },
        Variant::Insertion { pos: pos, seq: seq } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar() {
                match c {
                    Cigar::Ins(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        },
        Variant::Deletion { pos: pos, len: varlen } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar() {
                match c {
                    Cigar::Del(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        }
    }
}


pub struct ObservationMatrix<'a> {
    inner: VecDeque<u64>,
    reads: BTreeMap<u32, Vec<&'a bam::Record>>,
    variants: VecDeque<Variant>
}


impl<'a> ObservationMatrix<'a> {
    pub fn shrink_left(&mut self, k: u32) {
        self.variants.drain(..k);
        for b in self.inner.iter_mut() {
            *b = *b & (2.pow(self.ncols()) - 1);
        }
    }

    pub fn extend_right<I: Iterator<Item=Variant>>(&mut self, variants: I) {
        let k = variants.len();
        for b in self.inner.iter_mut() {
            *b = *b << k;
        }
        for (read, b) in self.reads.iter().zip(self.inner.iter_mut()) {
            for (i, variant) in variants.iter().rev().enumerate() {
                if supports_variant(read, variant) {
                    b |= 1 << i;
                }
            }
        }
        self.variants.extend(variants);
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u32) {
        let reads = self.reads.split_off(interval_end);
        self.inner.drain(..self.reads.len());
        self.reads = reads;
    }

    /// Add read, while considering given interval end.
    pub fn push_read(&mut self, read: &bam::Record, interval_end:u32) {
        let end_pos = read.cigar().end_pos();
        if end_pos >= interval_end {
            // only insert if end_pos is larger than the interval end
            self.reads.entry(end_pos).or_insert_with(|| Vec::new()).push(read);
            self.inner.push(0);
        }
    }

    pub fn ncols(&self) -> u32 {
        self.variants.len() as u32
    }

    pub fn print_haplotypes<O: io::Write>(
        &self,
        gene: &Gene,
        offset: u32,
        window_len: u32,
        refseq: &[u8],
        fasta_writer: &mut fasta::Writer<O>
    ) -> Result<(), Box<Error>> {
        let variants = self.variants.iter().collect_vec();
        // count haplotypes
        let mut haplotypes = VecMap::new();
        for b in self.inner {
            haplotypes.entry(b).or_insert(0) += 1;
        }
        for (haplotype, count) in haplotypes.iter() {
            // build haplotype sequence
            let mut seq = Vec::with_capacity(window_len);
            let freq = *count as f64 / self.inner.len();
            let mut i = offset;
            let mut j = 0;
            while i < offset + window_len {
                // TODO what happens if an insertion starts upstream of window and overlaps it
                while i == variants[j].pos {
                    if bitvector_is_set(haplotype, j) {
                        match variants[j] {
                            Variant::SNV { alt: a } => {
                                seq.push(a);
                                i += 1;
                            },
                            Variant::Insertion { seq: s } => {
                                seq.extend(s);
                                i += 1;
                            },
                            Variant::Deletion { len: l } => i += l
                        }
                        break;
                    }
                }
                seq.push(refseq[i - gene.start()]);
                i += 1
            }
            // restrict to window len (it could be that we insert too much above)
            seq = seq[..window_len];

            fasta_writer.write(
                format!("{}:offset={},af={:.2}", gene.id, offset, freq),
                None,
                &seq
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
    fasta_reader.read(&gene.chrom, gene.start(), gene.end(), refseq);
    for transcript in gene.transcripts {
        let observations = ObservationMatrix::new();

        for exon in transcript.exons {
            let mut offset = exon.start;

            while offset + window_len < exon.end {
                // advance window to next position
                let (added_vars, deleted_vars) = variant_buffer.fetch(
                    &gene.chrom.as_bytes(), offset, offset + window_len
                )?;
                let (added_reads, deleted_reads) = read_buffer.fetch(
                    &gene.chrom.as_bytes(), offset, offset + window_len
                )?;

                // delete rows
                observations.cleanup_reads(offset + window_len);

                // delete columns
                observations.shrink_left(deleted_vars);

                // add new reads
                for read in read_buffer.iter() {
                    observations.push_read(read, offset + window_len);
                }

                // collect variants
                let variants = variant_buffer.iter().skip(
                    variant_buffer.len() - added_vars
                ).map(|rec| Variant::new(rec).unwrap()).flatten();
                // add columns
                observations.extend_right(variants);

                // print haplotypes
                observations.print_haplotypes(gene, offset, window_len, refseq, fasta_writer)?;

                offset += 3;
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
    let phase_last_gene = || -> Result<(), Box<Error>> {
        if let Some(gene) = gene {
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
                phase_last_gene()?;

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
            "transcript" if gene.is_some() => {
                // register new transcript
                gene.unwrap().transcripts.push(
                    Transcript::new(record.attributes().get("transcript_id").expect(
                        "missing transcript_id attribute in GTF"
                    ))
                );
            },
            "exon" if gene.is_some() => {
                // register exon
                gene.unwrap()
                    .transcripts.last()
                    .expect("no transcript record before exon in GTF")
                    .exons.push(Interval::new(
                    *record.start() as u32,
                    *record.end() as u32
                ));
            },
            _ => continue
        }
    }
    phase_last_gene()?;
    Ok(())
}
