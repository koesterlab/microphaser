use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;
//use std::cmp;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;
use bio::utils::Strand;
use rust_htslib::{bam, bcf};
use rust_htslib::bam::record::Cigar;


use common::{Gene, Variant, Interval, Transcript};


pub fn bitvector_is_set(b: u64, k: usize) -> bool {
    (b & (1 << k)) != 0
}


pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        &Variant::SNV { pos, alt, .. } => {
            eprintln!("Variant pos: {}", variant.pos());
            eprintln!("Read pos: {}", read.pos());
            eprintln!("Read end: {}", read.seq().len() + (read.pos() as usize));
            eprintln!("Read to check support: {}", String::from_utf8_lossy(read.qname()));
//            for c in read.cigar().iter() {
//                match c {
//                    &Cigar::Del(_) => return Ok(false),
//                    _ => ()
//                }
//            }
            let b = match read.cigar().read_pos(pos, false, false) {
                Ok(None) => return Ok(false),
                Ok(Some(p)) => read.seq()[p as usize],
                _ => return Ok(false)
            };
//            let b = read.seq()[read.cigar().read_pos(pos, false, false)?.expect("bug: read does not enclose variant") as usize];
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


#[derive(Debug)]
pub struct Observation{
    read: bam::Record,
    haplotype: u64
}


impl Observation {
    pub fn update_haplotype(&mut self, i: usize, variant: &Variant) -> Result<(), Box<Error>> {
        eprintln!("Read pos {} ; variant pos {}", self.read.pos() as u32, variant.pos());
        if (self.read.pos() as u32) > variant.pos() {
            panic!("bug: read starts right of variant");
        }
        if supports_variant(&self.read, &variant)? {
//            eprintln!(
//                "Read {} supports the variant at {}",
//                String::from_utf8_lossy(self.read.qname()), variant.pos()
//            );
            self.haplotype |= 1 << i;
            //eprintln!("Haplotype {}", self.haplotype);
            //eprintln!("Haplotype {}", self.haplotype);
        }

        Ok(())
    }
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
        eprintln!("self.variants len {}", self.variants.len());
        eprintln!("range to drain: 0 - {}",k);
        self.variants.drain(..k);
        eprintln!("drained");
        let mask = 2u64.pow(self.ncols()) - 1;
        for obs in self.observations.values_mut().flatten() {
            obs.haplotype = obs.haplotype & mask;
        }
    }

    pub fn extend_right(
        &mut self, new_variants: Vec<Variant>
    ) -> Result<(), Box<Error>> {
        let k = new_variants.len();
        eprintln!("Extend variants!");
        eprintln!("New variants {}", k);
        if k > 0 {
            for obs in self.observations.values_mut().flatten() {
                obs.haplotype <<= k;
            }
        }
        for obs in self.observations.values_mut().flatten() {
            for (i, variant) in new_variants.iter().rev().enumerate() {
                obs.update_haplotype(i, variant)?;
            }
        }
        self.variants.extend(new_variants.into_iter());

        Ok(())
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u32, reverse: bool) {
        eprintln!("Number of reads(before removal): {}", self.observations.len());
        let observations = self.observations.split_off(&interval_end);
        if !reverse {
            self.observations = observations;
        }
        eprintln!("Number of reads after removal: {}",  self.observations.len());
    }

    pub fn contains(&mut self, read: &bam::Record) -> Result<bool, Box<Error>> {
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
    pub fn push_read(&mut self, read: bam::Record, interval_end:u32, interval_start:u32, reverse: bool) -> Result<(), Box<Error>> {
        let end_pos = read.cigar().end_pos()? as u32;
        let start_pos = read.pos() as u32;
        eprintln!("{} {} : {} {}", start_pos, end_pos, interval_start, interval_end);
        if end_pos >= interval_end && start_pos <= interval_start {
            eprintln!("Pushing read? {}-{}", start_pos, end_pos);
            // only insert if end_pos is larger than the interval end
            let mut obs = Observation { read: read, haplotype: 0 };
            for (i, variant) in self.variants.iter().enumerate() {
                obs.update_haplotype(i, variant)?;
            }
            if reverse {
                self.observations.entry(start_pos).or_insert_with(|| Vec::new()).push(obs);
            }
            else {
                self.observations.entry(end_pos).or_insert_with(|| Vec::new()).push(obs);
            }
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
        eprintln!("Observation length: {}", self.observations.values().len());
        for obs in self.observations.values().flatten() {
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let mut seq = Vec::with_capacity(window_len as usize);
        eprintln!("printing with offset: {}",offset);
        for (haplotype, count) in haplotypes.iter() {
            // VecMap forces usize as type for keys, but our haplotypes as u64
            let haplotype = haplotype as u64;
            // build haplotype sequence
            seq.clear();
            let mut n_somatic = 0;
            let mut n_variants = 0;
            let freq = *count as f64 / self.nrows() as f64;
            let mut i = offset;
            let mut j = 0;
            let mut window_end = offset + window_len;
            if variants.is_empty() {
                seq.extend(&refseq[(offset - gene.start()) as usize..(offset + window_len - gene.start()) as usize]);
            } else {
                while i < window_end {
                    //eprintln!("window_end: {}", window_end);
                    //eprintln!("i: {}", i);
                    // TODO what happens if a deletion starts upstream of window and overlaps it
                    while j < variants.len() && i == variants[j].pos() {
                        //eprintln!("j: {}", j);
                        if bitvector_is_set(haplotype, j) {
                            match variants[j] {
                                &Variant::SNV { alt, .. } => {
                                    seq.push(alt.to_ascii_lowercase());
                                    i += 1;
                                },
                                &Variant::Insertion { seq: ref s, .. } => {
                                    //eprintln!("Insertion Sequence length {}", s.len());
                                    seq.extend(s.to_ascii_lowercase().into_iter());
                                    i += 1;
                                    window_end -= (s.len() as u32) -1;
                                },
                                &Variant::Deletion { len, .. } => {
                                    i += len;
                                    window_end += len;
                                }
                            }
                            if !variants[j].is_germline() {
                                n_somatic += 1;
                            }
                            n_variants += 1;
                            j += 1;
                        } else {
                            j += 1;
                        }
                    }
                    seq.push(refseq[(i - gene.start()) as usize]);
                    //eprintln!("Sequence length: {}", seq.len());
                    i += 1
                }
            }
            eprintln!("Writing to fasta");
            fasta_writer.write(
                &format!(
                    "{}:{{\"offset\":{},\"af\":{:.2},\"variants\":{},\"somatic\":{}}}",
                    transcript.id, offset, freq, n_variants, n_somatic
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
    let mut variant_tree = BTreeMap::new();
    let mut read_tree = BTreeMap::new();
    eprintln!("Start Phasing");
    read_buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;

    //read_buffer.iter().map(|rec| read_tree.entry(rec.pos() as u32).or_insert_with(|| Vec::new()).push(rec));

    for rec in read_buffer.iter() {
        read_tree.entry(rec.pos() as u32).or_insert_with(|| Vec::new()).push(rec.clone())
    }
    eprintln!("Reads Tree len: {}", read_tree.len());

    let max_read_len = 101;
//    for (_p, v) in read_tree {
//        for r in v {
//             max_read_len = cmp::max(max_read_len, r.seq().len() as u32)
//        }
//    }

    let (_addedvars, _deletedvars) = variant_buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;
    let _vars = variant_buffer.iter_mut().map(|rec| variant_tree.insert(rec.pos(), Variant::new(rec).unwrap())).collect_vec();

//    for rec in variant_buffer.iter_mut() {
//        variant_tree.entry(rec.pos()).or_insert(Variant::new(rec).unwrap());
//    }

//    let count_vars = |range| btree.range(range).map(|var| var.1.len()).sum();
//    let genevars = buffer.iter().collect_vec();
//    for v in genevars {
//        btree.insert(v.pos(),*v);
//    }

    for transcript in &gene.transcripts {
        eprintln!("Transcript strand orientation: {:?}", transcript.strand);
        let mut observations = ObservationMatrix::new();
        let mut frameshifts = BTreeMap::new();
        frameshifts.insert(0, 0);
        //let mut next_exon_offset = 0;
        let mut exon_rest = 0;
        let mut last_window_vars = 0;
        for exon in &transcript.exons {
        eprintln!("Exon Start: {}", exon.start);
        eprintln!("Exon Ende: {}", exon.end);
            //eprintln!("Exon offset: {}", next_exon_offset);
            let current_exon_offset = match exon_rest {
                0 => 0,
                _ => 3 - exon_rest
            };
            exon_rest = 0;
            let mut offset = if transcript.strand == Strand::Reverse {
                exon.end - window_len - current_exon_offset
            } else {
                exon.start + current_exon_offset
            };
            let mut old_offset = offset;
            eprintln!("Variants left from previous Exon: {}", last_window_vars);
            observations.shrink_left(last_window_vars);
            last_window_vars = 0;
            loop {
                let valid = match transcript.strand {
                    Strand::Reverse => offset >= exon.start,
                    Strand::Forward => offset + window_len < exon.end,
                    Strand::Unknown => offset + window_len < exon.end
                };
                if !valid {
                    break;
                }
                //eprintln!("Old Offset {}", old_offset);
                //eprintln!("Offset {} + Windowlen: {}", offset, offset + window_len);
                // advance window to next position

                let nvars = variant_tree.range(offset..(offset + window_len)).map(|var| var.1).flatten().count();//.len()).sum();
                last_window_vars = nvars;
                eprintln!("Variants in window: {}",nvars);
                let mut added_vars = if offset == old_offset {
                    nvars
                } else if offset > old_offset {
                    variant_tree.range((old_offset + window_len)..(offset + window_len)).map(|var| var.1).flatten().count()
                } else {
                    variant_tree.range(offset..old_offset).map(|var| var.1).flatten().count()
                };

                let mut deleted_vars = if offset == old_offset {
                    0
                } else if offset > old_offset {
                    variant_tree.range(old_offset..offset).map(|var| var.1).flatten().count()
                } else {
                    variant_tree.range((offset + window_len)..(old_offset + window_len)).map(|var| var.1).flatten().count()
                };

//                let (added_vars, deleted_vars) = variant_buffer.fetch(
//                    &gene.chrom.as_bytes(), offset , offset + window_len -1
//                )?; // -1 because gtf is 1-based, bcf is 0-based

//                if transcript.strand == Strand::Reverse {
//                    read_buffer.fetch(
//                        &gene.chrom.as_bytes(), offset, offset + window_len
//                    )?;
//                }
//                else {
//                    read_buffer.fetch(
//                        &gene.chrom.as_bytes(), offset, offset + window_len
//                    )?;
//                }


                eprintln!("Offset: {} - max_read_len - window_len {}", offset, (max_read_len - window_len));
                let mut reads = if transcript.strand == Strand::Reverse {
                    if offset == exon.end - window_len - current_exon_offset {
                        read_tree.range((offset - (max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                    else {
                        read_tree.range((offset - (max_read_len - window_len))..(offset - (max_read_len - window_len) + 1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                }
                else {
                    if offset == exon.start + current_exon_offset {
                        read_tree.range((offset-(max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten().collect_vec()//offset+1?
                    }
                    else {
                        read_tree.range((offset-1)..(offset+1)).map(|rec| rec.1).flatten().collect_vec()//offset+1?
                    }
                };

//                let mut reads = match transcript.strand {
//                    Strand::Reverse => if offset == exon.end - window_len {
//                        read_tree.range((offset - (max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten()
//                    }
//                    else {
//                        read_tree.range((offset - (max_read_len - window_len))..(offset - (max_read_len - window_len) + 1)).map(|rec| rec.1).flatten()
//                    },
//                    _ => if offset == exon.start {
//                        read_tree.range((offset-(max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten()
//                    }
//                    else {
//                        read_tree.range((offset-1)..(offset+1)).map(|rec| rec.1).flatten()
//                    }
//                };

                eprintln!("Variants added: {}", added_vars);
                eprintln!("Variants deleted: {}", deleted_vars);
                {
                    // delete rows
                    if transcript.strand == Strand::Reverse {
                        observations.cleanup_reads(offset, true);
                    }
                    else {
                        observations.cleanup_reads(offset + window_len, false);
                    }
                    // delete columns
                    observations.shrink_left(deleted_vars);

                    // add new reads
                    eprintln!("Reads: {}", reads.len());
                    for read in reads {
                        //if !observations.contains(read)? {
                            //eprintln!("Add a read to observations! {}", read.pos() as u32);
                        if transcript.strand == Strand::Reverse {
                            observations.push_read(read.clone(), offset + window_len, offset, true)?;
                        }
                        else {
                            observations.push_read(read.clone(), offset + window_len, offset, false)?;
                        }
                        //}
                    }

//                    for (i, v) in variant_tree.range(offset..(offset + window_len)) {
//                        eprintln!("Variant pos {}", i);
//                        eprintln!("Variant: {}", v.flatten());
//                    }

                    // collect variants
                    //eprintln!("Var range len {}", variant_tree.range.len());
                    let variants = match transcript.strand {
                        Strand::Reverse => variant_tree.range_mut(offset..(offset + window_len)).rev()//.skip(
                            //nvars - added_vars)
                        .map(|var| var.1.clone()).flatten().skip(nvars - added_vars).collect_vec(),
                        _ => variant_tree.range_mut(offset..(offset + window_len)//).skip(
                            //nvars - added_vars
                        ).map(|var| var.1.clone()).flatten().skip(nvars - added_vars).collect_vec()
                    };
                    eprintln!("Variants(after deleting and adding): {}", variants.len());
                    // determine frameshifts
                    for variant in &variants {
                        eprintln!("Variants!");
                        eprintln!("Variant Pos: {}", variant.pos());
                        let s = variant.frameshift();
                        if s > 0 {
                            let previous = frameshifts.values().map(|prev| prev + s).collect_vec();
                            for s_ in previous {
                                frameshifts.insert(variant.end_pos(), s + s_);
                            }
                        }
                    }

                    // add columns

                    observations.extend_right(variants)?;

                    for (_, &frameshift) in frameshifts.range(..offset) {
                        //eprintln!("Offset {}", offset);
                        let coding_shift = match transcript.strand {
                                Strand::Forward => offset - exon.start,
                                Strand::Reverse => exon.end - offset,
                                Strand::Unknown => offset - exon.start//exon.end - offset
                        };
                        eprintln!("Offset - Exonstart % 3: {}", coding_shift % 3);
                        eprintln!("Shift: {}", frameshift + current_exon_offset);
                        if coding_shift % 3 == frameshift + current_exon_offset {
                            // print haplotypes
                            eprintln!("Should print haplotypes");
                            observations.print_haplotypes(
                                gene, transcript, offset, window_len, refseq, fasta_writer
                            )?;
                            exon_rest = match transcript.strand {
                                Strand::Forward => exon.end - (offset + window_len),
                                Strand::Reverse => offset - exon.start,
                                Strand::Unknown => exon.end - (offset + window_len)
                            };
                            eprintln!("Exon Rest {}", exon_rest);
//                            if rest < 3 && rest > 0 {
//                                next_exon_offset = 3 - rest;
//                                //break;
//                            }
                        }
                    }
                    //eprintln!("Variants in window: {}",nvars);
                    //offset += 1;
                    old_offset = offset;
                    match transcript.strand {
                        Strand::Reverse => offset -= 1,
                        Strand::Forward => offset += 1,
                        Strand::Unknown => offset += 1
                    }
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
    let mut variant_buffer = bcf::RecordBuffer::new(bcf_reader);//bcf::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence



    let mut gene = None;
    let mut phase_last_gene = | gene: Option<Gene>| -> Result<(), Box<Error>> {
        if let Some(ref gene) = gene {
            if gene.biotype == "protein_coding" {
//                let (_addedvars, _deletedvars) = variant_buffer.buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;
//                for variant in variant_buffer.buffer.iter() {
//                    variant_buffer.btree.entry(variant.pos()).or_insert(variant);
//                }
//                variant_buffer.buffer.iter().map(|rec| btree.insert(rec.pos(), rec));
                phase_gene(
                    &gene, fasta_reader, &mut read_buffer,
                    &mut variant_buffer, fasta_writer,
                    window_len,
                    &mut refseq
                )?;
            }
        }
        Ok(())
    };

    for record in gtf_reader.records() {
        let record = record?;
        match record.feature_type() {
            "gene" => {
                // first, phase the last gene
                phase_last_gene(gene)?;
                gene = Some(Gene::new(
                    record.attributes().get("gene_id").expect("missing gene_id in GTF"),
                    record.seqname(),
                    Interval::new(*record.start() as u32 - 1, *record.end() as u32),
                    record.attributes().get("gene_biotype").expect("missing gene_biotype in GTF")
                ));
            },
            "transcript" => {
                // register new transcript
                gene.as_mut()
                    .expect("no gene record before transcript in GTF").transcripts.push(
                    Transcript::new(record.attributes().get("transcript_id").expect(
                        "missing transcript_id attribute in GTF"
                    ),
                    record.strand().expect(
                    "missing strand information in GTF"
                    ))
                );
            },
            "exon" => {
                // register exon
                gene.as_mut().expect("no gene record before exon in GTF")
                    .transcripts.last_mut()
                    .expect("no transcript record before exon in GTF")
                    .exons.push(Interval::new(
                    *record.start() as u32 - 1,
                    *record.end() as u32
                ));
            },
            "start_codon" => {
                if record.strand() == Some(Strand::Forward) {
                    gene.as_mut().expect("no gene record before start_codon in GTF")
                        .transcripts.last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons.last_mut()
                        .expect("no exon record before start codon in GTF")
                        .start = *record.start() as u32 - 1;
                }
                else {
                    gene.as_mut().expect("no gene record before start_codon in GTF")
                        .transcripts.last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons.last_mut()
                        .expect("no exon record before start codon in GTF")
                        .end = *record.end() as u32;
                }
            },
            "stop_codon" => {
                if record.strand() == Some(Strand::Forward) {
                    gene.as_mut().expect("no gene record before stop_codon in GTF")
                        .transcripts.last_mut()
                        .expect("no transcript record before stop codon in GTF")
                        .exons.last_mut()
                        .expect("no exon record before stop codon in GTF")
                        .end = *record.end() as u32;
                }
                else {
                    eprintln!("stop_codon_start {}",*record.start() as u32 -1);
                    gene.as_mut().expect("no gene record before stop_codon in GTF")
                        .transcripts.last_mut()
                        .expect("no transcript record before stop codon in GTF")
                        .exons.last_mut()
                        .expect("no exon record before stop codon in GTF")
                        .start = *record.start() as u32 -1;
                }
            }
            _ => continue
        }
    }
    phase_last_gene(gene)?;

    Ok(())
}
