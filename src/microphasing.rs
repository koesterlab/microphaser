use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;
use std::fs;

use csv;
use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;
use bio::utils::Strand;
use rust_htslib::{bam, bcf};
use rust_htslib::bam::record::Cigar;


use common::{Gene, Variant, Interval, Transcript, PhasingStrand};


pub fn bitvector_is_set(b: u64, k: usize) -> bool {
    (b & (1 << k)) != 0
}

pub fn switch_ascii_case(c: u8) -> u8 {
    if c.is_ascii_uppercase() {
        c.to_ascii_lowercase()}
    else {
        c.to_ascii_uppercase()}
}

pub fn switch_ascii_case_vec(v: &Vec<u8>) -> Vec<u8> {
    let c = v[0];
    if c.is_ascii_uppercase() {
        v.to_ascii_lowercase()}
    else {
        v.to_ascii_uppercase()}
}

pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        &Variant::SNV { pos, alt, .. } => {
            debug!("Variant pos: {}", variant.pos());
            debug!("Read pos: {}", read.pos());
            debug!("Read end: {}", read.seq().len() + (read.pos() as usize));
            debug!("Read to check support: {}", String::from_utf8_lossy(read.qname()));
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

#[derive(Debug,Serialize)]
pub struct IDRecord{
    id: String,
    transcript: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    offset: u32,
    freq: f64,
    nvar: u32,
    nsomatic: u32,
    nvariant_sites: u32,
    strand: String,
    positions: String,
    aa_change: String
}

impl IDRecord{
    pub fn update(&self, rec: &IDRecord, offset: u32, seq: Vec<u8>) -> Self {
        let mut shaid = sha1::Sha1::new();
        // generate unique haplotype ID containing position, transcript and sequence
        let id = format!("{:?}{}{}", &seq, &self.transcript, offset);
        shaid.update(id.as_bytes());
        let fasta_id = format!("{}{}", &shaid.digest().to_string()[..15], self.strand.chars().next().unwrap());
        let mut positions = self.positions.clone();
        positions.push_str(&rec.positions);
        let mut aa_change = self.aa_change.clone();
        aa_change.push_str(&rec.aa_change);
        debug!("nvars {} {}",self.nvar,rec.nvar);
        IDRecord{id: fasta_id, transcript: self.transcript.to_owned(), gene_id: self.gene_id.to_owned(), gene_name: self.gene_name.to_owned(), chrom: self.chrom.to_owned(),
            offset: offset + self.offset, freq: self.freq * rec.freq, nvar: self.nvar + rec.nvar, nsomatic: self.nsomatic + rec.nsomatic, nvariant_sites: self.nvariant_sites + rec.nvariant_sites,
            strand: self.strand.to_owned(), positions: positions, aa_change: aa_change}
    }
}


#[derive(Debug)]
pub struct HaplotypeSeq {
    sequence: Vec<u8>,
    record: IDRecord
}

//impl HaplotypeSeq {
//    pub fn new() ->
//    pub fn sequence(&self) -> Vec<u8> {
//        &self.sequence
//    }
//}



#[derive(Debug)]
pub struct Observation{
    read: bam::Record,
    haplotype: u64
}


impl Observation {
    pub fn update_haplotype(&mut self, i: usize, variant: &Variant) -> Result<(), Box<Error>> {
        debug!("Read pos {} ; variant pos {}", self.read.pos() as u32, variant.pos());
        if (self.read.pos() as u32) > variant.pos() {
            panic!("bug: read starts right of variant");
        }
        if supports_variant(&self.read, &variant)? {
            debug!(
                "Read {} supports the variant at {}",
                String::from_utf8_lossy(self.read.qname()), variant.pos()
            );
            self.haplotype |= 1 << i;
            debug!("Haplotype: {}", self.haplotype)
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


    /// Drain variants that are no longer in the window
    pub fn shrink_left(&mut self, k: usize) {
        debug!("self.variants length: {}", self.variants.len());
        debug!("range to drain: 0 - {}",k);
        self.variants.drain(..k);
        debug!("drained!");
        let mask = 2u64.pow(self.ncols()) - 1;
        for obs in self.observations.values_mut().flatten() {
            obs.haplotype = obs.haplotype & mask;
        }
    }

    /// Add new variants
    pub fn extend_right(
        &mut self, new_variants: Vec<Variant>
    ) -> Result<(), Box<Error>> {
        let k = new_variants.len();
        debug!("Extend variants!");
        debug!("New variants {}", k);
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
        debug!("Number of reads(before removal): {}", self.observations.len());
        let observations = self.observations.split_off(&interval_end);
        if !reverse {
            self.observations = observations;//self.observations.split_off(&interval_end);
        }
        debug!("Number of reads(after removal): {}",  self.observations.len());
    }

    /// Check if read has already been added to observations - deprecated
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
        debug!("Read Start: {}, Read End: {} - Window Start: {}, Window End {}", start_pos, end_pos, interval_start, interval_end);
        if end_pos >= interval_end && start_pos <= interval_start {
            // only insert if end_pos is larger than the interval end
            let mut obs = Observation { read: read, haplotype: 0 };
            for (i, variant) in self.variants.iter().enumerate() {
                obs.update_haplotype(i, variant)?;
            }
            // debug!("Read: {}; haplotype: {}", String::from_utf8_lossy(obs.read.qname()), obs.haplotype);
            let pos = match reverse {
                true => start_pos,
                false => end_pos
            };
            self.observations.entry(pos).or_insert_with(|| Vec::new()).push(obs);
//            if reverse {
//                self.observations.entry(start_pos).or_insert_with(|| Vec::new()).push(obs);
//            }
//            else {
//                self.observations.entry(end_pos).or_insert_with(|| Vec::new()).push(obs);
//            }
        }
        Ok(())
    }

    pub fn ncols(&self) -> u32 {
        self.variants.len() as u32
    }

    pub fn nrows(&self) -> usize {
        self.observations.values().map(|o| o.len()).sum()
    }

    pub fn print_haplotypes<O: io::Write> (
        &self,
        gene: &Gene,
        transcript: &Transcript,
        offset: u32,
        exon_end: u32,
        exon_start: u32,
        window_len: u32,
        refseq: &[u8],
        fasta_writer: &mut fasta::Writer<O>,
        tsv_writer: &mut csv::Writer<fs::File>,
        prot_writer: &mut fasta::Writer<fs::File>,
        only_relevant: bool
    ) -> Result<(Vec<HaplotypeSeq>), Box<Error>> {
        let variants = self.variants.iter().collect_vec();
        // count haplotypes
        let mut haplotypes: VecMap<usize> = VecMap::new();
        for obs in self.observations.values().flatten() {
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let mut seq = Vec::with_capacity(window_len as usize);
        debug!("Gene Start {}, Gene End {}", gene.start(), gene.end());
        debug!("Printing at offset: {}",offset);
        debug!("refseq length {}", refseq.len());

        // Strand orientation
        let strand = match transcript.strand {
            PhasingStrand::Reverse => "Reverse",
            PhasingStrand::Forward => "Forward"
        };

        let mut haplotypes_vec = Vec::new();

        for (haplotype, count) in haplotypes.iter() {
            // VecMap forces usize as type for keys, but our haplotypes as u64
            let haplotype = haplotype as u64;
            debug!("Offset: {}", offset);
            debug!("Haplotype: {} ; count: {}", haplotype, count);
            debug!("Variants len: {}", variants.len());
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
                    debug!("window_end: {}", window_end);
                    debug!("i: {}", i);
                    debug!("j: {}", j);
                    // TODO what happens if a deletion starts upstream of window and overlaps it
                    
                    while j < variants.len() && i == variants[j].pos() {
                        debug!("j: {}, variantpos: {}", j, variants[j].pos());
                        if bitvector_is_set(haplotype, j) {
                            debug!("Haplotype: {} ; j: {}", haplotype, j);
                            if (j + 1) < variants.len() && i == variants[j+1].pos() { 
                                j += 1;
                            }
                            match variants[j] {
                                // if SNV, we push the alternative base instead of the normal one, and change the case of the letter for visualisation
                                &Variant::SNV { alt, .. } => {
                                    seq.push(switch_ascii_case(alt));
                                    i += 1;
                                },
                                // if insertion, we insert the new bases (with changed case) and decrease the window-end, since we added bases and made the sequence longer
                                &Variant::Insertion { seq: ref s, .. } => {
                                    seq.extend(switch_ascii_case_vec(s).into_iter());
                                    i += 1;
                                    window_end -= (s.len() as u32) -1;
                                },
                                // if deletion, we push the remaining base and increase the index to jump over the deleted bases. Then, we increase the window-end since we lost bases and need to fill up to 27.
                                &Variant::Deletion { len, .. } => {
                                    seq.push(refseq[(i - gene.start()) as usize]);
                                    i += len + 1;
                                    window_end += len + 1;
                                }
                            }

                            // counting somatic variants
                            if !variants[j].is_germline() {
                                n_somatic += 1;
                            }
                            // counting total number of variants
                            n_variants += 1;
                            j += 1;
                        } else {
                            j += 1;
                        }
                    }
                    // if no variant, just push the reference sequence
                    seq.push(refseq[(i - gene.start()) as usize]);
                    i += 1
                }
            }



            let mut shaid = sha1::Sha1::new();
            // generate unique haplotype ID containing position, transcript and sequence
            let mut id = format!("{:?}{}{}", &seq, transcript.id, offset);
            shaid.update(id.as_bytes());
            let fasta_id = format!("{}{}", &shaid.digest().to_string()[..15], strand.chars().next().unwrap());
            // gathering meta information on haplotype
            let mut n_variantsites = 0;
            let mut c = 0;
            let mut p_changes = String::new();
            let mut var_pos = String::new();
            while c<variants.len() as u32 {
                if c == 0 {
                    var_pos = format!("{}",variants[0].pos());
                    p_changes = variants[0].prot_change();
                    n_variantsites += 1;
                }
                else {
                    var_pos = format!("{}|{}", var_pos, variants[c as usize].pos());
                    p_changes = format!("{}|{}",p_changes,variants[c as usize].prot_change());
                    if !(variants[c as usize].pos() == variants[(c - 1) as usize].pos()) {
                        n_variantsites += 1;
                    }
                }
                c += 1
            }
            let record = IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(), gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(), chrom: gene.chrom.to_owned(), offset: offset, freq: freq, nvar: n_variants, nsomatic: n_somatic, nvariant_sites: n_variantsites as u32, strand: strand.to_string(), positions: var_pos.to_owned(), aa_change: p_changes.to_owned()};
            
            // make haplotype record to carry over to next exon
            let rest = exon_end - window_end;
            debug!("Rest: {}", rest);
            let start = offset - exon_start;
            debug!("Start: {}", start);
            // if there are split-codon bases left at the end of the exon, add them to the sequence
            let mut hap_seq = HaplotypeSeq {sequence: Vec::new(), record: IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(), gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(), chrom: gene.chrom.to_owned(), offset: offset, freq: freq, nvar: n_variants, nsomatic: n_somatic, nvariant_sites: n_variantsites as u32, strand: strand.to_string(), positions: var_pos.to_owned(), aa_change: p_changes.to_owned()}};
            if rest < 3 {
                let mut newseq = Vec::new();
                newseq.extend(&seq[3 as usize..window_len as usize]);
                newseq.extend_from_slice(&refseq[(window_end - gene.start()) as usize..(window_end + rest - gene.start()) as usize]);
                hap_seq = HaplotypeSeq {sequence: newseq, record: IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(), gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(), chrom: gene.chrom.to_owned(), offset: offset, freq: freq, nvar: n_variants, nsomatic: n_somatic, nvariant_sites: n_variantsites as u32, strand: strand.to_string(), positions: var_pos.to_owned(), aa_change: p_changes.to_owned()}};
            }
            // if there are split-codon bases at the start of an exon, add them to the sequence
            if start < 3 {
                let mut newseq = refseq[(offset - start - gene.start()) as usize..(offset - gene.start()) as usize].to_vec();
                debug!("Starting_sequence: {:?}", String::from_utf8_lossy(&seq));
                newseq.extend(&seq[..(window_len - 3) as usize]);
                hap_seq = HaplotypeSeq {sequence: newseq, record: IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(), gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(), chrom: gene.chrom.to_owned(), offset: offset, freq: freq, nvar: n_variants, nsomatic: n_somatic, nvariant_sites: n_variantsites as u32, strand: strand.to_string(), positions: var_pos.to_owned(), aa_change: p_changes.to_owned()}};
            }

            haplotypes_vec.push(hap_seq);

            debug!("relevant_check: {}, nvar: {}, freq: {} ", !(only_relevant), record.nvar > 0, record.freq < 1.00);
            debug!("is_relevant: {}", !(only_relevant) ||  record.freq < 1.00 || record.nvar > 0);
            // collect "background" haplotypes (only wildtype) for similarity comparison with neopeptides
            if !(only_relevant) && record.nvar == 0 && record.freq == 1.00 {
                prot_writer.write(&format!("{}", record.id), None, &seq[..window_len as usize])?;
                }
            // filter relevant haplotypes : variant haplotypes and their wildtype counterparts
            if record.freq < 1.00 || record.nvar > 0 {
                fasta_writer.write(&format!("{}", record.id), None, &seq[..window_len as usize])?;
                tsv_writer.serialize(record)?;
            }



        }
        Ok(haplotypes_vec)
    }
}


//pub fn close_splicing_gap<O: io::Write> (
//    prev_hap_vec: &Vec<HaplotypeSeq>,
//    hap_vec: &Vec<HaplotypeSeq>,
//    fasta_writer: &mut fasta::Writer<O>,
//    tsv_writer: &mut csv::Writer<fs::File>,
//    prot_writer: &mut fasta::Writer<fs::File>
//) -> Result<(), Box<Error>> {
////    for prev_hapseq in prev_hap_vec:
////        
//    Ok(())
//}


pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    gene: &Gene,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::buffer::RecordBuffer,
    fasta_writer: &mut fasta::Writer<O>,
    tsv_writer: &mut csv::Writer<fs::File>,
    prot_writer: &mut fasta::Writer<fs::File>,
    window_len: u32,
    refseq: &mut Vec<u8>,
    only_relevant: bool
) -> Result<(), Box<Error>> {
    // if an exon is near to the gene end, a deletion could cause refseq to overflow, so we increase the length of refseq
    let end_overflow = 100;
    fasta_reader.read(&gene.chrom, gene.start() as u64, (gene.end() + end_overflow) as u64, refseq)?;
    let mut variant_tree = BTreeMap::new();
    let mut read_tree = BTreeMap::new();
    debug!("Start Phasing");
    read_buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;

    // load read buffer into BTree
    for rec in read_buffer.iter() {
        read_tree.entry(rec.pos() as u32).or_insert_with(|| Vec::new()).push(rec.clone())
    }
    debug!("Reads Tree length: {}", read_tree.len());

    let max_read_len = 101;

    // load variant buffer into BTree
    let (_addedvars, _deletedvars) = variant_buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;
    let _vars = variant_buffer.iter_mut().map(|rec| variant_tree.insert(rec.pos(), Variant::new(rec).unwrap())).collect_vec();


    for transcript in &gene.transcripts {
        debug!("is coding: {}", transcript.is_coding());
        if !(transcript.is_coding()) {
            continue
        }
        debug!("Transcript strand orientation: {:?}", transcript.strand);
        let mut observations = ObservationMatrix::new();
        let mut frameshifts = BTreeMap::new();
        frameshifts.insert(0, 0);
        // Possible rest of an exon that does not form a complete codon yet
        let mut exon_rest = 0;
        // Haplotypes from the end of the last exon
        let mut prev_hap_vec: Vec<HaplotypeSeq> = Vec::new();
        let mut hap_vec: Vec<HaplotypeSeq> = Vec::new();
        // Possible variants that are still in the BTree after leaving the last exon (we do not drain variants if we leave an exon)
        let mut last_window_vars = 0;
        for exon in &transcript.exons {
        debug!("Exon Start: {}", exon.start);
        debug!("Exon End: {}", exon.end);
            // Possible offset at the exon start, first nucleotides could be part of a codon started in the previous exon
            let current_exon_offset = match exon_rest {
                0 => 0,
                _ => 3 - exon_rest
            };
            exon_rest = 0;
            let mut offset = if transcript.strand == PhasingStrand::Reverse {
                exon.end - window_len - current_exon_offset
            } else {
                exon.start + current_exon_offset
            };
            let mut old_offset = offset;
            debug!("Variants left from previous Exon: {}", last_window_vars);
            observations.shrink_left(last_window_vars);
            last_window_vars = 0;
            loop {
                let valid = match transcript.strand {
                    PhasingStrand::Reverse => offset >= exon.start,
                    PhasingStrand::Forward => offset + window_len < exon.end
                };
                if !valid {
                    break;
                }
                debug!("Offset {}, old offset {}", offset, old_offset);
                // advance window to next position
                let nvars = variant_tree.range(offset..(offset + window_len)).map(|var| var.1).flatten().count();
                // store number of variants in window in case it is the last window for this exon
                last_window_vars = nvars;
                debug!("Variants in window: {}",nvars);
                // first window in the exon, all variants found are newly added
                let mut added_vars = if offset == old_offset {
                    nvars
                // if we advance the window (forward or reverse), just the newly added variants are counted
                // forward orientation
                } else if offset > old_offset {
                    variant_tree.range((old_offset + window_len)..(offset + window_len)).map(|var| var.1).flatten().count()
                // reverse orientation
                } else {
                    variant_tree.range(offset..old_offset).map(|var| var.1).flatten().count()
                };

                // first window in the exon, no variants are deleted
                let mut deleted_vars = if offset == old_offset {
                    0
                // if we advance the window (forward or reverse), we will delete all variants that drop out of the window bounds
                // forward orientation
                } else if offset > old_offset {
                    variant_tree.range(old_offset..offset).map(|var| var.1).flatten().count()
                // reverse orientation
                } else {
                    variant_tree.range((offset + window_len)..(old_offset + window_len)).map(|var| var.1).flatten().count()
                };


                debug!("Offset: {} - max_read_len - window_len {}", offset, (max_read_len - window_len));
                let mut reads = if transcript.strand == PhasingStrand::Reverse {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.end - window_len - current_exon_offset {
                        read_tree.range((offset - (max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                    // while advancing the window (reverse orientation), we only add reads that end in the range between old and new window end, so we don't count any read twice 
                    else {
                        read_tree.range((offset - (max_read_len - window_len))..(offset - (max_read_len - window_len) + 1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                }
                else {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.start + current_exon_offset {
                        read_tree.range((offset-(max_read_len - window_len))..(offset+1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                    // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
                    else {
                        read_tree.range((offset-1)..(offset+1)).map(|rec| rec.1).flatten().collect_vec()
                    }
                };



                debug!("Variants added: {}", added_vars);
                debug!("Variants deleted: {}", deleted_vars);
                {
                    // delete rows
                    let reverse = match transcript.strand {
                        PhasingStrand::Reverse => true,
                        PhasingStrand::Forward => false};
                    if reverse {
                        observations.cleanup_reads(offset, reverse);
                    }
                    else {
                        observations.cleanup_reads(offset + window_len, reverse);
                    }
                    // delete columns
                    observations.shrink_left(deleted_vars);

                    // add new reads
                    debug!("Reads: {}", reads.len());
                    for read in reads {
                        observations.push_read(read.clone(), offset + window_len, offset, reverse)?;
//                        if transcript.strand == Strand::Reverse {
//                            observations.push_read(read.clone(), offset + window_len, offset, true)?;
//                        }
//                        else {
//                            observations.push_read(read.clone(), offset + window_len, offset, false)?;
//                        }
                    }


                    // collect variants
                    let variants = match transcript.strand {
                        PhasingStrand::Reverse => variant_tree.range_mut(offset..(offset + window_len)).rev()
                        .map(|var| var.1.clone()).flatten().skip(nvars - added_vars).collect_vec(),
                        PhasingStrand::Forward => variant_tree.range_mut(offset..(offset + window_len)
                        ).map(|var| var.1.clone()).flatten().skip(nvars - added_vars).collect_vec()
                    };
                    debug!("Variants(after deleting and adding): {}", variants.len());
                    // determine frameshifts
                    for variant in &variants {
                        debug!("Variants!");
                        debug!("Variant Pos: {}", variant.pos());
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
                        // possible shift if exon starts with the rest of a split codon (splicing)
                        let coding_shift = match transcript.strand {
                                PhasingStrand::Forward => offset - exon.start,
                                PhasingStrand::Reverse => exon.end - offset
                        };
                        debug!("Offset - Exonstart % 3: {}", coding_shift % 3);
                        debug!("Shift: {}", frameshift + current_exon_offset);
                        if coding_shift % 3 == frameshift + current_exon_offset {
                            // print haplotypes
                            debug!("Should print haplotypes");
//                            hap_vec = observations.print_haplotypes(
//                                gene, transcript, offset, exon.end, exon.start, window_len, refseq, fasta_writer, tsv_writer, prot_writer, only_relevant
//                            ).unwrap();
                            // possible unfinished codon at the end of an exon that continues at the start of the next exon
                            exon_rest = match transcript.strand {
                                PhasingStrand::Forward => exon.end - (offset + window_len),
                                PhasingStrand::Reverse => offset - exon.start
                            };
                            debug!("Exon Rest {}", exon_rest);
                            if exon_rest < 3 {
                                prev_hap_vec = observations.print_haplotypes(
                                gene, transcript, offset, exon.end, exon.start, window_len, refseq, fasta_writer, tsv_writer, prot_writer, only_relevant
                            ).unwrap(); }
                            else {
                                hap_vec = observations.print_haplotypes(
                                gene, transcript, offset, exon.end, exon.start, window_len, refseq, fasta_writer, tsv_writer, prot_writer, only_relevant
                            ).unwrap(); }
                        }
                    }

                    let at_splice_side = match transcript.strand {
                                PhasingStrand::Forward => offset - current_exon_offset == exon.start,
                                PhasingStrand::Reverse => offset + window_len + current_exon_offset == exon.end 
                    };


                    if at_splice_side {
                        let first_hap_vec = match transcript.strand {
                                PhasingStrand::Forward => &hap_vec,
                                PhasingStrand::Reverse => &prev_hap_vec 
                        };
                        let sec_hap_vec = match transcript.strand {
                                PhasingStrand::Forward => &prev_hap_vec,
                                PhasingStrand::Reverse => &hap_vec
                        };
                        //close_splicing_gap(&prev_hap_vec, &hap_vec, fasta_writer, tsv_writer, prot_writer);
                        //let it = prev_hap_vec.cartesian_product(hap_vec);
                        //let mut splice_sequences = Vec::new();
                        for hapseq in first_hap_vec{
                            let sequence = &hapseq.sequence;
                            let record = &hapseq.record;
                            for prev_hapseq in sec_hap_vec {
                                let mut prev_sequence = prev_hapseq.sequence.clone();
                                let prev_record = &prev_hapseq.record;
//                                match transcript.strand {
//                                    PhasingStrand::Forward => prev_sequence.extend(sequence),
//                                    PhasingStrand::Reverse => sequence.extend(prev_sequence)
//                                };
//                                let complete_seq = match transcript.strand {
//                                    PhasingStrand::Forward => prev_sequence,
//                                    PhasingStrand::Reverse => sequence
//                                };
                                prev_sequence.extend(sequence);
                                debug!("Complete Sequence : {:?}", String::from_utf8_lossy(&prev_sequence));
                                let mut splice_offset = 0;
                                while (splice_offset + window_len < prev_sequence.len() as u32) {
                                    let out_seq = &prev_sequence[splice_offset as usize..(splice_offset + window_len) as usize];
                                    debug!("Out Sequence : {:?}", String::from_utf8_lossy(&out_seq));
                                    let out_record = prev_record.update(record, splice_offset + 3, out_seq.to_vec());
                                    if !(only_relevant) && out_record.nvar == 0 && record.freq == 1.00 {
                                        prot_writer.write(&format!("{}", out_record.id), None, &out_seq[..window_len as usize])?;
                                    }
                                    // filter relevant haplotypes : variant haplotypes and their wildtype counterparts
                                    if record.freq < 1.00 || out_record.nvar > 0 {
                                        fasta_writer.write(&format!("{}", out_record.id), None, &out_seq[..window_len as usize])?;
                                        tsv_writer.serialize(out_record)?;
                                    }
                                    splice_offset += 3;
                                }
                                //splice_sequences.push(HaplotypeSeq {sequence: prev_sequence, frequency: complete_freq});
                            }
                        }
                        //prev_hap_vec = Vec::new();
                    }
                    //prev_hap_vec = hap_vec;

                    old_offset = offset;
                    match transcript.strand {
                        PhasingStrand::Reverse => offset -= 1,
                        PhasingStrand::Forward => offset += 1
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
    tsv_writer: &mut csv::Writer<fs::File>,
    prot_writer: &mut fasta::Writer<fs::File>,
    window_len: u32,
    only_relevant: bool
) -> Result<(), Box<Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader);
    let mut variant_buffer = bcf::buffer::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence
    debug!("refseq length {}", refseq.len());



    let mut gene = None;
    let mut phase_last_gene = | gene: Option<Gene>| -> Result<(), Box<Error>> {
        if let Some(ref gene) = gene {
            if gene.biotype == "protein_coding" {
                phase_gene(
                    &gene, fasta_reader, &mut read_buffer,
                    &mut variant_buffer, fasta_writer, tsv_writer, prot_writer,
                    window_len,
                    &mut refseq,
                    only_relevant
                )?;
            }
        }
        Ok(())
    };


    for record in gtf_reader.records() {
        let record = record?;
//        match record.strand() {
//            Some(Strand::Forward) => (),
//            Some(Strand::Reverse) => (),
//            _ => panic!("Strand unknown in GTF. Should be + or -")}
        match record.feature_type() {
            "gene" => {
                // first, phase the last gene
                phase_last_gene(gene)?;
                gene = Some(Gene::new(
                    record.attributes().get("gene_id").expect("missing gene_id in GTF"),
                    record.attributes().get("gene_name").expect("missing gene_name in GTF"),
                    record.seqname(),
                    Interval::new(*record.start() as u32 - 1, *record.end() as u32),
                    record.attributes().get("gene_biotype").expect("missing gene_biotype in GTF")
                ));
            },
            "transcript" => {
                // register new transcript
                debug!("Transcript found");
                gene.as_mut()
                    .expect("no gene record before transcript in GTF").transcripts.push(
                    Transcript::new(record.attributes().get("transcript_id").expect(
                        "missing transcript_id attribute in GTF"
                    ),
                    PhasingStrand::from(record.strand().expect(
                    "missing strand information in GTF"
                    )))
                );
            },
            "CDS" => {
                debug!("CDS found");
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
                    debug!("stop_codon_start {}",*record.start() as u32 -1);
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
