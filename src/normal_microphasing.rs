use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;

use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;
use bio::utils::Strand;
use rust_htslib::{bam, bcf};
use rust_htslib::bam::record::Cigar;


use crate::common::{Gene, Variant, Interval, Transcript, PhasingStrand};


pub fn bitvector_is_set(b: u64, k: usize) -> bool {
    (b & (1 << k)) != 0
}

pub fn switch_ascii_case(c: u8, r: u8) -> u8 {
    if r.is_ascii_uppercase() {
        c.to_ascii_lowercase()}
    else {
        c}
}

pub fn switch_ascii_case_vec(v: &Vec<u8>, r: u8) -> Vec<u8> {
    if r.is_ascii_uppercase() {
        v.to_ascii_lowercase()}
    else {
        v.to_ascii_uppercase()}
}

pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        &Variant::SNV { pos, alt, .. } => {
//            debug!("Variant pos: {}", variant.pos());
//            debug!("Read pos: {}", read.pos());
//            debug!("Read end: {}", read.seq().len() + (read.pos() as usize));
//            debug!("Read to check support: {}", String::from_utf8_lossy(read.qname()));
            let b = match read.cigar().read_pos(pos, false, false) {
                Ok(None) => return Ok(false),
                Ok(Some(p)) => read.seq()[p as usize],
                _ => return Ok(false)
            };
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
    nsomvariant_sites: u32,
    strand: String,
    somatic_positions: String,
    somatic_aa_change: String,
    germline_positions: String,
    germline_aa_change: String
}

impl IDRecord{
    pub fn update(&self, rec: &IDRecord, offset: u32, seq: Vec<u8>) -> Self {
        let mut shaid = sha1::Sha1::new();
        // generate unique haplotype ID containing position, transcript and sequence
        let id = format!("{:?}{}{}", &seq, &self.transcript, offset);
        shaid.update(id.as_bytes());
        let fasta_id = format!("{}{}", &shaid.digest().to_string()[..15], self.strand.chars().next().unwrap());
        let mut somatic_positions = self.somatic_positions.clone();
        somatic_positions.push_str(&rec.somatic_positions);
        let mut somatic_aa_change = self.somatic_aa_change.clone();
        somatic_aa_change.push_str(&rec.somatic_aa_change);
        let mut germline_positions = self.germline_positions.clone();
        germline_positions.push_str(&rec.germline_positions);
        let mut germline_aa_change = self.germline_aa_change.clone();
        germline_aa_change.push_str(&rec.germline_aa_change);
        debug!("nvars {} {}",self.nvar,rec.nvar);
        IDRecord{id: fasta_id, transcript: self.transcript.to_owned(), gene_id: self.gene_id.to_owned(), gene_name: self.gene_name.to_owned(), chrom: self.chrom.to_owned(),
            offset: offset + self.offset, freq: self.freq * rec.freq, nvar: self.nvar + rec.nvar, nsomatic: self.nsomatic + rec.nsomatic, nvariant_sites: self.nvariant_sites + rec.nvariant_sites, nsomvariant_sites: self.nsomvariant_sites + rec.nsomvariant_sites,
            strand: self.strand.to_owned(), somatic_positions: somatic_positions, somatic_aa_change: somatic_aa_change, germline_positions: germline_positions, germline_aa_change: germline_aa_change}
    }

    pub fn add_freq(&self, freq: f64) -> Self {
        let new_nvar = match freq{
            f if f > 0.0 => self.nvar - 1,
            _ => self.nvar
        };
        let new_somatic = match new_nvar < self.nsomatic {
            true => self.nsomatic - 1,
            false => self.nsomatic
        };
        IDRecord{id: self.id.to_owned(), transcript: self.transcript.to_owned(), gene_id: self.gene_id.to_owned(), gene_name: self.gene_name.to_owned(), chrom: self.chrom.to_owned(),
            offset: self.offset, freq: self.freq + freq, nvar: new_nvar, nsomatic: new_somatic , nvariant_sites: self.nvariant_sites, nsomvariant_sites: self.nsomvariant_sites,
            strand: self.strand.to_owned(), somatic_positions: self.somatic_positions.to_owned(), somatic_aa_change: self.somatic_aa_change.to_owned(), germline_positions: self.germline_positions.to_owned(), germline_aa_change: self.germline_aa_change.to_owned()}
    }

}


#[derive(Debug)]
pub struct HaplotypeSeq {
    sequence: Vec<u8>,
    record: IDRecord
}


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
        for obs in Itertools::flatten(self.observations.values_mut()) {
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
            for obs in Itertools::flatten(self.observations.values_mut()) {
                obs.haplotype <<= k;
            }
        }
        for obs in Itertools::flatten(self.observations.values_mut()) {
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
        fasta_writer: &mut fasta::Writer<O>
    ) -> Result<(Vec<HaplotypeSeq>), Box<Error>> {

        let variants_forward = self.variants.iter().collect_vec();
        let mut variants_reverse = variants_forward.clone();
        variants_reverse.reverse();
        let variants = match transcript.strand {
            PhasingStrand::Reverse => variants_reverse,
            PhasingStrand::Forward => variants_forward
        };
        // count haplotypes
        let mut haplotypes: VecMap<usize> = VecMap::new();
        for obs in Itertools::flatten(self.observations.values()) {
            debug!("obs {:?}",obs);
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
            debug!("window end: {}", window_end);
            // Profile for all variants: 0 - reference, 1 - germline, 2 - somatic
            let mut variant_profile = Vec::new();
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
                        if freq == 1.0 && !(variants[j].is_germline()) {
                            j += 1;
                            variant_profile.push(0);
                            continue;
                        }
                        if bitvector_is_set(haplotype, j) {
                            debug!("Haplotype: {} ; j: {}", haplotype, j);
                            if (j + 1) < variants.len() && i == variants[j+1].pos() { 
                                j += 1;
                            }
                            match variants[j] {
                                // if SNV, we push the alternative base instead of the normal one, and change the case of the letter for visualisation
                                &Variant::SNV { alt, .. } => {
                                    seq.push(switch_ascii_case(alt, refseq[(i - gene.start()) as usize]));
                                    i += 1;
                                },
                                // if insertion, we insert the new bases (with changed case) and decrease the window-end, since we added bases and made the sequence longer
                                &Variant::Insertion { seq: ref s, .. } => {
                                    seq.extend(switch_ascii_case_vec(s, refseq[(i - gene.start()) as usize]).into_iter());
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
                                variant_profile.push(2);
                            }
                            else {
                                variant_profile.push(1);
                            }
                            // counting total number of variants
                            n_variants += 1;
                            j += 1;
                        } else {
                            variant_profile.push(0);
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
            let id = format!("{:?}{}{}", &seq, transcript.id, offset);
            shaid.update(id.as_bytes());
            let fasta_id = format!("{}{}", &shaid.digest().to_string()[..15], strand.chars().next().unwrap());
            // gathering meta information on haplotype
            let mut n_variantsites = 0;
            let mut n_som_variantsites = 0;
            let mut c = 0;
            // protein changes
            let mut somatic_p_changes_vec = Vec::new();
            let mut germline_p_changes_vec = Vec::new();
            // position of the variant
            let mut somatic_var_pos_vec = Vec::new();
            let mut germline_var_pos_vec = Vec::new();
            // gather information iterating over the variants
            debug!("Variant profile len: {}", variant_profile.len());
            while c < variants.len() as u32 {
                if c == 0 {
                    if c < variant_profile.len() as u32 {
                        match variant_profile[0] {
                            // somatic
                            2 => {
                                somatic_var_pos_vec.push(variants[c as usize].pos().to_string());
                                somatic_p_changes_vec.push(variants[c as usize].prot_change());
                            },
                            // germline
                            1 => {
                                germline_var_pos_vec.push(variants[c as usize].pos().to_string());
                                germline_p_changes_vec.push(variants[c as usize].prot_change());
                            },
                            // not present in this haplotype
                            _ => {}
                        }
                    }
                    n_variantsites += 1;
                    if !(variants[0].is_germline()) {
                        n_som_variantsites += 1;
                    }
                }
                else {
                    if c < variant_profile.len() as u32 {
                        match variant_profile[c as usize] {
                            // somatic
                            2 => {
                                somatic_var_pos_vec.push(variants[c as usize].pos().to_string());
                                somatic_p_changes_vec.push(variants[c as usize].prot_change());
                            },
                            // germline
                            1 => {
                                germline_var_pos_vec.push(variants[c as usize].pos().to_string());
                                germline_p_changes_vec.push(variants[c as usize].prot_change());
                            },
                            // not present in this haplotype
                            _ => {}
                        }
                    }
                    if !(variants[c as usize].pos() == variants[(c - 1) as usize].pos()) {
                        n_variantsites += 1;
                        if !(variants[c as usize].is_germline()) {
                            n_som_variantsites += 1;
                        }
                    }
                }
                c += 1
            }

            let somatic_var_pos = somatic_var_pos_vec.join("|");
            let somatic_p_changes = somatic_p_changes_vec.join("|");
            let germline_var_pos = germline_var_pos_vec.join("|");
            let germline_p_changes = germline_p_changes_vec.join("|");

            // build the info record
            let record = IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(),
                gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(),
                chrom: gene.chrom.to_owned(), offset, freq, nvar: n_variants, nsomatic: n_somatic,
                nvariant_sites: n_variantsites as u32, nsomvariant_sites: n_som_variantsites as u32,
                strand: strand.to_string(), somatic_positions: somatic_var_pos.to_owned(), somatic_aa_change: somatic_p_changes.to_owned(),
                germline_positions: germline_var_pos.to_owned(), germline_aa_change: germline_p_changes.to_owned()};
            
            // make haplotype record to carry over to next exon
            debug!("exon end: {}, window end {}", exon_end, window_end);
            let rest = match window_end > exon_end {
                true => 0,
                false => exon_end - window_end
            };
            debug!("Rest: {}", rest);
            let start = offset - exon_start;
            debug!("Start: {}", start);

            let newseq = match rest < 3 {
                // if there are split-codon bases left at the end of the exon, add them to the sequence
                true => {
                    let mut s = Vec::new();
                    s.extend(&seq[3 as usize..window_len as usize]);
                    s.extend_from_slice(&refseq[(window_end - gene.start()) as usize..(window_end + rest - gene.start()) as usize]);
                    s
                }
                false => match start < 3 {
                     // if there are split-codon bases at the start of an exon, add them to the sequence
                    true => {
                        let mut s = refseq[(offset - start - gene.start()) as usize..(offset - gene.start()) as usize].to_vec();
                        s.extend(&seq[..(window_len - 3) as usize]);
                        s
                    }
                    false => Vec::new()
                }
            };

            let hap_seq = HaplotypeSeq {sequence: newseq,
                record: IDRecord {id: fasta_id.to_owned(), transcript: transcript.id.to_owned(),
                    gene_id: gene.id.to_owned(), gene_name: gene.name.to_owned(), chrom: gene.chrom.to_owned(),
                    offset, freq, nvar: n_variants, nsomatic: n_somatic,
                    nvariant_sites: n_variantsites as u32, nsomvariant_sites: n_som_variantsites as u32,
                    strand: strand.to_string(), somatic_positions: somatic_var_pos.to_owned(), somatic_aa_change: somatic_p_changes.to_owned(),
                    germline_positions: germline_var_pos.to_owned(), germline_aa_change: germline_p_changes.to_owned()}};

            haplotypes_vec.push(hap_seq);
            // print haplotypes
            fasta_writer.write(&format!("{}", record.id), None, &seq[..window_len as usize])?;
//            tsv_writer.serialize(record)?;
        }
        Ok(haplotypes_vec)
    }
}



pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    gene: &Gene,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::buffer::RecordBuffer,
    fasta_writer: &mut fasta::Writer<O>,
    window_len: u32,
    refseq: &mut Vec<u8>
) -> Result<(), Box<Error>> {
    // if an exon is near to the gene end, a deletion could cause refseq to overflow, so we increase the length of refseq
    let end_overflow = 100;
    fasta_reader.read(&gene.chrom, gene.start() as u64, (gene.end() + end_overflow) as u64, refseq)?;
    let mut variant_tree = BTreeMap::new();
    let mut read_tree = BTreeMap::new();
    debug!("Start Phasing");
    read_buffer.fetch(&gene.chrom.as_bytes(),gene.start(),gene.end())?;

    let mut max_read_len = 50 as u32;
    // load read buffer into BTree
    for rec in read_buffer.iter() {
        if rec.seq().len() as u32 > max_read_len {
        max_read_len = rec.seq().len() as u32;}
        read_tree.entry(rec.pos() as u32).or_insert_with(|| Vec::new()).push(rec.clone())
    }
    debug!("Reads Tree length: {}", read_tree.len());

    //let max_read_len = 101;

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
                    PhasingStrand::Forward => offset + window_len <= exon.end
                };
                if !valid {
                    break;
                }
                debug!("Offset {}, old offset {}", offset, old_offset);
                // advance window to next position
                let nvars = Itertools::flatten(variant_tree.range(offset..(offset + window_len)).map(|var| var.1)).count();
                // store number of variants in window in case it is the last window for this exon
                last_window_vars = nvars;
                debug!("Variants in window: {}",nvars);
                // first window in the exon, all variants found are newly added
                let added_vars = if offset == old_offset {
                    nvars
                // if we advance the window (forward or reverse), just the newly added variants are counted
                // forward orientation
                } else if offset > old_offset {
                    Itertools::flatten(variant_tree.range((old_offset + window_len)..(offset + window_len)).map(|var| var.1)).count()
                // reverse orientation
                } else {
                    Itertools::flatten(variant_tree.range(offset..old_offset).map(|var| var.1)).count()
                };

                // first window in the exon, no variants are deleted
                let deleted_vars = if offset == old_offset {
                    0
                // if we advance the window (forward or reverse), we will delete all variants that drop out of the window bounds
                // forward orientation
                } else if offset > old_offset {
                    Itertools::flatten(variant_tree.range(old_offset..offset).map(|var| var.1)).count()
                // reverse orientation
                } else {
                    Itertools::flatten(variant_tree.range((offset + window_len)..(old_offset + window_len)).map(|var| var.1)).count()
                };


                debug!("Offset: {} - max_read_len - window_len {}", offset, (max_read_len - window_len));
                let reads = if transcript.strand == PhasingStrand::Reverse {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.end - window_len - current_exon_offset {
                        Itertools::flatten(read_tree.range((offset - (max_read_len - window_len))..(offset+1)).map(|rec| rec.1)).collect_vec()
                    }
                    // while advancing the window (reverse orientation), we only add reads that end in the range between old and new window end, so we don't count any read twice 
                    else {
                        Itertools::flatten(read_tree.range((offset - (max_read_len - window_len))..(offset - (max_read_len - window_len) + 1)).map(|rec| rec.1)).collect_vec()
                    }
                }
                else {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.start + current_exon_offset {
                        Itertools::flatten(read_tree.range((offset-(max_read_len - window_len))..(offset+1)).map(|rec| rec.1)).collect_vec()
                    }
                    // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
                    else {
                        Itertools::flatten(read_tree.range((offset-1)..(offset+1)).map(|rec| rec.1)).collect_vec()
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
                        PhasingStrand::Reverse => Itertools::flatten(variant_tree.range_mut(offset..(offset + window_len)).rev()
                        .map(|var| var.1.clone())).skip(nvars - added_vars).collect_vec(),
                        PhasingStrand::Forward => Itertools::flatten(variant_tree.range_mut(offset..(offset + window_len)
                        ).map(|var| var.1.clone())).skip(nvars - added_vars).collect_vec()
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
                                gene, transcript, offset, exon.end, exon.start, window_len, refseq, fasta_writer
                            ).unwrap(); }
                            else {
                                hap_vec = observations.print_haplotypes(
                                gene, transcript, offset, exon.end, exon.start, window_len, refseq, fasta_writer
                            ).unwrap(); }
                        }
                    }

                    // check if the current offset is at a splice side
                    let at_splice_side = match transcript.strand {
                                PhasingStrand::Forward => offset - current_exon_offset == exon.start,
                                PhasingStrand::Reverse => offset + window_len + current_exon_offset == exon.end 
                    };

                    // at a splice side, merge the last sequence of the prev exon and the first sequence of the next exon
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
                        let mut output_map: BTreeMap<(u32, Vec<u8>), (Vec<u8>, IDRecord)> = BTreeMap::new();

                        // iterate over all combinations of splice side haplotypes
                        for hapseq in first_hap_vec{
                            let sequence = &hapseq.sequence;
                            let record = &hapseq.record;
                            for prev_hapseq in sec_hap_vec {
                                let mut prev_sequence = prev_hapseq.sequence.clone();
                                let prev_record = &prev_hapseq.record;
                                // paste together a sequence spanning the splice side
                                prev_sequence.extend(sequence);
                                debug!("Complete Sequence : {:?}", String::from_utf8_lossy(&prev_sequence));
                                // slide window over the spanning sequence
                                let mut splice_offset = 0;
                                while splice_offset + window_len <= prev_sequence.len() as u32 {
                                    let out_seq = &prev_sequence[splice_offset as usize..(splice_offset + window_len) as usize];
                                    let out_record = prev_record.update(record, splice_offset + 3, out_seq.to_vec());
                                    debug!("Out Sequence : {:?}", String::from_utf8_lossy(&out_seq));
                                    // check if the sequence was already printed at that position, i.e. the variant defining the different haplotypes left the window
                                    let id_tuple = (splice_offset, out_seq.to_vec());
                                    let old_freq = match output_map.get_mut(&id_tuple) {
                                        Some(x) => x.1.freq,
                                        None => 0.0
                                    };
                                    *output_map.entry(id_tuple).or_insert((out_seq.to_vec(), out_record)) = (out_seq.to_vec(), out_record.add_freq(old_freq));
                                    splice_offset += 3;
                                }
                            }
                        }
//                      let out_record = prev_record.update(record, splice_offset + 3, out_seq.to_vec());
                        for (_key, val) in output_map.iter() {
                            let out_record = &val.1;
                            let out_seq = &val.0;
                            fasta_writer.write(&format!("{}", out_record.id), None, &out_seq[..window_len as usize])?;
//                            tsv_writer.serialize(out_record)?;
                        }
                    }
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
    window_len: u32,
) -> Result<(), Box<Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader);
    let mut variant_buffer = bcf::buffer::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence
    debug!("refseq length {}", refseq.len());
    debug!("Stared Phasing");
    let mut gene = None;
    let mut phase_last_gene = | gene: Option<Gene>| -> Result<(), Box<Error>> {
        if let Some(ref gene) = gene {
            if gene.biotype == "protein_coding" {
                phase_gene(
                    &gene, fasta_reader, &mut read_buffer,
                    &mut variant_buffer, fasta_writer,
                    window_len,
                    &mut refseq,
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

