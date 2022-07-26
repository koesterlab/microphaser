use std::cmp;
use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fs;
use std::io;
use std::rc::Rc;

use csv;
use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bcf};

use crate::common::{Interval, Variant};

pub fn bitvector_is_set(b: u64, k: usize) -> bool {
    (b & (1 << k)) != 0
}

pub fn switch_ascii_case(c: u8, r: u8) -> u8 {
    if r.is_ascii_uppercase() {
        c.to_ascii_lowercase()
    } else {
        c
    }
}

pub fn switch_ascii_case_vec(v: &[u8], r: u8) -> Vec<u8> {
    if r.is_ascii_uppercase() {
        v.to_ascii_lowercase()
    } else {
        v.to_ascii_uppercase()
    }
}

pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<dyn Error>> {
    match *variant {
        Variant::SNV { pos, alt, .. } => {
            //            debug!("Variant pos: {}", variant.pos());
            //            debug!("Read pos: {}", read.pos());
            //            debug!("Read end: {}", read.seq().len() + (read.pos() as usize));
            //            debug!("Read to check support: {}", String::from_utf8_lossy(read.qname()));
            let b = match read.cigar().read_pos(pos as u32, false, false) {
                Ok(None) => return Ok(false),
                Ok(Some(p)) => read.seq()[p as usize],
                _ => return Ok(false),
            };
            Ok(b == alt)
        }
        Variant::Insertion { .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                if let Cigar::Ins(_) = *c {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        Variant::Deletion { .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                if let Cigar::Del(_) = *c {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}

#[derive(Debug, Serialize)]
pub struct IDRecord {
    id: String,
    chrom: String,
    offset: u64,
    freq: f64,
    depth: u32,
    nvar: u32,
    nsomatic: u32,
    nvariant_sites: u32,
    nsomvariant_sites: u32,
    variant_sites: String,
    somatic_positions: String,
    somatic_aa_change: String,
    germline_positions: String,
    germline_aa_change: String,
    normal_sequence: String,
    mutant_sequence: String,
}

#[derive(Debug)]
pub struct Chrom {
    pub id: String,
    pub interval: Interval,
}

#[derive(Debug)]
pub struct HaplotypeSeq {
    sequence: Vec<u8>,
    record: IDRecord,
}

#[derive(Debug)]
pub struct Observation {
    read: Rc<bam::Record>,
    haplotype: u64,
}

impl Observation {
    pub fn update_haplotype(&mut self, i: usize, variant: &Variant) -> Result<(), Box<dyn Error>> {
        debug!(
            "Read pos {} ; variant pos {}",
            self.read.pos() as u32,
            variant.pos()
        );
        if (self.read.pos() as u64) > variant.pos() {
            panic!("bug: read starts right of variant");
        }
        if supports_variant(&self.read, &variant)? {
            debug!(
                "Read {} supports the variant at {}",
                String::from_utf8_lossy(self.read.qname()),
                variant.pos()
            );
            self.haplotype |= 1 << i;
            debug!("Haplotype: {}", self.haplotype)
        }

        Ok(())
    }
}

pub struct ObservationMatrix {
    observations: BTreeMap<u64, Vec<Observation>>,
    variants: VecDeque<Variant>,
}

impl Default for ObservationMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl ObservationMatrix {
    pub fn new() -> Self {
        ObservationMatrix {
            observations: BTreeMap::new(),
            variants: VecDeque::new(),
        }
    }

    /// Drain variants that are no longer in the window
    pub fn shrink_left(&mut self, k: usize) {
        debug!("self.variants length: {}", self.variants.len());
        debug!("range to drain: 0 - {}", k);
        self.variants.drain(..k);
        debug!("drained!");
        let mask = 2u64.pow(self.ncols()) - 1;
        for obs in Iterator::flatten(self.observations.values_mut()) {
            obs.haplotype &= mask;
        }
    }

    /// Add new variants
    pub fn extend_right(&mut self, new_variants: Vec<Variant>) -> Result<(), Box<dyn Error>> {
        let k = new_variants.len();
        debug!("Extend variants!");
        debug!("New variants {}", k);
        if k > 0 {
            for obs in Iterator::flatten(self.observations.values_mut()) {
                obs.haplotype <<= k;
            }
        }
        for obs in Iterator::flatten(self.observations.values_mut()) {
            for (i, variant) in new_variants.iter().rev().enumerate() {
                obs.update_haplotype(i, variant)?;
            }
        }
        self.variants.extend(new_variants.into_iter());

        Ok(())
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u64) {
        debug!(
            "Number of reads(before removal): {}",
            self.observations.len()
        );
        let observations = self.observations.split_off(&interval_end);
        self.observations = observations; //self.observations.split_off(&interval_end);
        debug!(
            "Number of reads(after removal): {}",
            self.observations.len()
        );
    }

    /// Check if read has already been added to observations - deprecated
    pub fn contains(&mut self, read: &bam::Record) -> Result<bool, Box<dyn Error>> {
        let end_pos = read.cigar().end_pos() as u64;
        let qname = read.qname();
        if self.observations.contains_key(&end_pos) {
            for obs in self.observations.get(&end_pos).unwrap() {
                if obs.read.qname() == qname {
                    return Ok(true);
                }
            }
            return Ok(false);
        }
        Ok(false)
    }

    /// Add read, while considering given interval end. TODO: Think about reads that are not overlapping variant
    pub fn push_read(
        &mut self,
        read: Rc<bam::Record>,
        interval_end: u64,
        interval_start: u64,
    ) -> Result<(), Box<dyn Error>> {
        let end_pos = read.cigar().end_pos() as u64;
        let start_pos = read.pos() as u64;
        debug!(
            "Read Start: {}, Read End: {} - Window Start: {}, Window End {}",
            start_pos, end_pos, interval_start, interval_end
        );
        if end_pos >= interval_end && start_pos <= interval_start {
            // only insert if end_pos is larger than the interval end
            let mut obs = Observation { read, haplotype: 0 };
            for (i, variant) in self.variants.iter().enumerate() {
                obs.update_haplotype(i, variant)?;
            }
            // debug!("Read: {}; haplotype: {}", String::from_utf8_lossy(obs.read.qname()), obs.haplotype);
            let pos = end_pos;
            self.observations
                .entry(pos)
                .or_insert_with(Vec::new)
                .push(obs);
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
        chrom: &str,
        offset: u64,
        window_len: u64,
        refseq: &[u8],
        fasta_writer: &mut fasta::Writer<O>,
        tsv_writer: &mut csv::Writer<fs::File>,
        normal_writer: &mut fasta::Writer<fs::File>,
        only_relevant: bool,
    ) -> Result<(), Box<dyn Error>> {
        let variants = self.variants.iter().collect_vec();
        // count haplotypes
        let mut haplotypes: VecMap<usize> = VecMap::new();
        for obs in Iterator::flatten(self.observations.values()) {
            debug!("obs {:?}", obs);
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let mut seq = Vec::with_capacity(window_len as usize);
        let mut germline_seq = Vec::with_capacity(window_len as usize);
        debug!("Printing at offset: {}", offset);
        debug!("refseq length {}", refseq.len());

        for (haplotype, count) in haplotypes.iter() {
            // VecMap forces usize as type for keys, but our haplotypes as u64
            let haplotype = haplotype as u64;
            let mut indel = false;
            debug!("Offset: {}", offset);
            debug!("Haplotype: {} ; count: {}", haplotype, count);
            debug!("Variants len: {}", variants.len());
            // build haplotype sequence
            seq.clear();
            germline_seq.clear();
            let mut n_somatic = 0;
            let mut n_variants = 0;
            let freq = *count as f64 / self.nrows() as f64;
            let depth = self.nrows() as u32;
            let mut i = offset;
            let mut j = 0;
            let mut window_end = offset + window_len;
            // Profile for all variants: 0 - reference, 1 - germline, 2 - somatic
            let mut variant_profile = Vec::new();
            //let mut somatic_profile = Vec::new();
            if variants.len() < 2 {
                germline_seq.extend(&refseq[offset as usize..(offset + window_len) as usize]);
                seq.extend(&refseq[offset as usize..(offset + window_len) as usize]);
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
                            if (j + 1) < variants.len() && i == variants[j + 1].pos() {
                                j += 1;
                            }
                            match *variants[j] {
                                // if SNV, we push the alternative base instead of the normal one, and change the case of the letter for visualisation
                                Variant::SNV { alt, .. } => {
                                    match variants[j].is_germline() {
                                        true => germline_seq
                                            .push(switch_ascii_case(alt, refseq[i as usize])),
                                        false => germline_seq.push(refseq[i as usize]),
                                    }
                                    seq.push(switch_ascii_case(alt, refseq[i as usize]));
                                    i += 1;
                                }
                                // if insertion, we insert the new bases (with changed case) and decrease the window-end, since we added bases and made the sequence longer
                                Variant::Insertion { seq: ref s, .. } => {
                                    match variants[j].is_germline() {
                                        true => germline_seq.extend(
                                            switch_ascii_case_vec(s, refseq[i as usize])
                                                .into_iter(),
                                        ),
                                        false => indel = true,
                                    }
                                    seq.extend(
                                        switch_ascii_case_vec(s, refseq[i as usize]).into_iter(),
                                    );
                                    i += 1;
                                    window_end -= (s.len() as u64) - 1;
                                }
                                // if deletion, we push the remaining base and increase the index to jump over the deleted bases. Then, we increase the window-end since we lost bases and need to fill up to 27.
                                Variant::Deletion { len, .. } => {
                                    match variants[j].is_germline() {
                                        true => germline_seq.push(refseq[i as usize]),
                                        false => indel = true,
                                    }
                                    seq.push(refseq[i as usize]);
                                    i += len + 1;
                                    window_end += len + 1;
                                }
                            }

                            // counting somatic variants
                            if !variants[j].is_germline() {
                                n_somatic += 1;
                                variant_profile.push(2);
                            } else {
                                variant_profile.push(1);
                            }
                            // counting total number of variants
                            n_variants += 1;
                        } else {
                            variant_profile.push(0);
                        }
                        j += 1;
                    }
                    // if no variant, just push the reference sequence
                    seq.push(refseq[i as usize]);
                    germline_seq.push(refseq[i as usize]);
                    i += 1
                }
            }
            // for indels, do not use the corresponding normal, but search for one with small hamming distance
            if indel {
                germline_seq.clear();
            }
            let mut shaid = sha1::Sha1::new();
            // generate unique haplotype ID containing position, transcript and sequence
            let id = format!("{:?}{}", &seq, offset);
            shaid.update(id.as_bytes());
            let fasta_id = (&shaid.digest().to_string()[..15]).to_string();
            // normal sequence of haplotype
            let normal_peptide = match germline_seq.len() {
                0 => String::from_utf8_lossy(&germline_seq),
                _ => String::from_utf8_lossy(&germline_seq[..window_len as usize]),
            };
            // neopeptide sequence
            let neopeptide = String::from_utf8_lossy(&seq[..window_len as usize]);
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
            let mut variantsites_pos_vec = Vec::new();
            // gather information iterating over the variants
            debug!("Variant profile len: {}", variant_profile.len());
            while c < variants.len() as u32 {
                if c < variant_profile.len() as u32 {
                    match variant_profile[c as usize] {
                        // somatic
                        2 => {
                            somatic_var_pos_vec.push(variants[c as usize].pos().to_string());
                            somatic_p_changes_vec.push(variants[c as usize].prot_change());
                        }
                        // germline
                        1 => {
                            germline_var_pos_vec.push(variants[c as usize].pos().to_string());
                            germline_p_changes_vec.push(variants[c as usize].prot_change());
                        }
                        // not present in this haplotype
                        _ => {}
                    }
                    // check if variant position is already in the variant_site list
                    if c == 0 || variants[c as usize].pos() != variants[(c - 1) as usize].pos() {
                        n_variantsites += 1;
                        variantsites_pos_vec.push(variants[c as usize].pos().to_string());
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
            let variantsites_pos = variantsites_pos_vec.join("|");
            // build the info record
            let record = IDRecord {
                id: fasta_id.to_owned(),
                chrom: chrom.to_owned(),
                offset,
                freq,
                depth,
                nvar: n_variants,
                nsomatic: n_somatic,
                nvariant_sites: n_variantsites as u32,
                nsomvariant_sites: n_som_variantsites as u32,
                variant_sites: variantsites_pos.to_owned(),
                somatic_positions: somatic_var_pos.to_owned(),
                somatic_aa_change: somatic_p_changes.to_owned(),
                germline_positions: germline_var_pos.to_owned(),
                germline_aa_change: germline_p_changes.to_owned(),
                normal_sequence: normal_peptide.to_owned().to_string(),
                mutant_sequence: neopeptide.to_owned().to_string(),
            };

            debug!(
                "relevant_check: {}, nvar: {}, freq: {} ",
                !(only_relevant),
                record.nvar > 0,
                record.freq < 1.00
            );
            debug!(
                "is_relevant: {}",
                !(only_relevant) || record.freq < 1.00 || record.nvar > 0
            );
            // write neopeptides, information and matching normal peptide to files
            if record.nvariant_sites > 1 {
                fasta_writer.write(&record.id.to_string(), None, &seq[..window_len as usize])?;
                if !germline_seq.is_empty() {
                    normal_writer.write(
                        &record.id.to_string(),
                        None,
                        &germline_seq[..window_len as usize],
                    )?;
                }
                tsv_writer.serialize(record)?;
            }
        }
        Ok(())
    }
}

pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    sequence: &fasta::Sequence,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::buffer::RecordBuffer,
    fasta_writer: &mut fasta::Writer<O>,
    tsv_writer: &mut csv::Writer<fs::File>,
    normal_writer: &mut fasta::Writer<fs::File>,
    window_len: u64,
    refseq: &mut Vec<u8>,
    only_relevant: bool,
    unsupported_alleles_warning_only: bool,
) -> Result<(), Box<dyn Error>> {
    let mut chunk = 0;
    while chunk < sequence.len - 1000000 {
        fasta_reader.fetch(
            &sequence.name,
            chunk,
            cmp::min(chunk + 1000000, sequence.len - 1),
        )?;
        fasta_reader.read(refseq)?;
        let mut variant_tree = BTreeMap::new();
        let mut read_tree = BTreeMap::new();
        debug!("Start Phasing");
        read_buffer.fetch(&sequence.name.as_bytes(), chunk, chunk + 1000000)?;
        let chrom = &sequence.name;
        let mut max_read_len = 50;
        // load read buffer into BTree
        for rec in read_buffer.iter() {
            if rec.seq().len() as u64 > max_read_len {
                max_read_len = rec.seq().len() as u64;
            }
            read_tree
                .entry(rec.pos() as u64)
                .or_insert_with(Vec::new)
                .push(rec.clone())
        }
        debug!("Reads Tree length: {}", read_tree.len());

        // load variant buffer into BTree
        let (_addedvars, _deletedvars) =
            variant_buffer.fetch(&sequence.name.as_bytes(), chunk, chunk + 1000000)?;
        let _vars = variant_buffer
            .iter_mut()
            .map(|rec| variant_tree.insert(rec.pos() as u64, Variant::new(rec, unsupported_alleles_warning_only).unwrap()))
            .collect_vec();

        let mut observations = ObservationMatrix::new();
        let mut frameshifts = BTreeMap::new();
        frameshifts.insert(0, 0);
        let mut offset = chunk;
        let mut old_offset = offset;
        loop {
            let valid = offset + window_len <= chunk + 1000000;
            if !valid {
                break;
            }
            debug!("Offset {}, old offset {}", offset, old_offset);
            // advance window to next position
            let nvars = Iterator::flatten(
                variant_tree
                    .range(offset..(offset + window_len))
                    .map(|var| var.1),
            )
            .count();
            debug!("Variants in window: {}", nvars);
            // first window in the exon, all variants found are newly added
            let added_vars = if offset == old_offset {
                nvars
            // if we advance the window (forward or reverse), just the newly added variants are counted
            } else {
                Iterator::flatten(
                    variant_tree
                        .range((old_offset + window_len)..(offset + window_len))
                        .map(|var| var.1),
                )
                .count()
            };

            // first window, no variants are deleted
            let deleted_vars = if offset == old_offset {
                0
            } else {
                Iterator::flatten(variant_tree.range(old_offset..offset).map(|var| var.1)).count()
            };

            debug!(
                "Offset: {} - max_read_len - window_len {}",
                offset,
                (max_read_len - window_len)
            );
            let reads = {
                // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                if offset == 0 {
                    Iterator::flatten(read_tree.range(offset..(offset + 1)).map(|rec| rec.1))
                        .collect_vec()
                }
                // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
                else {
                    Iterator::flatten(read_tree.range(offset..(offset + 1)).map(|rec| rec.1))
                        .collect_vec()
                }
            };

            debug!("Variants added: {}", added_vars);
            debug!("Variants deleted: {}", deleted_vars);
            {
                // delete rows
                observations.cleanup_reads(offset + window_len);
                // delete columns
                observations.shrink_left(deleted_vars);

                // add new reads
                debug!("Reads: {}", reads.len());
                for read in reads {
                    observations.push_read(read.clone(), offset + window_len, offset)?;
                }

                // collect variants
                let variants = Iterator::flatten(
                    variant_tree
                        .range_mut(offset..(offset + window_len))
                        .map(|var| var.1.clone()),
                )
                .skip(nvars - added_vars)
                .collect_vec();
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
                    let coding_shift = offset;
                    debug!("Offset - Exonstart % 3: {}", coding_shift % 3);
                    debug!("Shift: {}", frameshift);
                    if coding_shift % 3 == frameshift {
                        // print haplotypes
                        debug!("Should print haplotypes");
                        observations
                            .print_haplotypes(
                                chrom,
                                offset,
                                window_len,
                                refseq,
                                fasta_writer,
                                tsv_writer,
                                normal_writer,
                                only_relevant,
                            )
                            .unwrap();
                    }
                }

                old_offset = offset;
                offset += 1;
            }
        }
        chunk += 1000000;
    }
    // fasta_reader.fetch(&sequence.name, 0 as u64, sequence.len - 1)?;
    // fasta_reader.read(refseq)?;
    // let mut variant_tree = BTreeMap::new();
    // let mut read_tree = BTreeMap::new();
    // debug!("Start Phasing");
    // read_buffer.fetch(&sequence.name.as_bytes(), 0, sequence.len as u32- 1)?;
    // let chrom = &sequence.name;
    // let mut max_read_len = 50 as u32;
    // // load read buffer into BTree
    // for rec in read_buffer.iter() {
    //     if rec.seq().len() as u32 > max_read_len {
    //     max_read_len = rec.seq().len() as u32;}
    //     read_tree.entry(rec.pos() as u32).or_insert_with(|| Vec::new()).push(rec.clone())
    // }
    // debug!("Reads Tree length: {}", read_tree.len());
    //
    // // load variant buffer into BTree
    // let (_addedvars, _deletedvars) = variant_buffer.fetch(&sequence.name.as_bytes(), 0, sequence.len as u32 - 1)?;
    // let _vars = variant_buffer.iter_mut().map(|rec| variant_tree.insert(rec.pos(), Variant::new(rec).unwrap())).collect_vec();
    //
    //
    //
    // let mut observations = ObservationMatrix::new();
    // let mut frameshifts = BTreeMap::new();
    // frameshifts.insert(0, 0);
    // let mut offset = 0;
    // let mut old_offset = offset;
    // loop {
    //     let valid = offset + window_len <= sequence.len  as u32 - 1;
    //     if !valid {
    //         break;
    //     }
    //     debug!("Offset {}, old offset {}", offset, old_offset);
    //     // advance window to next position
    //     let nvars = Iterator::flatten(variant_tree.range(offset..(offset + window_len)).map(|var| var.1)).count();
    //     debug!("Variants in window: {}",nvars);
    //     // first window in the exon, all variants found are newly added
    //     let added_vars = if offset == old_offset {
    //         nvars
    //     // if we advance the window (forward or reverse), just the newly added variants are counted
    //     } else {
    //         Iterator::flatten(variant_tree.range((old_offset + window_len)..(offset + window_len)).map(|var| var.1)).count()
    //     };
    //
    //     // first window, no variants are deleted
    //     let deleted_vars = if offset == old_offset {
    //         0
    //     } else {
    //         Iterator::flatten(variant_tree.range(old_offset..offset).map(|var| var.1)).count()
    //     };
    //
    //     debug!("Offset: {} - max_read_len - window_len {}", offset, (max_read_len - window_len));
    //     let reads = {
    //         // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
    //         if offset == 0 {
    //             Iterator::flatten(read_tree.range(offset..(offset+1)).map(|rec| rec.1)).collect_vec()
    //         }
    //         // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
    //         else {
    //             Iterator::flatten(read_tree.range((offset-1)..(offset+1)).map(|rec| rec.1)).collect_vec()
    //         }
    //     };
    //
    //
    //
    //     debug!("Variants added: {}", added_vars);
    //     debug!("Variants deleted: {}", deleted_vars);
    //     {
    //         // delete rows
    //         observations.cleanup_reads(offset + window_len);
    //         // delete columns
    //         observations.shrink_left(deleted_vars);
    //
    //         // add new reads
    //         debug!("Reads: {}", reads.len());
    //         for read in reads {
    //             observations.push_read(read.clone(), offset + window_len, offset)?;
    //         }
    //
    //
    //         // collect variants
    //         let variants = Iterator::flatten(variant_tree.range_mut(offset..(offset + window_len)
    //             ).map(|var| var.1.clone())).skip(nvars - added_vars).collect_vec();
    //         debug!("Variants(after deleting and adding): {}", variants.len());
    //         // determine frameshifts
    //         for variant in &variants {
    //             debug!("Variants!");
    //             debug!("Variant Pos: {}", variant.pos());
    //             let s = variant.frameshift();
    //             if s > 0 {
    //                 let previous = frameshifts.values().map(|prev| prev + s).collect_vec();
    //                 for s_ in previous {
    //                     frameshifts.insert(variant.end_pos(), s + s_);
    //                 }
    //             }
    //         }
    //
    //         // add columns
    //         observations.extend_right(variants)?;
    //
    //         for (_, &frameshift) in frameshifts.range(..offset) {
    //             // possible shift if exon starts with the rest of a split codon (splicing)
    //             let coding_shift = offset;
    //             debug!("Offset - Exonstart % 3: {}", coding_shift % 3);
    //             debug!("Shift: {}", frameshift);
    //             if coding_shift % 3 == frameshift {
    //                 // print haplotypes
    //                 debug!("Should print haplotypes");
    //                 observations.print_haplotypes(
    //                     chrom, offset, window_len, refseq, fasta_writer, tsv_writer, normal_writer, only_relevant
    //                 ).unwrap();
    //             }
    //         }
    //
    //         old_offset = offset;
    //         offset += 1;
    //     }
    // }
    Ok(())
}

pub fn phase<F: io::Read + io::Seek, O: io::Write>(
    fasta_reader: &mut fasta::IndexedReader<F>,
    bcf_reader: bcf::Reader,
    bam_reader: bam::IndexedReader,
    fasta_writer: &mut fasta::Writer<O>,
    tsv_writer: &mut csv::Writer<fs::File>,
    normal_writer: &mut fasta::Writer<fs::File>,
    window_len: u64,
    only_relevant: bool,
    unsupported_alleles_warning_only: bool,
) -> Result<(), Box<dyn Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader, false);
    let mut variant_buffer = bcf::buffer::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence
    debug!("refseq length {}", refseq.len());
    let seqs = fasta_reader.index.sequences();
    for s in seqs {
        phase_gene(
            &s,
            fasta_reader,
            &mut read_buffer,
            &mut variant_buffer,
            fasta_writer,
            tsv_writer,
            normal_writer,
            window_len,
            &mut refseq,
            only_relevant,
            unsupported_alleles_warning_only,
        )?;
    }
    Ok(())
}
