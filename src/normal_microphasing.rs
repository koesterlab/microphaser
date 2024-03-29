use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fs;
use std::io;
use std::rc::Rc;

use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;

use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bcf};

use bio_types::strand::Strand;

use crate::common::{Gene, Interval, PhasingStrand, Transcript, Variant};

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
            let b = match read.cigar().read_pos(pos as u32, false, false) {
                Ok(None) => return Ok(false),
                Ok(Some(p)) => read.seq()[p as usize],
                _ => return Ok(false),
            };
            Ok(b == alt)
        }
        Variant::Insertion { len, .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                if let Cigar::Ins(_) = *c {
                    match c.len() == len as u32 {
                        true => return Ok(true),
                        false => (),
                    }
                }
            }
            Ok(false)
        }
        Variant::Deletion { len, .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                if let Cigar::Del(_) = *c {
                    match c.len() == len as u32 {
                        true => return Ok(true),
                        false => (),
                    }
                }
            }
            Ok(false)
        }
    }
}

#[derive(Debug, Serialize)]
pub struct IDRecord {
    id: String,
    transcript: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    offset: u64,
    frame: u64,
    freq: f64,
    depth: u32,
    nvar: u32,
    nsomatic: u32,
    nvariant_sites: u32,
    nsomvariant_sites: u32,
    strand: String,
    variant_sites: String,
    somatic_positions: String,
    somatic_aa_change: String,
    germline_positions: String,
    germline_aa_change: String,
    peptide_sequence: String,
}

impl IDRecord {
    pub fn update(&self, rec: &IDRecord, offset: u64, seq: Vec<u8>) -> Self {
        let mut shaid = sha1::Sha1::new();
        // generate unique haplotype ID containing position, transcript and sequence
        let id = format!("{:?}{}{}", &seq, &self.transcript, offset);
        shaid.update(id.as_bytes());
        let fasta_id = format!(
            "{}{}",
            &shaid.digest().to_string()[..15],
            self.strand.chars().next().unwrap()
        );
        let mut somatic_positions = self.somatic_positions.clone();
        somatic_positions.push_str(&rec.somatic_positions);
        let mut somatic_aa_change = self.somatic_aa_change.clone();
        somatic_aa_change.push_str(&rec.somatic_aa_change);
        let mut germline_positions = self.germline_positions.clone();
        germline_positions.push_str(&rec.germline_positions);
        let mut germline_aa_change = self.germline_aa_change.clone();
        germline_aa_change.push_str(&rec.germline_aa_change);
        debug!("nvars {} {}", self.nvar, rec.nvar);
        IDRecord {
            id: fasta_id,
            transcript: self.transcript.to_owned(),
            gene_id: self.gene_id.to_owned(),
            gene_name: self.gene_name.to_owned(),
            chrom: self.chrom.to_owned(),
            offset: offset + self.offset,
            frame: self.frame,
            freq: self.freq * rec.freq,
            depth: self.depth,
            nvar: self.nvar + rec.nvar,
            nsomatic: self.nsomatic + rec.nsomatic,
            nvariant_sites: self.nvariant_sites + rec.nvariant_sites,
            nsomvariant_sites: self.nsomvariant_sites + rec.nsomvariant_sites,
            strand: self.strand.to_owned(),
            variant_sites: self.variant_sites.to_owned() + &rec.variant_sites,
            somatic_positions,
            somatic_aa_change,
            germline_positions,
            germline_aa_change,
            peptide_sequence: String::from_utf8(seq).unwrap(),
        }
    }

    pub fn add_freq(&self, freq: f64) -> Self {
        let new_nvar = match freq {
            f if f > 0.0 => self.nvar - 1,
            _ => self.nvar,
        };
        let new_somatic = match new_nvar < self.nsomatic {
            true => self.nsomatic - 1,
            false => self.nsomatic,
        };
        IDRecord {
            id: self.id.to_owned(),
            transcript: self.transcript.to_owned(),
            gene_id: self.gene_id.to_owned(),
            gene_name: self.gene_name.to_owned(),
            chrom: self.chrom.to_owned(),
            offset: self.offset,
            frame: self.frame,
            freq: self.freq + freq,
            depth: self.depth,
            nvar: new_nvar,
            nsomatic: new_somatic,
            nvariant_sites: self.nvariant_sites,
            nsomvariant_sites: self.nsomvariant_sites,
            strand: self.strand.to_owned(),
            variant_sites: self.variant_sites.to_owned(),
            somatic_positions: self.somatic_positions.to_owned(),
            somatic_aa_change: self.somatic_aa_change.to_owned(),
            germline_positions: self.germline_positions.to_owned(),
            germline_aa_change: self.germline_aa_change.to_owned(),
            peptide_sequence: self.peptide_sequence.to_owned(),
        }
    }
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
            self.read.pos() as u64,
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
    pub fn cleanup_reads(&mut self, interval_end: u64, reverse: bool) {
        debug!(
            "Number of reads(before removal): {}",
            self.observations.len()
        );
        let observations = self.observations.split_off(&interval_end);
        if !reverse {
            self.observations = observations; //self.observations.split_off(&interval_end);
        }
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
        reverse: bool,
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
            let pos = match reverse {
                true => start_pos,
                false => end_pos,
            };
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
        gene: &Gene,
        transcript: &Transcript,
        offset: u64,
        splice_end: u64,
        splice_pos: u64,
        splice_gap: u64,
        exon_end: u64,
        exon_start: u64,
        window_len: u64,
        refseq: &[u8],
        tsv_writer: &mut csv::Writer<fs::File>,
        fasta_writer: &mut fasta::Writer<O>,
        is_short_exon: bool,
        frame: u64,
    ) -> Result<Vec<HaplotypeSeq>, Box<dyn Error>> {
        let variants_forward = self.variants.iter().collect_vec();
        let mut variants_reverse = variants_forward.clone();
        variants_reverse.reverse();
        let variants = match transcript.strand {
            PhasingStrand::Reverse => variants_reverse,
            PhasingStrand::Forward => variants_forward,
        };
        // count haplotypes
        let mut haplotypes: VecMap<usize> = VecMap::new();
        for obs in Iterator::flatten(self.observations.values()) {
            debug!("obs {:?}", obs);
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let splice_window_len = splice_end - offset;
        let mut seq = Vec::with_capacity(splice_window_len as usize);
        debug!("Gene Start {}, Gene End {}", gene.start(), gene.end());
        debug!("Printing at offset: {}", offset);
        debug!("refseq length {}", refseq.len());

        // Strand orientation
        let strand = match transcript.strand {
            PhasingStrand::Reverse => "Reverse",
            PhasingStrand::Forward => "Forward",
        };

        //frameshift
        debug!("Frame {}", frame);
        let mut haplotypes_vec = Vec::new();
        // If there are no reads covering the window, fill with normal sequence and mark - important for gaps in coverage (see frameshifts)
        if haplotypes.is_empty() {
            haplotypes.insert(0, 0);
        }
        for (haplotype, count) in haplotypes.iter() {
            // VecMap forces usize as type for keys, but our haplotypes as u64
            let haplotype = haplotype as u64;
            debug!("Offset: {}", offset);
            debug!("Haplotype: {} ; count: {}", haplotype, count);
            debug!("Variants len: {}", variants.len());
            // build haplotype sequence
            seq.clear();
            let mut insertion = false;
            let mut n_somatic = 0;
            let mut n_variants = 0;
            let freq = *count as f64 / self.nrows() as f64;
            let depth = self.nrows() as u32;
            let mut i = offset;
            let mut j = 0;
            let mut window_end = splice_end;
            debug!("window end: {}", window_end);
            // Profile for all variants: 0 - reference, 1 - germline, 2 - somatic
            let mut variant_profile = Vec::new();
            if variants.is_empty() {
                seq.extend(
                    &refseq[(offset - gene.start()) as usize..(window_end - gene.start()) as usize],
                );
            } else {
                while i < window_end {
                    debug!("window_end: {}", window_end);
                    debug!("i: {}", i);
                    debug!("j: {}", j);
                    // TODO what happens if a deletion starts upstream of window and overlaps it

                    while j < variants.len() && i == variants[j].pos() {
                        debug!("j: {}, variantpos: {}", j, variants[j].pos());
                        if (freq - 1.0).abs() < f64::EPSILON && !(variants[j].is_germline()) {
                            j += 1;
                            variant_profile.push(0);
                            continue;
                        }
                        if bitvector_is_set(haplotype, j) {
                            debug!("Haplotype: {} ; j: {}", haplotype, j);
                            if (j + 1) < variants.len() && i == variants[j + 1].pos() {
                                j += 1;
                            }
                            match *variants[j] {
                                // if SNV, we push the alternative base instead of the normal one, and change the case of the letter for visualisation
                                Variant::SNV { alt, .. } => {
                                    seq.push(switch_ascii_case(
                                        alt,
                                        refseq[(i - gene.start()) as usize],
                                    ));
                                    i += 1;
                                }
                                // if insertion, we insert the new bases (with changed case) and decrease the window-end, since we added bases and made the sequence longer
                                Variant::Insertion { seq: ref s, .. } => {
                                    seq.extend(
                                        switch_ascii_case_vec(
                                            s,
                                            refseq[(i - gene.start()) as usize],
                                        )
                                        .into_iter(),
                                    );
                                    insertion = true;
                                    i += 1;
                                }
                                // if deletion, we push the remaining base and increase the index to jump over the deleted bases. Then, we increase the window-end since we lost bases and need to fill up to 27.
                                Variant::Deletion { len, .. } => {
                                    seq.push(refseq[(i - gene.start()) as usize]);
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
                    seq.push(refseq[(i - gene.start()) as usize]);
                    i += 1
                }
            }
            let this_window_len = match seq.len() < window_len as usize {
                true => seq.len() as u64,
                false => window_len,
            };
            debug!("this window len {}", this_window_len);
            let peptide = match splice_pos {
                1 => String::from_utf8_lossy(&seq[splice_gap as usize..]),
                0 => match insertion {
                    true => String::from_utf8_lossy(&seq),
                    false => String::from_utf8_lossy(&seq[..this_window_len as usize]),
                },
                _ => String::from_utf8_lossy(&seq),
            };
            let stop_gain = match transcript.strand {
                PhasingStrand::Forward => {
                    peptide.starts_with("TGA")
                        || peptide.starts_with("TAG")
                        || peptide.starts_with("TAA")
                }
                PhasingStrand::Reverse => {
                    peptide.ends_with("TCA") || peptide.ends_with("CTA") || peptide.ends_with("TTA")
                }
            };
            if stop_gain && splice_pos != 2 {
                // if the peptide is not in the correct reading frame because of leftover bases, we do not care about the stop codon since it is not in the ORF
                debug!("Peptide with STOP codon: {}", peptide);
                continue;
            }

            let mut shaid = sha1::Sha1::new();
            // generate unique haplotype ID containing position, transcript and sequence
            let id = format!("{:?}{}{}", &seq, transcript.id, offset);
            shaid.update(id.as_bytes());
            let fasta_id = format!(
                "{}{}",
                &shaid.digest().to_string()[..15],
                strand.chars().next().unwrap()
            );
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
                transcript: transcript.id.to_owned(),
                gene_id: gene.id.to_owned(),
                gene_name: gene.name.to_owned(),
                chrom: gene.chrom.to_owned(),
                offset,
                frame,
                freq,
                depth,
                nvar: n_variants,
                nsomatic: n_somatic,
                nvariant_sites: n_variantsites as u32,
                nsomvariant_sites: n_som_variantsites as u32,
                strand: strand.to_string(),
                variant_sites: variantsites_pos.to_owned(),
                somatic_positions: somatic_var_pos.to_owned(),
                somatic_aa_change: somatic_p_changes.to_owned(),
                germline_positions: germline_var_pos.to_owned(),
                germline_aa_change: germline_p_changes.to_owned(),
                peptide_sequence: peptide.to_owned().to_string(),
            };

            // make haplotype record to carry over to next exon
            debug!("exon end: {}, window end {}", exon_end, window_end);
            let rest = match window_end > exon_end {
                true => 0,
                false => exon_end - window_end,
            };
            debug!("Rest: {}", rest);
            let start = offset - exon_start;
            debug!("Start: {}", start);

            let newseq = &seq;

            let hap_seq = HaplotypeSeq {
                sequence: newseq.to_vec(),
                record: IDRecord {
                    id: fasta_id.to_owned(),
                    transcript: transcript.id.to_owned(),
                    gene_id: gene.id.to_owned(),
                    gene_name: gene.name.to_owned(),
                    chrom: gene.chrom.to_owned(),
                    offset,
                    frame,
                    freq,
                    depth,
                    nvar: n_variants,
                    nsomatic: n_somatic,
                    nvariant_sites: n_variantsites as u32,
                    nsomvariant_sites: n_som_variantsites as u32,
                    strand: strand.to_string(),
                    variant_sites: variantsites_pos.to_owned(),
                    somatic_positions: somatic_var_pos.to_owned(),
                    somatic_aa_change: somatic_p_changes.to_owned(),
                    germline_positions: germline_var_pos.to_owned(),
                    germline_aa_change: germline_p_changes.to_owned(),
                    peptide_sequence: String::from_utf8_lossy(&newseq).to_string(),
                },
            };

            haplotypes_vec.push(hap_seq);
            // print haplotypes
            if !(is_short_exon) {
                match splice_pos {
                    1 => fasta_writer.write(
                        &record.id.to_string(),
                        None,
                        &seq[splice_gap as usize..],
                    )?,
                    0 => fasta_writer.write(
                        &record.id.to_string(),
                        None,
                        &seq[..window_len as usize],
                    )?,
                    _ => {}
                };
                tsv_writer.serialize(record)?;
            }
        }
        Ok(haplotypes_vec)
    }
}

pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    gene: &Gene,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::buffer::RecordBuffer,
    tsv_writer: &mut csv::Writer<fs::File>,
    fasta_writer: &mut fasta::Writer<O>,
    window_len: u64,
    refseq: &mut Vec<u8>,
    unsupported_allele_warning_only: bool,
) -> Result<(), Box<dyn Error>> {
    // if an exon is near to the gene end, a deletion could cause refseq to overflow, so we increase the length of refseq
    let end_overflow = 100;
    fasta_reader.fetch(
        &gene.chrom,
        gene.start() as u64,
        (gene.end() + end_overflow) as u64,
    )?;
    fasta_reader.read(refseq)?;
    let mut variant_tree = BTreeMap::new();
    let mut read_tree = BTreeMap::new();
    debug!("Start Phasing");
    read_buffer.fetch(&gene.chrom.as_bytes(), gene.start(), gene.end())?;

    let mut max_read_len = 0;
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
        variant_buffer.fetch(&gene.chrom.as_bytes(), gene.start(), gene.end())?;
    let _vars = variant_buffer
        .iter_mut()
        .map(|rec| {
            variant_tree.insert(
                rec.pos() as u64,
                Variant::new(rec, unsupported_allele_warning_only).unwrap(),
            )
        })
        .collect_vec();

    for transcript in &gene.transcripts {
        debug!("is coding: {}", transcript.is_coding());
        if !(transcript.is_coding()) {
            continue;
        }
        let exon_number = transcript.exons.len();
        debug!("Number of Exons in Transcript: {}", exon_number);
        debug!("Transcript strand orientation: {:?}", transcript.strand);
        let mut observations = ObservationMatrix::new();
        let mut frameshifts = BTreeMap::new();
        match transcript.strand {
            PhasingStrand::Forward => frameshifts.insert(0, 0),
            PhasingStrand::Reverse => frameshifts.insert(gene.end(), 0),
        };
        // Possible rest of an exon that does not form a complete codon yet
        let mut exon_rest = 0;

        // Haplotypes from the end of the last exon
        let mut prev_hap_vec: Vec<HaplotypeSeq> = Vec::new();
        let mut hap_vec: Vec<HaplotypeSeq> = Vec::new();
        // Possible variants that are still in the BTree after leaving the last exon (we do not drain variants if we leave an exon)
        let mut last_window_vars = 0;
        for (exon_count, exon) in transcript.exons.iter().enumerate() {
            // if all frames are closed, finish the transcript
            if frameshifts.is_empty() {
                break;
            }
            debug!("Exon Start: {}", exon.start);
            debug!("Exon End: {}", exon.end);

            if exon.start > exon.end {
                continue;
            }
            let is_last_exon = exon_count as usize == exon_number - 1;
            let is_first_exon = exon_count == 0;
            debug!("Exon Length: {}", exon.end - exon.start);
            let exon_len = exon.end - exon.start;
            debug!("Exon Rest: {}", exon_rest);
            // Possible offset at the exon start, first nucleotides could be part of a codon started in the previous exon
            let current_exon_offset = match exon_rest {
                0 => 0,
                _ => 3 - exon_rest,
            };
            debug!("Exon Offset: {}", current_exon_offset);
            let is_short_exon = match exon_len < 3 {
                true => true,
                false => {
                    window_len >= exon_len - current_exon_offset - (3 - current_exon_offset) % 3
                }
            };
            debug!("Short Exon: {}", is_short_exon);
            // if the exon is shorter than the window, we need to fix the window len for this exon
            let mut exon_window_len = match is_short_exon {
                false => window_len,
                true => (exon_len - current_exon_offset) - ((exon_len - current_exon_offset) % 3),
            };
            if exon_window_len == 0 {
                exon_window_len = exon_len
            }
            // let mut exon_window_len = match is_short_exon {
            //     false => window_len,
            //     true => (exon_len - current_exon_offset) - ((exon_len - current_exon_offset) % 3),
            // };
            // if exon_window_len == 0 {
            //     exon_window_len = exon_len
            // }
            exon_rest = 0;
            let mut offset = if transcript.strand == PhasingStrand::Reverse {
                exon.end - exon_window_len - current_exon_offset
            } else {
                exon.start + current_exon_offset
            };
            let mut reached_end = false;
            debug!("Exon window lenght: {}", exon_window_len);
            debug!("Starting Offset of the Exon: {}", offset);
            let mut old_offset = offset;
            let mut old_end = old_offset + exon_window_len;
            debug!("Variants left from previous Exon: {}", last_window_vars);
            observations.shrink_left(last_window_vars);
            last_window_vars = 0;
            let mut is_first_exon_window = true;
            loop {
                // if all frames are closed, finish the transcript
                if frameshifts.is_empty() {
                    break;
                }
                let valid = match transcript.strand {
                    PhasingStrand::Reverse => offset >= exon.start,
                    PhasingStrand::Forward => offset + exon_window_len <= exon.end,
                };
                if !valid {
                    break;
                }
                if max_read_len < exon_window_len {
                    break;
                }

                let rest = match transcript.strand {
                    PhasingStrand::Forward => exon.end - (offset + exon_window_len),
                    PhasingStrand::Reverse => offset - exon.start,
                };

                let is_last_exon_window = rest < 3;
                if is_last_exon_window {
                    debug!("Last exon window");
                }
                let (splice_side_offset, splice_end, splice_gap, splice_pos) =
                    match transcript.strand {
                        PhasingStrand::Forward => {
                            if is_short_exon {
                                (
                                    offset - current_exon_offset,
                                    offset + exon_window_len + rest,
                                    current_exon_offset + rest,
                                    2,
                                )
                            } else if is_first_exon_window {
                                if is_last_exon_window {
                                    (
                                        offset - current_exon_offset,
                                        offset + exon_window_len + rest,
                                        current_exon_offset + rest,
                                        2,
                                    )
                                } else {
                                    (
                                        offset - current_exon_offset,
                                        offset + exon_window_len,
                                        current_exon_offset,
                                        1,
                                    )
                                }
                            } else if is_last_exon_window {
                                (offset, offset + exon_window_len + rest, rest, 0)
                            } else {
                                (offset, offset + exon_window_len, 0, 0)
                            }
                        }
                        PhasingStrand::Reverse => {
                            if is_short_exon {
                                (
                                    offset - rest,
                                    offset + exon_window_len + current_exon_offset,
                                    current_exon_offset + rest,
                                    2,
                                )
                            } else if is_first_exon_window {
                                (
                                    offset,
                                    offset + exon_window_len + current_exon_offset,
                                    current_exon_offset,
                                    0,
                                )
                            } else if is_last_exon_window {
                                (offset - rest, offset + exon_window_len, rest, 1)
                            } else {
                                (offset, offset + exon_window_len, 0, 0)
                            }
                        }
                    };
                debug!(
                    "Splice_offset {}, Offset {}, old offset {}",
                    splice_side_offset, offset, old_offset
                );
                debug!("WinStart {}, WinEnd {}", splice_side_offset, splice_end);

                // advance window to next position
                let nvars = Iterator::flatten(
                    variant_tree
                        .range(splice_side_offset..splice_end)
                        .map(|var| var.1),
                )
                .count();
                // store number of variants in window in case it is the last window for this exon
                last_window_vars = nvars;
                debug!("Variants in window: {}", nvars);
                // first window in the exon, all variants found are newly added
                let added_vars = if is_first_exon_window {
                    nvars
                } else if is_short_exon {
                    0
                // has already added last variants in a previous iterations
                } else if reached_end {
                    debug!("End of Exon reached");
                    0
                // if we advance the window (forward or reverse), just the newly added variants are counted
                // forward orientation
                } else if splice_side_offset > old_offset {
                    debug!(
                        "variant range: {:?}",
                        variant_tree.range(old_end..splice_end)
                    );
                    Iterator::flatten(variant_tree.range((old_end)..(splice_end)).map(|var| var.1))
                        .count()
                // reverse orientation
                } else {
                    Iterator::flatten(
                        variant_tree
                            .range(splice_side_offset..old_offset)
                            .map(|var| var.1),
                    )
                    .count()
                };

                // first window in the exon, no variants are deleted
                let deleted_vars = if offset == old_offset || is_short_exon {
                    0
                // if we advance the window (forward or reverse), we will delete all variants that drop out of the window bounds
                // forward orientation
                } else if splice_side_offset > old_offset {
                    Iterator::flatten(
                        variant_tree
                            .range(old_offset..splice_side_offset)
                            .map(|var| var.1),
                    )
                    .count()
                // reverse orientation
                } else {
                    Iterator::flatten(
                        variant_tree
                            .range((splice_end)..old_end) //(old_offset + exon_window_len))
                            .map(|var| var.1),
                    )
                    .count()
                };

                // Final variants have been added in this iteration
                if is_last_exon_window {
                    reached_end = true;
                }

                debug!("Offset + wlen: {}", offset + exon_window_len);
                debug!("Old offset + wlen: {}", old_offset + exon_window_len);
                debug!(
                    "Offset: {} - max_read_len {} - window_len {}",
                    offset, max_read_len, exon_window_len
                );
                debug!(
                    "Offset: {} ; (Offset - (max_read_len - window_len)) = {}",
                    splice_side_offset,
                    splice_side_offset - (max_read_len - exon_window_len)
                );
                let reads = if transcript.strand == PhasingStrand::Reverse {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.end - exon_window_len - current_exon_offset {
                        debug!("First exon window");
                        Iterator::flatten(
                            read_tree
                                .range(
                                    (splice_side_offset - (max_read_len - exon_window_len))
                                        ..(splice_side_offset + 1),
                                )
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                    // while advancing the window (reverse orientation), we only add reads that end in the range between old and new window end, so we don't count any read twice
                    else {
                        Iterator::flatten(
                            read_tree
                                .range(
                                    (splice_side_offset - (max_read_len - exon_window_len))
                                        ..(splice_side_offset + 1),
                                )
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                } else {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.start + current_exon_offset {
                        Iterator::flatten(
                            read_tree
                                .range(
                                    (splice_side_offset - (max_read_len - exon_window_len))
                                        ..(splice_side_offset + 1),
                                )
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                    // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
                    else {
                        Iterator::flatten(
                            read_tree
                                .range((splice_side_offset)..(splice_side_offset + 1))
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                };

                debug!("Variants added: {}", added_vars);
                debug!("Variants deleted: {}", deleted_vars);
                {
                    // delete rows
                    let reverse = match transcript.strand {
                        PhasingStrand::Reverse => true,
                        PhasingStrand::Forward => false,
                    };
                    if reverse {
                        observations.cleanup_reads(splice_side_offset, reverse);
                    } else {
                        observations.cleanup_reads(splice_end, reverse);
                    }
                    // delete columns
                    observations.shrink_left(deleted_vars);

                    // add new reads
                    debug!("Reads: {}", reads.len());
                    for read in reads {
                        observations.push_read(
                            read.clone(),
                            splice_end,
                            splice_side_offset,
                            reverse,
                        )?;
                    }

                    // collect variants
                    let variants = match transcript.strand {
                        PhasingStrand::Reverse => Iterator::flatten(
                            variant_tree
                                .range_mut(splice_side_offset..splice_end)
                                .rev()
                                .map(|var| var.1.clone()),
                        )
                        .skip(nvars - added_vars)
                        .collect_vec(),
                        PhasingStrand::Forward => Iterator::flatten(
                            variant_tree
                                .range_mut(splice_side_offset..splice_end)
                                .map(|var| var.1.clone()),
                        )
                        .skip(nvars - added_vars)
                        .collect_vec(),
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
                                frameshifts.insert(variant.end_pos(), s_);
                            }
                        }
                    }

                    // add columns
                    observations.extend_right(variants)?;

                    let mut stopped_frameshift = 3;
                    let mut active_frameshifts = match transcript.strand {
                        PhasingStrand::Forward => frameshifts.range(..offset),
                        PhasingStrand::Reverse => frameshifts.range(offset + exon_window_len..),
                    };
                    debug!("Active frameshifts: {:?}", active_frameshifts);
                    let mut frameshift_count = 0;
                    let mut main_orf = false;
                    for (&key, &frameshift) in &mut active_frameshifts {
                        // possible shift if exon starts with the rest of a split codon (splicing)
                        debug!("Frameshift: {}", frameshift);
                        if frameshift == 0 {
                            main_orf = true;
                        }
                        frameshift_count += 1;
                        let coding_shift = match transcript.strand {
                            PhasingStrand::Forward => offset - exon.start,
                            PhasingStrand::Reverse => exon.end - offset,
                        };
                        let has_frameshift = frameshift > 0;
                        debug!("Coding Shift: {}", coding_shift);
                        debug!("Current Exon Offset: {}", current_exon_offset);
                        debug!("Coding Shift % 3: {}", coding_shift % 3);
                        debug!("Shift: {}", (frameshift + current_exon_offset) % 3);
                        if coding_shift % 3 == (frameshift + current_exon_offset) % 3
                            || is_short_exon
                        {
                            // print haplotypes
                            debug!("Should print haplotypes");
                            if !has_frameshift {
                                //follows the main ORF - not including frameshifts!
                                // possible unfinished codon at the end of an exon that continues at the start of the next exon
                                exon_rest = match transcript.strand {
                                    PhasingStrand::Forward => exon.end - (offset + exon_window_len),
                                    PhasingStrand::Reverse => offset - exon.start,
                                };
                                if exon_window_len < 3 {
                                    exon_rest = exon_window_len;
                                }
                            }
                            debug!("Exon Rest {}", exon_rest);
                            let haplotype_results = observations
                                .print_haplotypes(
                                    gene,
                                    transcript,
                                    splice_side_offset,
                                    splice_end,
                                    splice_pos,
                                    splice_gap,
                                    exon.end,
                                    exon.start,
                                    exon_window_len,
                                    refseq,
                                    tsv_writer,
                                    fasta_writer,
                                    is_short_exon,
                                    frameshift,
                                )
                                .unwrap();
                            if haplotype_results.is_empty() {
                                debug!("EMPTY");
                                stopped_frameshift = key;
                                debug!("Offset {}", offset);
                            }
                            if exon_rest < 3 && ((!is_short_exon) || is_first_exon) {
                                prev_hap_vec = haplotype_results;
                            } else {
                                hap_vec = haplotype_results;
                            }
                        }
                    }
                    if frameshift_count == 0 || !main_orf {
                        // no active frameshifts
                        frameshifts.clear();
                        break;
                    }
                    frameshifts.remove(&stopped_frameshift);
                    debug!("frameshifts: {:?}", frameshifts);
                    // if all frames are closed, finish the transcript
                    if frameshifts.is_empty() {
                        break;
                    }
                    // check if the current offset is at a splice side
                    let at_splice_side = match transcript.strand {
                        PhasingStrand::Forward => offset - current_exon_offset == exon.start,
                        PhasingStrand::Reverse => {
                            offset + exon_window_len + current_exon_offset == exon.end
                        }
                    };
                    is_first_exon_window = false;
                    // at a splice side, merge the last sequence of the prev exon and the first sequence of the next exon
                    if at_splice_side && (!is_first_exon) {
                        let first_hap_vec = match transcript.strand {
                            PhasingStrand::Forward => &hap_vec,
                            PhasingStrand::Reverse => &prev_hap_vec,
                        };
                        let sec_hap_vec = match transcript.strand {
                            PhasingStrand::Forward => &prev_hap_vec,
                            PhasingStrand::Reverse => &hap_vec,
                        };
                        debug!("At Splice Side");
                        debug!("{:?}", first_hap_vec);
                        debug!("{:?}", sec_hap_vec);
                        let mut output_map: BTreeMap<(u64, Vec<u8>), (Vec<u8>, IDRecord)> =
                            BTreeMap::new();

                        //Test: Vector for new hap_seq if we are in a short exon
                        let mut new_hap_vec = Vec::new();

                        // iterate over all combinations of splice side haplotypes
                        for hapseq in first_hap_vec {
                            let sequence = &hapseq.sequence;
                            let record = &hapseq.record;
                            for prev_hapseq in sec_hap_vec {
                                let mut prev_sequence = prev_hapseq.sequence.clone();
                                let prev_record = &prev_hapseq.record;
                                // paste together a sequence spanning the splice side
                                prev_sequence.extend(sequence);
                                debug!(
                                    "Complete Sequence : {:?}",
                                    String::from_utf8_lossy(&prev_sequence)
                                );

                                //Test: Keep even wildtype records for merging if we are in a short exon
                                if is_short_exon {
                                    debug!("Exon is shorter than window - merge");
                                    let new_hap_seq = HaplotypeSeq {
                                        sequence: prev_sequence.clone(),
                                        record: prev_record.update(
                                            record,
                                            0,
                                            prev_sequence.to_vec(),
                                        ),
                                    };
                                    new_hap_vec.push(new_hap_seq);
                                    debug!("New HapVec {:?}", new_hap_vec);
                                }
                                // slide window over the spanning sequence
                                let mut splice_offset = 3;
                                if (transcript.strand == PhasingStrand::Reverse) && (exon_rest < 3)
                                {
                                    splice_offset += exon_rest;
                                }
                                let mut end_offset = 3;
                                if is_last_exon_window {
                                    end_offset = 0;
                                }
                                if (prev_sequence.len() as u64) < (2 * window_len) {
                                    // short exon as first or last window:
                                    if transcript.strand == PhasingStrand::Forward {
                                        // take complete first exon
                                        splice_offset = 0;
                                    } else {
                                        //take complete first exon
                                        end_offset = 0;
                                    }
                                }
                                while splice_offset + window_len
                                    <= (prev_sequence.len() - end_offset) as u64
                                {
                                    let out_seq = &prev_sequence[splice_offset as usize
                                        ..(splice_offset + window_len) as usize];
                                    let out_record =
                                        prev_record.update(record, splice_offset, out_seq.to_vec());
                                    debug!(
                                        "Out Sequence : {:?}",
                                        String::from_utf8_lossy(&out_seq)
                                    );
                                    // check if the sequence was already printed at that position, i.e. the variant defining the different haplotypes left the window
                                    let id_tuple = (splice_offset, out_seq.to_vec());
                                    let old_freq = match output_map.get_mut(&id_tuple) {
                                        Some(x) => x.1.freq,
                                        None => 0.0,
                                    };
                                    *output_map
                                        .entry(id_tuple)
                                        .or_insert((out_seq.to_vec(), out_record)) =
                                        (out_seq.to_vec(), out_record.add_freq(old_freq));
                                    splice_offset += 3;
                                }
                            }
                        }
                        if is_short_exon && !is_last_exon {
                            prev_hap_vec = new_hap_vec;
                        } else {
                            for (_key, val) in output_map.iter() {
                                let out_record = &val.1;
                                let out_seq = &val.0;
                                fasta_writer.write(
                                    &out_record.id.to_string(),
                                    None,
                                    &out_seq[..window_len as usize],
                                )?;
                                tsv_writer.serialize(out_record)?;
                            }
                        }
                    }
                    old_offset = splice_side_offset;
                    old_end = splice_end;
                    match transcript.strand {
                        PhasingStrand::Reverse => offset -= 1,
                        PhasingStrand::Forward => offset += 1,
                    }
                    // if all frames are closed, finish the transcript
                    if frameshifts.is_empty() {
                        break;
                    }
                }
                // if all frames are closed, finish the transcript
                if frameshifts.is_empty() {
                    break;
                }
                if is_short_exon {
                    debug!("Exon Rest (End Of Loop): {}", exon_rest);
                    break;
                }
                debug!("Exon Rest (End Of Loop): {}", exon_rest);
            }
        }
        // if all frames are closed, finish the transcript
        if frameshifts.is_empty() {
            continue;
        }
    }
    Ok(())
}

pub fn phase<F: io::Read + io::Seek, G: io::Read, O: io::Write>(
    fasta_reader: &mut fasta::IndexedReader<F>,
    gtf_reader: &mut gff::Reader<G>,
    bcf_reader: bcf::Reader,
    bam_reader: bam::IndexedReader,
    tsv_writer: &mut csv::Writer<fs::File>,
    fasta_writer: &mut fasta::Writer<O>,
    window_len: u64,
    unsupported_allele_warning_only: bool,
) -> Result<(), Box<dyn Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader, false);
    let mut variant_buffer = bcf::buffer::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence
    debug!("refseq length {}", refseq.len());
    debug!("Stared Phasing");
    let mut gene = None;
    let mut start_codon_found = false;
    let mut phase_last_gene = |gene: &Gene| -> Result<(), Box<dyn Error>> {
        if gene.biotype == "protein_coding" {
            phase_gene(
                &gene,
                fasta_reader,
                &mut read_buffer,
                &mut variant_buffer,
                tsv_writer,
                fasta_writer,
                window_len,
                &mut refseq,
                unsupported_allele_warning_only,
            )?;
        }
        Ok(())
    };
    let mut last_chrom: String = "not_yet_set".into();
    let mut last_start: u64 = 0;
    for record in gtf_reader.records() {
        debug!("New Record!");
        let record = record?;
        match record.feature_type() {
            "gene" => {
                // first, phase the last gene
                if let Some(ref g) = gene {
                    phase_last_gene(g)?;
                    last_chrom = g.chrom.to_owned();
                    last_start = g.start();
                }
                debug!("Gene found");
                let gene_name = record
                    .attributes()
                    .get("gene_name")
                    .expect("missing gene_name in GTF");
                if last_chrom == record.seqname() {
                    assert!(last_start <= *record.start(), "Your GTF file is not sorted correctly. Gene {} starts at {}, while previous gene record started at {}.", gene_name, *record.start(), last_start);
                }
                gene = Some(Gene::new(
                    record
                        .attributes()
                        .get("gene_id")
                        .expect("missing gene_id in GTF"),
                    gene_name,
                    record.seqname(),
                    Interval::new(
                        *record.start() as u64 - 1,
                        *record.end() as u64,
                        record.frame(),
                    ),
                    record
                        .attributes()
                        .get("gene_biotype")
                        .expect("missing gene_biotype in GTF"),
                ));
            }
            "transcript" => {
                // register new transcript
                debug!("Transcript found");
                start_codon_found = false;
                gene.as_mut()
                    .expect("no gene record before transcript in GTF")
                    .transcripts
                    .push(Transcript::new(
                        record
                            .attributes()
                            .get("transcript_id")
                            .expect("missing transcript_id attribute in GTF"),
                        record
                            .attributes()
                            .get("transcript_biotype")
                            .expect("missing transcript_biotype in GTF"),
                        PhasingStrand::from(
                            record.strand().expect("missing strand information in GTF"),
                        ),
                    ));
            }
            "CDS" => {
                debug!("CDS found");
                // register exon
                gene.as_mut()
                    .expect("no gene record before exon in GTF")
                    .transcripts
                    .last_mut()
                    .expect("no transcript record before exon in GTF")
                    .exons
                    .push(Interval::new(
                        *record.start() as u64 - 1,
                        *record.end() as u64,
                        record.frame(),
                    ));
            }
            /*             "CDS" => {
                debug!("CDS found");
                // register exon
                gene.as_mut()
                    .expect("no gene record before exon in GTF")
                    .transcripts
                    .last_mut()
                    .expect("no transcript record before exon in GTF")
                    .exons
                    .last_mut()
                    .expect("no exon record before start codon in GTF")
                    .update(*record.start() as u32 -1,
                        *record.end() as u32,
                        record.frame())
                        .unwrap();
            } */
            "start_codon" => {
                if start_codon_found {
                    continue;
                }
                start_codon_found = true;
                if record.strand() == Some(Strand::Forward) {
                    gene.as_mut()
                        .expect("no gene record before start_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before start codon in GTF")
                        .start = *record.start() as u64 - 1;
                } else {
                    gene.as_mut()
                        .expect("no gene record before start_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before start codon in GTF")
                        .end = *record.end() as u64;
                }
            }
            _ => continue,
        }
    }
    if let Some(ref g) = gene {
        phase_last_gene(g)?;
    }

    Ok(())
}
