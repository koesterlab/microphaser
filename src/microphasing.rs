use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fs;
use std::io;

use csv;
use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;
use bio::io::gff;

use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bcf};

use bio_types::strand::Strand;

use crate::common::{Gene, Variant, Interval, Transcript, PhasingStrand};


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

pub fn switch_ascii_case_vec(v: &Vec<u8>, r: u8) -> Vec<u8> {
    if r.is_ascii_uppercase() {
        v.to_ascii_lowercase()
    } else {
        v.to_ascii_uppercase()
    }
}

pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<dyn Error>> {
    match variant {
        &Variant::SNV { pos, alt, .. } => {
            let b = match read.cigar().read_pos(pos, false, false) {
                Ok(None) => return Ok(false),
                Ok(Some(p)) => read.seq()[p as usize],
                _ => return Ok(false),
            };
            Ok(b == alt)
        }
        &Variant::Insertion { .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                match c {
                    &Cigar::Ins(_) => return Ok(true),
                    _ => (),
                }
            }
            Ok(false)
        }
        &Variant::Deletion { .. } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar().iter() {
                match c {
                    &Cigar::Del(_) => return Ok(true),
                    _ => (),
                }
            }
            Ok(false)
        }
    }
}

#[derive(Debug, Serialize, Clone)]
pub struct IDRecord {
    id: String,
    transcript: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    offset: u32,
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
    normal_sequence: String,
    mutant_sequence: String,
}

impl IDRecord {
    pub fn update(&self, rec: &IDRecord, offset: u32, wt_seq: Vec<u8>, mt_seq: Vec<u8>) -> Self {
        debug!("Start updating record");
        let mut shaid = sha1::Sha1::new();
        // generate unique haplotype ID containing position, transcript and sequence
        let id = format!("{:?}{}{}", &mt_seq, &self.transcript, offset);
        shaid.update(id.as_bytes());
        let fasta_id = format!(
            "{}{}",
            &shaid.digest().to_string()[..15],
            self.strand.chars().next().unwrap()
        );

        let somatic_positions = self.somatic_positions.split("|");
        let somatic_aa_change: Vec<&str> = self.somatic_aa_change.split("|").collect();
        let other_somatic_aa_change: Vec<&str> = rec.somatic_aa_change.split("|").collect();
        let germline_positions = self.germline_positions.split("|");
        let germline_aa_change: Vec<&str> = self.germline_aa_change.split("|").collect();
        let other_germline_aa_change: Vec<&str> = rec.germline_aa_change.split("|").collect();

        let mut s_p_vec = Vec::new();
        let mut g_p_vec = Vec::new();
        let mut s_aa_vec = Vec::new();
        let mut g_aa_vec = Vec::new();

        let mut nvariants = 0;
        let mut nsomatic = 0;

        let mut c = 0;

        for p in somatic_positions {
            debug!("{}", p);
            if p == "" {
                break;
            }
            if self.offset + offset <= p.parse::<u32>().unwrap() {
                s_p_vec.push(p.to_string());
                s_aa_vec.push(somatic_aa_change[c]);
                nsomatic += 1;
                nvariants += 1;
            }
            c += 1;
        }
        c = 0;
        for p in rec.somatic_positions.split("|") {
            debug!("{}", p);
            if p == "" {
                break;
            }
            if rec.offset >= p.parse::<u32>().unwrap() - offset {
                debug!("hey");
                s_p_vec.push(p.to_string());
                s_aa_vec.push(other_somatic_aa_change[c]);
                nsomatic += 1;
                nvariants += 1;
                debug!("{}", nsomatic);
            }
            c += 1
        }
        c = 0;
        for p in germline_positions {
            debug!("{}", p);
            if p == "" {
                break;
            }
            if self.offset + offset <= p.parse::<u32>().unwrap() {
                g_p_vec.push(p.to_string());
                g_aa_vec.push(germline_aa_change[c]);
                nvariants += 1;
            }
            c += 1;
        }
        c = 0;
        for p in rec.germline_positions.split("|") {
            if p == "" {
                break;
            }
            if rec.offset >= p.parse::<u32>().unwrap() - offset {
                g_p_vec.push(p.to_string());
                g_aa_vec.push(other_germline_aa_change[c]);
                nvariants += 1;
            }
            c += 1;
        }

        debug!("nvars {} {}", self.nvar, rec.nvar);
        IDRecord {
            id: fasta_id,
            transcript: self.transcript.to_owned(),
            gene_id: self.gene_id.to_owned(),
            gene_name: self.gene_name.to_owned(),
            chrom: self.chrom.to_owned(),
            offset: offset + self.offset,
            freq: self.freq * rec.freq,
            depth: self.depth,
            nvar: nvariants,
            nsomatic: nsomatic,
            nvariant_sites: self.nvariant_sites + rec.nvariant_sites,
            nsomvariant_sites: self.nsomvariant_sites + rec.nsomvariant_sites,
            strand: self.strand.to_owned(),
            variant_sites: self.variant_sites.to_owned() + &rec.variant_sites,
            somatic_positions: s_p_vec.join("|"),
            somatic_aa_change: s_aa_vec.join("|"),
            germline_positions: g_p_vec.join("|"),
            germline_aa_change: g_aa_vec.join("|"),
            normal_sequence: String::from_utf8(wt_seq).unwrap(),
            mutant_sequence: String::from_utf8(mt_seq).unwrap(),
        }
    }

    pub fn add_freq(&self, freq: f64) -> Self {
        let new_nvar = match self.nvar == 0 {
            true => self.nvar,
            false => match freq {
                f if f > 0.0 => self.nvar - 1,
                _ => self.nvar,
            },
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
            normal_sequence: self.normal_sequence.to_owned(),
            mutant_sequence: self.mutant_sequence.to_owned(),
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
    read: bam::Record,
    haplotype: u64,
}

impl Observation {
    pub fn update_haplotype(&mut self, i: usize, variant: &Variant) -> Result<(), Box<dyn Error>> {
        debug!(
            "Read name {} ; Read pos {} ; variant pos {}",
            String::from_utf8_lossy(self.read.qname()),
            self.read.pos() as u32,
            variant.pos()
        );
        if (self.read.pos() as u32) > variant.pos() {
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
    observations: BTreeMap<u32, Vec<Observation>>,
    variants: VecDeque<Variant>,
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
        for obs in Itertools::flatten(self.observations.values_mut()) {
            obs.haplotype = obs.haplotype & mask;
        }
    }

    /// Add new variants
    pub fn extend_right(&mut self, new_variants: Vec<Variant>) -> Result<(), Box<dyn Error>> {
        let k = new_variants.len();
        debug!("Extend variants!");
        debug!("New variants {}", k);
        if k > 0 {
            for obs in Itertools::flatten(self.observations.values_mut()) {
                obs.haplotype <<= k;
                debug!("{}", String::from_utf8_lossy(obs.read.qname()));
                debug!("{}", obs.haplotype)
            }
        }
        for obs in Itertools::flatten(self.observations.values_mut()) {
            for (i, variant) in new_variants.iter().rev().enumerate() {
                debug!("Checking new variant support in existing Reads");
                obs.update_haplotype(i, variant)?;
            }
        }
        self.variants.extend(new_variants.into_iter());

        Ok(())
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u32, reverse: bool) {
        debug!(
            "Number of reads(before removal): {}",
            self.observations.len()
        );
        let mut observations = self.observations.split_off(&interval_end);
        for obs in Itertools::flatten(observations.values_mut()) {
            debug!("Removed {}", String::from_utf8_lossy(obs.read.qname()));
        }
        for obs in Itertools::flatten(self.observations.values_mut()) {
            debug!("Kept {}", String::from_utf8_lossy(obs.read.qname()));
        }
        if !reverse {
            self.observations = observations; //self.observations.split_off(&interval_end);
        }
        debug!(
            "Number of reads(after removal): {}",
            self.observations.len()
        );
    }

    /// Check if read has already been added to observations
    pub fn contains(&mut self, read: &bam::Record) -> Result<bool, Box<dyn Error>> {
        let pos = read.pos() as u32;
        let qname = read.qname();
        if self.observations.contains_key(&pos) {
            debug!("Read Pos {} already there", pos);
            for obs in self.observations.get(&pos).unwrap() {
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
        read: bam::Record,
        interval_end: u32,
        interval_start: u32,
        reverse: bool,
    ) -> Result<(), Box<dyn Error>> {
        let end_pos = read.cigar().end_pos() as u32;
        let start_pos = read.pos() as u32;
        debug!(
            "Read Start: {}, Read End: {} - Window Start: {}, Window End {}",
            start_pos, end_pos, interval_start, interval_end
        );
        debug!(
            "Read {} to be added",
            String::from_utf8_lossy(read.qname())
        );
        if end_pos >= interval_end && start_pos <= interval_start && !(self.contains(&read).unwrap()) {
            // only insert if end_pos is larger than the interval end
            let mut obs = Observation {
                read: read,
                haplotype: 0,
            };
            for (i, variant) in self.variants.iter().rev().enumerate() {
                debug!("Checking variant support for new reads!");
                obs.update_haplotype(i, variant)?;
            }
            let pos = match reverse {
                true => start_pos,
                false => end_pos,
            };

            self.observations
                .entry(pos)
                .or_insert_with(|| Vec::new())
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
        offset: u32,
        exon_end: u32,
        exon_start: u32,
        window_len: u32,
        refseq: &[u8],
        fasta_writer: &mut fasta::Writer<O>,
        tsv_writer: &mut csv::Writer<fs::File>,
        normal_writer: &mut fasta::Writer<fs::File>,
        is_short_exon: bool
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
        for obs in Itertools::flatten(self.observations.values()) {
            debug!("obs {:?}", obs);
            debug!("obs haplotype:  {}", obs.haplotype);
            *haplotypes.entry(obs.haplotype as usize).or_insert(0) += 1;
        }
        let mut seq = Vec::with_capacity(window_len as usize);
        let mut germline_seq = Vec::with_capacity(window_len as usize);
        debug!("Gene Start {}, Gene End {}", gene.start(), gene.end());
        debug!("Printing at offset: {}", offset);
        debug!("refseq length {}", refseq.len());

        // Strand orientation
        let strand = match transcript.strand {
            PhasingStrand::Reverse => "Reverse",
            PhasingStrand::Forward => "Forward",
        };

        let mut haplotypes_vec = Vec::new();

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
            if variants.is_empty() {
                debug!("HA");
                debug!("Test {}", offset - gene.start());
                germline_seq.extend(
                    &refseq[(offset - gene.start()) as usize
                        ..(offset + window_len - gene.start()) as usize],
                );
                seq.extend(
                    &refseq[(offset - gene.start()) as usize
                        ..(offset + window_len - gene.start()) as usize],
                );
                debug!("HAHA");
            //continue;
        } else {
            while i < window_end {
                debug!("window_end: {}", window_end);
                debug!("i: {}", i);
                debug!("j: {}", j);
                // TODO what happens if a deletion starts upstream of window and overlaps it

                while j < variants.len() && i == variants[j].pos() {
                    debug!("j: {}, variantpos: {}", j, variants[j].pos());
                    let bit_pos = match transcript.strand {
                        PhasingStrand::Reverse => j,
                        PhasingStrand::Forward => variants.len() - 1 -j,
                    };
                    debug!("frameshift: {}", variants[j].frameshift());
                    if bitvector_is_set(haplotype, bit_pos) {
                        debug!("Haplotype: {} ; j: {}", haplotype, j);
                        // if (j + 1) < variants.len() && i == variants[j + 1].pos() {
                        //     j += 1;
                        // }
                        match variants[j] {
                            // if SNV, we push the alternative base instead of the normal one, and change the case of the letter for visualisation
                            &Variant::SNV { alt, .. } => {
                                debug!("Variant: SNV");
                                match variants[j].is_germline() {
                                    true => germline_seq.push(switch_ascii_case(
                                        alt,
                                        refseq[(i - gene.start()) as usize],
                                    )),
                                    false => {
                                        germline_seq.push(refseq[(i - gene.start()) as usize])
                                    }
                                }
                                seq.push(switch_ascii_case(
                                    alt,
                                    refseq[(i - gene.start()) as usize],
                                ));
                                i += 1;
                            }
                            // if insertion, we insert the new bases (with changed case) and decrease the window-end, since we added bases and made the sequence longer
                            &Variant::Insertion { seq: ref s, .. } => {
                                debug!("Variant: INS");
                                match variants[j].is_germline() {
                                    true => germline_seq.extend(
                                        switch_ascii_case_vec(
                                            s,
                                            refseq[(i - gene.start()) as usize],
                                        )
                                        .into_iter(),
                                    ),
                                    false => indel = true,
                                }
                                seq.extend(
                                    switch_ascii_case_vec(
                                        s,
                                        refseq[(i - gene.start()) as usize],
                                    )
                                    .into_iter(),
                                );
                                //indel = true;
                                i += 1;
                                if strand == "Forward" {
                                    window_end -= (s.len() as u32) - 1;
                                } else {
                                    window_end -= s.len() as u32 - variants[j].frameshift() - 1;
                                }
                            }
                            // if deletion, we push the remaining base and increase the index to jump over the deleted bases. Then, we increase the window-end since we lost bases and need to fill up to 27.
                            &Variant::Deletion { len, .. } => {
                                debug!("Variant: DEL");
                                match variants[j].is_germline() {
                                    true => {
                                        germline_seq.push(refseq[(i - gene.start()) as usize])
                                    }
                                    false => indel = true,
                                }
                                seq.push(refseq[(i - gene.start()) as usize]);
                                //indel = true;
                                i += len + 1;
                                // if strand == "Forward" {
                                //     window_end += len;
                                // } else {
                                //     window_end += len - variants[j].frameshift();
                                // }
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

                        j += 1;
                    } else {
                        variant_profile.push(0);
                        j += 1;
                    }
                }
                // if no variant, just push the reference sequence
                if i < window_end {
                    seq.push(refseq[(i - gene.start()) as usize]);
                    germline_seq.push(refseq[(i - gene.start()) as usize]);
                    i += 1;
                }
                debug!("Sequence: {:?}", seq);
            }
            debug!("all variants {}; som variants: {}", n_variants, n_somatic);
        }
            // for indels, do not use the corresponding normal, but search for one with small hamming distance
            if indel {
                debug!("indel");
                germline_seq.clear();
            }
            if strand == "Reverse" {
                debug!("{}", seq.len());
                debug!("reverse");
                // if seq.len() < window_len as usize {
                //     let mut add_seq = refseq[(offset as usize - (window_len as usize - seq.len()) - gene.start() as usize)..(offset - gene.start()) as usize].to_vec();
                //     let mut add_seq_g = add_seq.clone();
                //     add_seq.extend(seq);
                //     seq = add_seq;
                //     add_seq_g.extend(germline_seq);
                //     germline_seq = add_seq_g;
                // }
                if seq.len() > window_len as usize {
                    seq = seq[(seq.len() - window_len as usize)..].to_vec();
                    if germline_seq.len() > window_len as usize {
                        germline_seq = germline_seq[(germline_seq.len() - window_len as usize)..].to_vec();
                    }
                }
            }
            let this_window_len = match seq.len() < window_len as usize {
                true => seq.len() as u32,
                false => window_len
            };
            debug!("germline_seq: {:?} mut_seq: {:?}", germline_seq, seq);
            let mut shaid = sha1::Sha1::new();
            // generate unique haplotype ID containing position, transcript and sequence
            let id = format!("{:?}{}{}", &seq, transcript.id, offset);
            shaid.update(id.as_bytes());
            let fasta_id = format!(
                "{}{}",
                &shaid.digest().to_string()[..15],
                strand.chars().next().unwrap()
            );
            // normal sequence of haplotype
            let normal_peptide = match germline_seq.len() {
                0 => String::from_utf8_lossy(&germline_seq),
                _ => String::from_utf8_lossy(&germline_seq[..this_window_len as usize]),
            };
            // neopeptide sequence
            let neopeptide = String::from_utf8_lossy(&seq[..this_window_len as usize]);
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
                        },
                        // germline
                        1 => {
                            germline_var_pos_vec.push(variants[c as usize].pos().to_string());
                            germline_p_changes_vec.push(variants[c as usize].prot_change());
                        },
                        // not present in this haplotype
                        _ => {}
                    }
                    // check if variant position is already in the variant_site list
                    if c == 0 {
                        n_variantsites += 1;
                        variantsites_pos_vec.push(variants[c as usize].pos().to_string());
                        if !(variants[c as usize].is_germline()) {
                            n_som_variantsites += 1;
                        }
                    }
                    else if !(variants[c as usize].pos() == variants[(c - 1) as usize].pos()) {
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
                offset: offset,
                freq: freq,
                depth: depth,
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
                normal_sequence: normal_peptide.to_owned().to_string(),
                mutant_sequence: neopeptide.to_owned().to_string(),
            };

            // make haplotype record to carry over to next exon
            let rest = match window_end > exon_end {
                true => 0,
                false => exon_end - window_end,
            };
            debug!("Rest: {}", rest);
            let start = offset - exon_start;
            debug!("Start: {}", start);


            //TODO: Activate this part as cleaner code
            let (mutseq, normseq) = match rest < 3 {
                // if there are split-codon bases left at the end of the exon, add them to the sequence
                true => match start < 3 {
                    true => {
                        let mut mutseq = refseq
                            [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
                            .to_vec();
                        mutseq.extend(&seq);
                        mutseq.extend_from_slice(
                            &refseq[(window_end - gene.start()) as usize
                                ..(window_end + rest - gene.start()) as usize],
                        );
                        let mut normseq = refseq
                            [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
                            .to_vec();
                        match indel == false {//normal_peptide.len() > 0 {
                            true => {
                                normseq.extend(&germline_seq);
                                normseq.extend_from_slice(
                                    &refseq[(window_end - gene.start()) as usize
                                        ..(window_end + rest - gene.start()) as usize],
                                );
                            },
                            false => normseq.clear(),
                        }
                        (mutseq, normseq)
                    },
                    false => {
                        let mut mutseq = Vec::new();
                        mutseq.extend(&seq[3 as usize..this_window_len as usize]);
                        mutseq.extend_from_slice(
                            &refseq[(window_end - gene.start()) as usize
                                ..(window_end + rest - gene.start()) as usize],
                        );
                        let mut normseq = Vec::new();
                        if normal_peptide.len() > 0 {
                            normseq.extend(&germline_seq[3 as usize..this_window_len as usize]);
                            normseq.extend_from_slice(
                                &refseq[(window_end - gene.start()) as usize
                                    ..(window_end + rest - gene.start()) as usize],
                            )
                        }
                        (mutseq, normseq)
                    }
                }
                false => match start < 3 {
                     // if there are split-codon bases at the start of an exon, add them to the sequence
                    true => {
                        let mut mutseq = refseq
                            [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
                            .to_vec();
                        mutseq.extend(&seq[..(this_window_len - 3) as usize]);
                        let mut normseq = refseq
                            [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
                            .to_vec();
                        match indel == false {//normal_peptide.len() > 0 {
                            true => {
                                normseq.extend(&germline_seq[..(this_window_len - 3) as usize]);
                            },
                            false => normseq.clear(),
                        }
                        (mutseq, normseq)
                    }
                    false => (Vec::new(), Vec::new())
                }
            };

            // let new_normal = match rest < 3 {
            //     // if there are split-codon bases left at the end of the exon, add them to the sequence
            //     true => {
            //         let mut s = Vec::new();
            //         s.extend(&germline_seq[3 as usize..window_len as usize]);
            //         s.extend_from_slice(&refseq[(window_end - gene.start()) as usize..(window_end + rest - gene.start()) as usize]);
            //         s
            //     }
            //     false => match start < 3 {
            //          // if there are split-codon bases at the start of an exon, add them to the sequence
            //         true => {
            //             let mut s = refseq[(offset - start - gene.start()) as usize..(offset - gene.start()) as usize].to_vec();
            //             s.extend(&germline_seq[..(window_len - 3) as usize]);
            //             s
            //         }
            //         false => Vec::new()
            //     }
            // };

            let hap_seq = HaplotypeSeq {
                sequence: Vec::new(),
                record: IDRecord {
                    id: fasta_id.to_owned(),
                    transcript: transcript.id.to_owned(),
                    gene_id: gene.id.to_owned(),
                    gene_name: gene.name.to_owned(),
                    chrom: gene.chrom.to_owned(),
                    offset: offset,
                    freq: freq,
                    nvar: n_variants,
                    depth: depth,
                    nsomatic: n_somatic,
                    nvariant_sites: n_variantsites as u32,
                    nsomvariant_sites: n_som_variantsites as u32,
                    strand: strand.to_string(),
                    variant_sites: variantsites_pos.to_owned(),
                    somatic_positions: somatic_var_pos.to_owned(),
                    somatic_aa_change: somatic_p_changes.to_owned(),
                    germline_positions: germline_var_pos.to_owned(),
                    germline_aa_change: germline_p_changes.to_owned(),
                    normal_sequence: String::from_utf8_lossy(&normseq).to_string(),
                    mutant_sequence: String::from_utf8_lossy(&mutseq).to_string(),
                },
            };
            // if the exon is so short it will only be used to be combined with other exons, don't split off anything

            // if there are split-codon bases left at the end of the exon, add them to the sequence
            // if rest < 3 {
            //     let mut newseq = Vec::new();
            //     newseq.extend(&seq[3 as usize..window_len as usize]);
            //     newseq.extend_from_slice(
            //         &refseq[(window_end - gene.start()) as usize
            //             ..(window_end + rest - gene.start()) as usize],
            //     );
            //     let mut new_normal = Vec::new();
            //     if normal_peptide.len() > 0 {
            //         new_normal.extend(&germline_seq[3 as usize..window_len as usize]);
            //         new_normal.extend_from_slice(
            //             &refseq[(window_end - gene.start()) as usize
            //                 ..(window_end + rest - gene.start()) as usize],
            //         )
            //     }
            //     hap_seq = HaplotypeSeq {
            //         sequence: Vec::new(),
            //         record: IDRecord {
            //             id: fasta_id.to_owned(),
            //             transcript: transcript.id.to_owned(),
            //             gene_id: gene.id.to_owned(),
            //             gene_name: gene.name.to_owned(),
            //             chrom: gene.chrom.to_owned(),
            //             offset: offset,
            //             freq: freq,
            //             depth: depth,
            //             nvar: n_variants,
            //             nsomatic: n_somatic,
            //             nvariant_sites: n_variantsites as u32,
            //             nsomvariant_sites: n_som_variantsites as u32,
            //             strand: strand.to_string(),
            //             variant_sites: variantsites_pos.to_owned(),
            //             somatic_positions: somatic_var_pos.to_owned(),
            //             somatic_aa_change: somatic_p_changes.to_owned(),
            //             germline_positions: germline_var_pos.to_owned(),
            //             germline_aa_change: germline_p_changes.to_owned(),
            //             normal_sequence: String::from_utf8_lossy(&new_normal).to_string(),
            //             mutant_sequence: String::from_utf8_lossy(&newseq).to_string(),
            //         },
            //     };
            // }
            // // if there are split-codon bases at the start of an exon, add them to the sequence
            // if start < 3 {
            //     let mut newseq = refseq
            //         [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
            //         .to_vec();
            //     newseq.extend(&seq[..(window_len - 3) as usize]);
            //     let mut new_normal = refseq
            //         [(offset - start - gene.start()) as usize..(offset - gene.start()) as usize]
            //         .to_vec();
            //     match normal_peptide.len() > 0 {
            //         true => new_normal.extend(&germline_seq[..(window_len - 3) as usize]),
            //         false => new_normal.clear(),
            //     }
            //     debug!("Starting_sequence: {:?}", String::from_utf8_lossy(&seq));
            //     hap_seq = HaplotypeSeq {
            //         sequence: Vec::new(),
            //         record: IDRecord {
            //             id: fasta_id.to_owned(),
            //             transcript: transcript.id.to_owned(),
            //             gene_id: gene.id.to_owned(),
            //             gene_name: gene.name.to_owned(),
            //             chrom: gene.chrom.to_owned(),
            //             offset: offset,
            //             freq: freq,
            //             depth: depth,
            //             nvar: n_variants,
            //             nsomatic: n_somatic,
            //             nvariant_sites: n_variantsites as u32,
            //             nsomvariant_sites: n_som_variantsites as u32,
            //             strand: strand.to_string(),
            //             variant_sites: variantsites_pos.to_owned(),
            //             somatic_positions: somatic_var_pos.to_owned(),
            //             somatic_aa_change: somatic_p_changes.to_owned(),
            //             germline_positions: germline_var_pos.to_owned(),
            //             germline_aa_change: germline_p_changes.to_owned(),
            //             normal_sequence: String::from_utf8_lossy(&new_normal).to_string(),
            //             mutant_sequence: String::from_utf8_lossy(&newseq).to_string(),
            //         },
            //     };
            // }

            haplotypes_vec.push(hap_seq);

            // write neopeptides, information and matching normal peptide to files
            if record.nsomatic > 0 && !(is_short_exon){
                fasta_writer.write(&format!("{}", record.id), None, &seq[..this_window_len as usize])?;
                if germline_seq.len() > 0 {
                    normal_writer.write(
                        &format!("{}", record.id),
                        None,
                        &germline_seq[..this_window_len as usize],
                    )?;
                }
                tsv_writer.serialize(record)?;
            }
        }
        debug!("{:?}", haplotypes_vec);
        Ok(haplotypes_vec)
    }
}

pub fn phase_gene<F: io::Read + io::Seek, O: io::Write>(
    gene: &Gene,
    fasta_reader: &mut fasta::IndexedReader<F>,
    read_buffer: &mut bam::RecordBuffer,
    variant_buffer: &mut bcf::buffer::RecordBuffer,
    fasta_writer: &mut fasta::Writer<O>,
    tsv_writer: &mut csv::Writer<fs::File>,
    normal_writer: &mut fasta::Writer<fs::File>,
    window_len: u32,
    refseq: &mut Vec<u8>,
) -> Result<(), Box<dyn Error>> {
    // if an exon is near to the gene end, a deletion could cause refseq to overflow, so we increase the length of refseq
    let end_overflow = 100;
    fasta_reader.fetch(
        &gene.chrom,
        gene.start() as u64,
        (gene.end() + end_overflow) as u64
    )?;
    fasta_reader.read(refseq)?;
    let mut variant_tree = BTreeMap::new();
    let mut read_tree = BTreeMap::new();
    debug!("Start Phasing");
    read_buffer.fetch(&gene.chrom.as_bytes(), gene.start(), gene.end())?;

    let mut max_read_len = 0 as u32;
    // load read buffer into BTree
    for rec in read_buffer.iter() {
        if rec.seq().len() as u32 > max_read_len {
            max_read_len = rec.seq().len() as u32;
        }
        read_tree
            .entry(rec.pos() as u32)
            .or_insert_with(|| Vec::new())
            .push(rec.clone())
    }

    // if max_read_len is zero, there are no reads - break
    // if max_read_len == 0 {
    //     return Ok(())
    // } TODO: Think about ending here

    // if the fixed window_len is to large, set it to read length

    debug!("Reads Tree length: {}", read_tree.len());

    // load variant buffer into BTree
    let (_addedvars, _deletedvars) =
        variant_buffer.fetch(&gene.chrom.as_bytes(), gene.start(), gene.end())?;
    let _vars = variant_buffer
        .iter_mut()
        .map(|rec| variant_tree.insert(rec.pos(), Variant::new(rec).unwrap()))
        .collect_vec();


    for transcript in &gene.transcripts {
        debug!("is coding: {}", transcript.is_coding());
        if !(transcript.is_coding()) {
            continue;
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
        // Variable showing exon count
        let mut exon_number = 0;
        for exon in &transcript.exons {
            debug!("Exon Start: {}", exon.start);
            debug!("Exon End: {}", exon.end);
            if exon.start > exon.end {
                continue;
            }
            exon_number += 1;
            debug!("Exon Length: {}", exon.end - exon.start);
            // Possible offset at the exon start, first nucleotides could be part of a codon started in the previous exon
            let exon_len = exon.end - exon.start;
            debug!("Exon Rest: {}", exon_rest);
            let current_exon_offset = match exon_number { 
                1 => exon.frame, //u32::from_str(exon.frame).unwrap(),
                _ => match exon_rest {
                    0 => 0,
                    _ => 3 - exon_rest,
                },
            };
            debug!("Exon Offset: {}", current_exon_offset);
            let is_short_exon = window_len > exon_len - current_exon_offset;
            // if the exon is shorter than the window, we need to fix the window len for this exon
            let mut exon_window_len = match is_short_exon {
                false => window_len,
                true => (exon_len - current_exon_offset) - ((exon_len - current_exon_offset) % 3),
            };
            if exon_window_len == 0 {
                exon_window_len = exon_len
            }
            exon_rest = 0;
            let mut offset = if transcript.strand == PhasingStrand::Reverse {
                exon.end - exon_window_len - current_exon_offset
            } else {
                exon.start + current_exon_offset
            };
            debug!("Exon window lenght: {}", exon_window_len);
            debug!("Starting Offset of the Exon: {}", offset);
            let mut old_offset = offset;
            debug!("Variants left from previous Exon: {}", last_window_vars);
            observations.shrink_left(last_window_vars);
            last_window_vars = 0;
            loop {
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
                debug!("Offset {}, old offset {}", offset, old_offset);
                // advance window to next position
                let nvars = Itertools::flatten(
                    variant_tree
                        .range(offset..(offset + exon_window_len))
                        .map(|var| var.1),
                )
                .count();
                // store number of variants in window in case it is the last window for this exon
                last_window_vars = nvars;
                debug!("Variants in window: {}", nvars);
                // first window in the exon, all variants found are newly added
                let added_vars = if offset == old_offset {
                    nvars
                // if we advance the window (forward or reverse), just the newly added variants are counted
                // forward orientation
                } else if offset > old_offset {
                    Itertools::flatten(
                        variant_tree
                            .range((old_offset + exon_window_len)..(offset + exon_window_len))
                            .map(|var| var.1),
                    )
                    .count()
                // reverse orientation
                } else {
                    Itertools::flatten(variant_tree.range(offset..old_offset).map(|var| var.1))
                        .count()
                };

                // first window in the exon, no variants are deleted
                let deleted_vars = if offset == old_offset {
                    0
                // if we advance the window (forward or reverse), we will delete all variants that drop out of the window bounds
                // forward orientation
                } else if offset > old_offset {
                    Itertools::flatten(variant_tree.range(old_offset..offset).map(|var| var.1))
                        .count()
                // reverse orientation
                } else {
                    Itertools::flatten(
                        variant_tree
                            .range((offset + exon_window_len)..(old_offset + exon_window_len))
                            .map(|var| var.1),
                    )
                    .count()
                };
                debug!("Offset + wlen: {}", offset + exon_window_len);
                debug!("Old offset + wlen: {}", old_offset + exon_window_len);
                debug!(
                    "Offset: {} - max_read_len {} - window_len {}",
                    offset,
                    max_read_len,
                    exon_window_len
                );
                let reads = if transcript.strand == PhasingStrand::Reverse {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    debug!("Offset: {} ; (Offset - (max_read_len - window_len)) = {}", offset, offset - (max_read_len - exon_window_len));
                    if offset == exon.end - exon_window_len - current_exon_offset {
                        debug!("First exon window");
                        Itertools::flatten(
                            read_tree
                                .range((offset - (max_read_len - exon_window_len))..(offset + 1))
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                    // while advancing the window (reverse orientation), we only add reads that end in the range between old and new window end, so we don't count any read twice
                    else {
                        Itertools::flatten(
                            read_tree
                                .range(
                                    (offset - (max_read_len - exon_window_len))
                                        ..(offset + 1),
                                )
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                } else {
                    // at the first window of the exon, we add all reads (including those starting before the window start) that enclose the window
                    if offset == exon.start + current_exon_offset {
                        Itertools::flatten(
                            read_tree
                                .range((offset - (max_read_len - exon_window_len))..(offset + 1))
                                .map(|rec| rec.1),
                        )
                        .collect_vec()
                    }
                    // while advancing the window, we only add reads that start in the range between old and new window, so we don't count any read twice
                    else {
                        Itertools::flatten(
                            read_tree.range((offset)..(offset + 1)).map(|rec| rec.1),
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
                        observations.cleanup_reads(offset, reverse);
                    } else {
                        observations.cleanup_reads(offset + exon_window_len, reverse);
                    }
                    // delete columns
                    observations.shrink_left(deleted_vars);

                    // add new reads
                    debug!("Reads: {}", reads.len());
                    for read in reads {
                        observations.push_read(
                            read.clone(),
                            offset + exon_window_len,
                            offset,
                            reverse,
                        )?;
                    }

                    // collect variants
                    let variants = match transcript.strand {
                        PhasingStrand::Reverse => Itertools::flatten(
                            variant_tree
                                .range_mut(offset..(offset + exon_window_len))
                                .rev()
                                .map(|var| var.1.clone()),
                        )
                        .skip(nvars - added_vars)
                        .collect_vec(),
                        PhasingStrand::Forward => Itertools::flatten(
                            variant_tree
                                .range_mut(offset..(offset + exon_window_len))
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
                            PhasingStrand::Reverse => exon.end - offset,
                        };
                        debug!("Offset - Exonstart % 3: {}", coding_shift % 3);
                        debug!("Shift: {}", frameshift + current_exon_offset);
                        if coding_shift % 3 == frameshift + current_exon_offset {
                            // print haplotypes
                            debug!("Should print haplotypes");
                            // possible unfinished codon at the end of an exon that continues at the start of the next exon
                            exon_rest = match transcript.strand {
                                PhasingStrand::Forward => exon.end - (offset + exon_window_len),
                                PhasingStrand::Reverse => offset - exon.start,
                            };
                            debug!("Exon Rest {}", exon_rest);
                            if exon_rest < 3 && (!is_short_exon) {
                                prev_hap_vec = observations
                                    .print_haplotypes(
                                        gene,
                                        transcript,
                                        offset,
                                        exon.end,
                                        exon.start,
                                        exon_window_len,
                                        refseq,
                                        fasta_writer,
                                        tsv_writer,
                                        normal_writer,
                                        is_short_exon
                                    )
                                    .unwrap();
                            } else {
                                hap_vec = observations
                                    .print_haplotypes(
                                        gene,
                                        transcript,
                                        offset,
                                        exon.end,
                                        exon.start,
                                        exon_window_len,
                                        refseq,
                                        fasta_writer,
                                        tsv_writer,
                                        normal_writer,
                                        is_short_exon
                                    )
                                    .unwrap();
                            }
                        }
                    }

                    // check if the current offset is at a splice side
                    debug!("{}", offset - current_exon_offset);
                    let at_splice_side = match transcript.strand {
                        PhasingStrand::Forward => offset - current_exon_offset == exon.start,
                        PhasingStrand::Reverse => {
                            offset + exon_window_len + current_exon_offset == exon.end
                        }
                    };

                    // at a splice side, merge the last sequence of the prev exon and the first sequence of the next exon
                    if at_splice_side {
                        debug!("SpliceSide");
                        let first_hap_vec = match transcript.strand {
                            PhasingStrand::Forward => &hap_vec,
                            PhasingStrand::Reverse => &prev_hap_vec,
                        };
                        let sec_hap_vec = match transcript.strand {
                            PhasingStrand::Forward => &prev_hap_vec,
                            PhasingStrand::Reverse => &hap_vec,
                        };
                        debug!("first_hap_vec: {:?}", first_hap_vec);
                        debug!("sec_hap_vec: {:?}", sec_hap_vec);

                        let mut output_map: BTreeMap<
                            (u32, Vec<u8>, Vec<u8>),
                            (Vec<u8>, IDRecord, Vec<u8>),
                        > = BTreeMap::new();

                        //Test: Vector for new hap_seq if we are in a short exon
                        let mut new_hap_vec = Vec::new();

                        // iterate over all combinations of splice side haplotypes
                        for hapseq in first_hap_vec {
                            //let sequence = &hapseq.sequence;
                            let record = &hapseq.record;
                            let wt_sequence = &record.normal_sequence;
                            let mt_sequence = &record.mutant_sequence;
                            for prev_hapseq in sec_hap_vec {
                                let prev_record = &prev_hapseq.record;
                                let prev_wt_sequence = &prev_record.normal_sequence;
                                let prev_mt_sequence = &prev_record.mutant_sequence;
                                // combine the normal sequence
                                let new_wt = format!("{}{}", prev_wt_sequence, wt_sequence);
                                let new_wt_sequence = new_wt.as_bytes();
                                // combine all possibilites of mutated sequences
                                let mut new_mt_sequences = Vec::new();
                                if wt_sequence != mt_sequence {
                                    new_mt_sequences
                                        .push(format!("{}{}", prev_wt_sequence, mt_sequence));
                                    if prev_wt_sequence != prev_mt_sequence {
                                        new_mt_sequences
                                            .push(format!("{}{}", prev_mt_sequence, wt_sequence));
                                        new_mt_sequences
                                            .push(format!("{}{}", prev_mt_sequence, mt_sequence));
                                    }
                                } else if prev_wt_sequence != prev_mt_sequence {
                                    new_mt_sequences
                                        .push(format!("{}{}", prev_mt_sequence, mt_sequence));
                                }
                                debug!(
                                    "Complete WT Sequence : {:?}",
                                    String::from_utf8_lossy(&new_wt_sequence)
                                );

                                //Test: Keep even wildtype records for merging if we are in a short exon
                                if is_short_exon {
                                    debug!("Exon is shorter than window - merge");
                                    let new_hap_seq =  HaplotypeSeq {
                                        sequence: Vec::new(),
                                        record: prev_record.update(
                                            record,
                                            0,
                                            new_wt_sequence.to_vec(),
                                            new_wt_sequence.to_vec()
                                        )
                                    };
                                    new_hap_vec.push(new_hap_seq);
                                    debug!("New HapVec {:?}", new_hap_vec );
                                }

                                // slide window over the spanning sequence
                                for new_mt in new_mt_sequences {
                                    let new_mt_sequence = new_mt.as_bytes();
                                    debug!(
                                        "Complete MT Sequence : {:?}",
                                        String::from_utf8_lossy(&new_mt_sequence)
                                    );

                                    //Test: Merge short exon to previous window and save as hap_vec
                                    if is_short_exon {
                                        debug!("Exon is shorter than window - merge");
                                        let new_hap_seq =  HaplotypeSeq {
                                            sequence: Vec::new(),
                                            record: prev_record.update(
                                                record,
                                                0,
                                                new_mt_sequence.to_vec(),
                                                new_wt_sequence.to_vec()
                                            )
                                        };
                                        new_hap_vec.push(new_hap_seq);
                                        debug!("New HapVec {:?}", new_hap_vec );
                                        continue;
                                    }

                                    let mut splice_offset = 0;
                                    debug!("MT_Seq len {}", new_mt_sequence.len() as u32);
                                    debug!("WT_Seq len {}", new_wt_sequence.len() as u32);
                                    while splice_offset + window_len <= new_mt_sequence.len() as u32
                                    {
                                        debug!("splice offset: {}", splice_offset);
                                        debug!("splice offset + windowlen: {}", splice_offset + window_len);
                                        // check if wildtype sequence is shorter because of indels in on of the sequences
                                        let out_wt_seq = match splice_offset + window_len <= new_wt_sequence.len() as u32 {
                                            true => &new_wt_sequence[splice_offset as usize
                                                ..(splice_offset + window_len) as usize],
                                            false => &[]
                                        };
                                        let out_mt_seq = &new_mt_sequence[splice_offset as usize
                                            ..(splice_offset + window_len) as usize];

                                        debug!(
                                            "Out MT Sequence : {:?}",
                                            String::from_utf8_lossy(&out_mt_seq)
                                        );
                                        debug!(
                                            "Out WT Sequence : {:?}",
                                            String::from_utf8_lossy(&out_wt_seq)
                                        );
                                        // non mutated sites
                                        if out_wt_seq == out_mt_seq {
                                            debug!("equal");
                                            splice_offset += 3;
                                            continue;
                                        }
                                        let out_record = prev_record.update(
                                            record,
                                            splice_offset + 3,
                                            out_wt_seq.to_vec(),
                                            out_mt_seq.to_vec(),
                                        );
                                        debug!("prevRecord: {:?}", prev_record);
                                        debug!("afterRecord: {:?}", record);
                                        debug!("Record: {:?}", out_record);
                                        // check if the sequence was already printed at that position, i.e. the variant defining the different haplotypes left the window
                                        let id_tuple = (
                                            splice_offset,
                                            out_mt_seq.to_vec(),
                                            out_wt_seq.to_vec(),
                                        );
                                        let old_freq = match output_map.get_mut(&id_tuple) {
                                            Some(x) => x.1.freq,
                                            None => 0.0,
                                        };
                                        *output_map.entry(id_tuple).or_insert((
                                            out_mt_seq.to_vec(),
                                            out_record,
                                            out_wt_seq.to_vec(),
                                        )) = (
                                            out_mt_seq.to_vec(),
                                            out_record.add_freq(old_freq),
                                            out_wt_seq.to_vec(),
                                        );
                                        splice_offset += 3;
                                    }
                                }
                            }
                        }
                        if is_short_exon {
                            prev_hap_vec = new_hap_vec;
                        }
                        else {
                            for (_key, val) in output_map.iter() {
                                let out_record = &val.1;
                                let out_mt_seq = &val.0;
                                let out_wt_seq = &val.2;

                                debug!("Out Sequence 2: {:?}", String::from_utf8_lossy(&out_mt_seq));
                                // filter relevant haplotypes : variant haplotypes and their wildtype counterparts
                                if out_record.nsomatic > 0 {
                                    fasta_writer.write(
                                        &format!("{}", out_record.id),
                                        None,
                                        &out_mt_seq[..window_len as usize],
                                    )?;
                                    if out_wt_seq != &[] {
                                        normal_writer.write(
                                            &format!("{}", out_record.id),
                                            None,
                                            &out_wt_seq[..window_len as usize],
                                        )?;
                                    }
                                    tsv_writer.serialize(out_record)?;
                                }
                            }
                            if is_short_exon {
                                prev_hap_vec = new_hap_vec;
                            }
                        }
                        debug!("Prevhapvec {:?}", prev_hap_vec );
                    }
                    old_offset = offset;
                    match transcript.strand {
                        PhasingStrand::Reverse => offset -= 1,
                        PhasingStrand::Forward => offset += 1,
                    }
                    debug!("New offset: {}", offset);
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
    normal_writer: &mut fasta::Writer<fs::File>,
    window_len: u32,
) -> Result<(), Box<dyn Error>> {
    let mut read_buffer = bam::RecordBuffer::new(bam_reader, false);
    let mut variant_buffer = bcf::buffer::RecordBuffer::new(bcf_reader);
    let mut refseq = Vec::new(); // buffer for reference sequence
    debug!("refseq length {}", refseq.len());

    let mut gene = None;
    let mut phase_last_gene = |gene: Option<Gene>| -> Result<(), Box<dyn Error>> {
        if let Some(ref gene) = gene {
            if gene.biotype == "protein_coding" {
                phase_gene(
                    &gene,
                    fasta_reader,
                    &mut read_buffer,
                    &mut variant_buffer,
                    fasta_writer,
                    tsv_writer,
                    normal_writer,
                    window_len,
                    &mut refseq,
                )?;
            }
        }
        Ok(())
    };
    for record in gtf_reader.records() {
        debug!("New Record!");
        let record = record?;
        match record.feature_type() {
            "gene" => {
                // first, phase the last gene
                phase_last_gene(gene)?;
                debug!("Gene found");
                gene = Some(Gene::new(
                    record
                        .attributes()
                        .get("gene_id")
                        .expect("missing gene_id in GTF"),
                    record
                        .attributes()
                        .get("gene_name")
                        .expect("missing gene_name in GTF"),
                    record.seqname(),
                    Interval::new(*record.start() as u32 - 1, *record.end() as u32, record.frame()),
                    record
                        .attributes()
                        .get("gene_biotype")
                        .expect("missing gene_biotype in GTF"),
                ));
            }
            "transcript" => {
                // register new transcript
                debug!("Transcript found");
                gene.as_mut()
                    .expect("no gene record before transcript in GTF")
                    .transcripts
                    .push(Transcript::new(
                        record
                            .attributes()
                            .get("transcript_id")
                            .expect("missing transcript_id attribute in GTF"),
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
                        *record.start() as u32 - 1,
                        *record.end() as u32,
                        record.frame()
                    ));
            }
            "start_codon" => {
                if record.strand() == Some(Strand::Forward) {
                    gene.as_mut()
                        .expect("no gene record before start_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before start codon in GTF")
                        .start = *record.start() as u32 - 1;
                } else {
                    gene.as_mut()
                        .expect("no gene record before start_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before start codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before start codon in GTF")
                        .end = *record.end() as u32;
                }
            }
            "stop_codon" => {
                if record.strand() == Some(Strand::Forward) {
                    gene.as_mut()
                        .expect("no gene record before stop_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before stop codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before stop codon in GTF")
                        .end = *record.end() as u32;
                } else {
                    debug!("stop_codon_start {}", *record.start() as u32 - 1);
                    gene.as_mut()
                        .expect("no gene record before stop_codon in GTF")
                        .transcripts
                        .last_mut()
                        .expect("no transcript record before stop codon in GTF")
                        .exons
                        .last_mut()
                        .expect("no exon record before stop codon in GTF")
                        .start = *record.start() as u32 - 1;
                }
            }
            _ => continue,
        }
    }
    phase_last_gene(gene)?;

    Ok(())
}
