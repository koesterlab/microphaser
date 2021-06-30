use std::error::Error;
use std::fs;
use std::io;

use std::collections::{BTreeMap, HashMap, HashSet};

use alphabets::dna;
use bio::alphabets;
use bio::io::fasta;

extern crate bincode;
use bincode::{deserialize_from, serialize_into};

use crate::common::IDRecord;

/* #[derive(Deserialize, debug, Serialize, Clone)]
pub struct IDRecord{
    id: String,
    transcript: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    offset: u32,
    frame: u32,
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
    mutant_sequence: String
} */

fn make_pairs() -> HashMap<&'static [u8], &'static [u8]> {
    // data structure for mapping codons to amino acids
    let grouped = vec![
        ("I", vec!["ATT", "ATC", "ATA"]),
        ("L", vec!["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]),
        ("V", vec!["GTT", "GTC", "GTA", "GTG"]),
        ("F", vec!["TTT", "TTC"]),
        ("M", vec!["ATG"]),
        ("C", vec!["TGT", "TGC"]),
        ("A", vec!["GCT", "GCC", "GCA", "GCG"]),
        ("G", vec!["GGT", "GGC", "GGA", "GGG"]),
        ("P", vec!["CCT", "CCC", "CCA", "CCG"]),
        ("T", vec!["ACT", "ACC", "ACA", "ACG"]),
        ("S", vec!["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
        ("Y", vec!["TAT", "TAC"]),
        ("W", vec!["TGG"]),
        ("Q", vec!["CAA", "CAG"]),
        ("N", vec!["AAT", "AAC"]),
        ("H", vec!["CAT", "CAC"]),
        ("E", vec!["GAA", "GAG"]),
        ("D", vec!["GAT", "GAC"]),
        ("K", vec!["AAA", "AAG"]),
        ("R", vec!["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
        ("X", vec!["TAA", "TAG", "TGA"]),
    ];
    let mut map = HashMap::new();
    for (aa, codons) in grouped.into_iter() {
        for codon in codons {
            map.insert(codon.as_bytes(), aa.as_bytes());
        }
    }
    map
}

fn to_aminoacid(v: &[u8]) -> Result<&'static [u8], ()> {
    // translate a codon triplet to amino acid
    let map = make_pairs();
    match map.get(v) {
        Some(aa) => Ok(aa),
        None => Err(()),
    }
}

fn to_protein(s: &[u8], mut frame: i32) -> Result<Vec<u8>, ()> {
    let case_seq = s.to_ascii_uppercase();
    let mut r = case_seq.to_vec();
    // reverse complement the sequence if the transcript is in reverse orientation
    if frame < 0 {
        r = dna::revcomp(&case_seq.to_vec());
        frame = frame * (-1)
    }
    let mut p = vec![];
    let mut i = frame as usize - 1;
    // iterate over all codons in the sequence and translate them
    while i < r.len() - 2 {
        let sub = &r[i..(i + 3)];
        let aa = to_aminoacid(sub)?;
        p.extend_from_slice(aa);
        i += 3;
    }
    Ok(p)
}

pub fn build<F: io::Read>(
    reference_reader: fasta::Reader<F>,
    binary_writer: fs::File,
    fasta_writer: &mut fasta::Writer<fs::File>,
    peptide_length: usize,
) -> Result<(), Box<dyn Error>> {
    // build hashSet from reference peptide sequences
    let mut ref_set = HashSet::new();
    for record in (reference_reader).records() {
        let record = record?;
        let id = record.id();
        let seq = record.seq();
        // check if peptide is in forward or reverse orientation
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1,
        };
        debug!("{}", String::from_utf8_lossy(seq));
        debug!("{}", String::from_utf8_lossy(&seq.to_ascii_uppercase()));
        let base_length = peptide_length * 3;
        let mut i = 0;
        while i + base_length <= seq.len() {
            let pepseq = to_protein(&seq[i..(i + base_length)], frame).unwrap();
            fasta_writer.write(&format!("{}", id), None, &pepseq)?;
            ref_set.insert(pepseq);
            i += 3;
        }
        // let mut i = 0;
        // while i + peptide_length < pepseq.len() {
        //     let sequence = &pepseq[i..(i+peptide_length)];
        //     ref_set.insert(sequence);
        //     i += 1;
        // }
    }
    // save as binary
    serialize_into(binary_writer, &ref_set)?;
    debug!("Reference is done!");
    Ok(())
}

pub fn compute_ml(alt_depths: &Vec<f64>, depth: &Vec<u32>) -> Result<f64, Box<dyn Error>> {
    let mut ad: f64 = 0.0;
    for d in alt_depths {
        debug!("Alt depth: {}", d);
        ad += d;
    }
    let depth_sum: u32 = depth.iter().sum();
    debug!("Observed alt depth sum: {}", ad);
    debug!("Depth * number of peptides: {}", depth_sum);
    let max_likelihood = ad / depth_sum as f64;
    debug!("ML: {}", max_likelihood);
    Ok(max_likelihood)
}

pub fn filter<F: io::Read, O: io::Write>(
    reference_reader: fs::File,
    tsv_reader: &mut csv::Reader<F>,
    fasta_writer: &mut fasta::Writer<O>,
    normal_writer: &mut fasta::Writer<fs::File>,
    tsv_writer: &mut csv::Writer<fs::File>,
    removed_writer: &mut csv::Writer<fs::File>,
    removed_fasta_writer: &mut fasta::Writer<fs::File>,
    peptide_length: usize,
) -> Result<(), Box<dyn Error>> {
    // load refernce HashSet from file
    let ref_set: HashSet<Vec<u8>> = deserialize_from(reference_reader).unwrap();
    let mut current = (String::from(""), String::from(""), String::from(""));
    let mut current_variant = (String::from(""), String::from(""), String::from(""));
    let mut region_sites = (String::from(""), String::from(""));
    let mut frequencies: BTreeMap<(u64, String, String), Vec<f64>> = BTreeMap::new();
    let mut depth: BTreeMap<(u64, String, String), Vec<u32>> = BTreeMap::new();
    let mut records: BTreeMap<(u64, String, String), Vec<(IDRecord, String, String)>> =
        BTreeMap::new();
    let mut seen_peptides = HashSet::new();
    let mut stop_gained = BTreeMap::new();
    // get peptide info from info.tsv table (including sequences)
    for record in tsv_reader.records() {
        let record = record?;
        debug!("Current Record {:?}", record);
        let row: IDRecord = record.deserialize(None)?;
        let id = &row.id;
        let somatic_positions = &row.somatic_positions;
        // get position of somatic variant
        let som_pos = match somatic_positions.is_empty() {
            // no somatic variant in peptide, but still a neopeptide
            // -> downstream of frameshift, keep complete sequence
            true => 0,
            false => match somatic_positions.contains("|") {
                true => 0,
                false => somatic_positions.parse::<usize>().unwrap(),
            },
        };
        let orientation = *&row.strand.as_str();
        let offset = *&row.offset as usize;
        let mt_seq = &row.mutant_sequence.as_bytes();
        let wt_seq = &row.normal_sequence.as_bytes();
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1,
        };
        let neopeptide = to_protein(mt_seq, frame).unwrap();
        let wt_peptide = match wt_seq.len() == 0 {
            true => vec![],
            false => to_protein(wt_seq, frame).unwrap(),
        };

        let mut i = 0;

        // if Stop Codon in peptide, remove downstream peptides
        let check = (row.transcript.clone(), row.frame);
        //check.push_str(&row.frame.to_string());
        if stop_gained.contains_key(&check) {
            let downstream_of_stop = match orientation {
                "Forward" => offset > *stop_gained.get(&check).unwrap(),
                "Reverse" => offset < *stop_gained.get(&check).unwrap(),
                _ => false,
            };
            if downstream_of_stop {
                continue;
            }
        }
        if neopeptide.contains(&"X".as_bytes()[0]) && ((row.freq - 1.0).abs() < f64::EPSILON || row.frame > 0) {
            let tuple = (row.transcript.clone(), row.frame);
            stop_gained.insert(tuple, offset);
        }
        let current_neo = neopeptide.clone();
        // iterate over the peptide in variable-length windows
        while i + peptide_length <= current_neo.len() {
            let n_peptide = &current_neo[i..(i + peptide_length)];
            // Terminate at stop codon
            if n_peptide.contains(&"X".as_bytes()[0]) {
                break;
            }
            let w_peptide = match wt_peptide.len() >= i + peptide_length {
                true => &wt_peptide[i..(i + peptide_length)],
                false => &wt_peptide,
            };
            debug!(
                "neopeptide {}, wtpeptide {}",
                &String::from_utf8_lossy(n_peptide),
                &String::from_utf8_lossy(w_peptide)
            );
            // skip smaller peptides not containing a somatic variant
            if w_peptide.len() == 0 && som_pos > 0 {
                match orientation {
                    "Forward" => {
                        if ((i + peptide_length) * 3) + offset <= som_pos {
                            i += 1;
                            continue;
                        }
                    }
                    "Reverse" => {
                        if (neopeptide.len() - (i + peptide_length)) * 3 + offset > som_pos {
                            i += 1;
                            continue;
                        }
                    }
                    _ => (),
                };
            }
            i += 1;
            // remove self-similar peptides
            if n_peptide == w_peptide {
                continue;
            }
            // check if we already saw this peptide in this transcript
            let transcript = &row.transcript;
            let sites = &row.variant_sites;
            let current_sites = (transcript.to_string(), sites.to_string());
            let vars = &row.somatic_positions;
            let germline_vars = &row.germline_positions;
            if (
                transcript.to_string(),
                vars.to_string(),
                germline_vars.to_string(),
            ) == current
            {
                debug!("Current: {:?}", current);
                if seen_peptides.contains(&String::from_utf8_lossy(n_peptide).to_string()) {
                    //let row2 = row.clone();
                    //removed_writer.serialize(row2)?;
                    debug!("Already Seen: {}", &String::from_utf8_lossy(n_peptide));
                    continue;
                }
            } else {
                current = (
                    transcript.to_string(),
                    vars.to_string(),
                    germline_vars.to_string(),
                );
                seen_peptides = HashSet::new();
            }
            if current_variant == ("".to_string(), "".to_string(), "".to_string()) {
                current_variant = (
                    transcript.to_string(),
                    vars.to_string(),
                    germline_vars.to_string(),
                );
            }
            debug!("{}", &String::from_utf8_lossy(n_peptide));
            let current_peptide = String::from_utf8_lossy(n_peptide);
            seen_peptides.insert(current_peptide.to_string());
            let mut row2 = row.clone();
            //row2.id = id.replace("F", &i.to_string()).replace("R", &i.to_string());
            let counter_string = format!("{}_", &i.to_string());
            let new_id = counter_string + &row2.id;
            row2.id = new_id;
            let frameshift = row2.frame;
            let current_freq = row2.freq;
            let current_depth = row2.depth;
            let value_tuple = (
                row2,
                String::from_utf8_lossy(n_peptide).to_string(),
                String::from_utf8_lossy(w_peptide).to_string(),
            );
            //let active_variants = (vars.to_string(), germline_vars.to_string());
            if current_sites != region_sites {
                //current != current_variant { //som_pos != current_variant {
                debug!("Printing records");
                for (key, entries) in &records {
                    let ml = compute_ml(&frequencies.get(key).unwrap(), &depth.get(key).unwrap())
                        .unwrap();
                    for (row, np, wp) in entries {
                        let n_peptide = np.as_bytes();
                        let w_peptide = wp.as_bytes();
                        let mut out_row = row.clone();
                        out_row.freq = match out_row.depth == 0 {
                            true => 0.0,
                            false => ml,
                        };
                        debug!("Handling Peptide {}", &String::from_utf8_lossy(n_peptide));
                        // check if the somatic peptide is present in the reference normal peptidome
                        match ref_set.contains(n_peptide) {
                            true => {
                                removed_fasta_writer.write(
                                    &format!("{}", out_row.id),
                                    None,
                                    &n_peptide,
                                )?;
                                removed_writer.serialize(out_row)?;
                                debug!(
                                    "Removed Peptide due to germline similar: {}",
                                    &String::from_utf8_lossy(n_peptide)
                                );
                            }
                            false => {
                                fasta_writer.write(&format!("{}", out_row.id), None, &n_peptide)?;
                                //if we don't have a matching normal, do not write an empty entry to the output
                                if w_peptide.len() > 0 {
                                    normal_writer.write(
                                        &format!("{}", out_row.id),
                                        None,
                                        &w_peptide,
                                    )?;
                                }
                                tsv_writer.serialize(out_row)?;
                            }
                        }
                    }
                }
                frequencies.clear();
                frequencies.insert(
                    (frameshift, vars.to_string(), germline_vars.to_string()),
                    vec![current_freq * current_depth as f64],
                );
                depth.clear();
                depth.insert(
                    (frameshift, vars.to_string(), germline_vars.to_string()),
                    vec![current_depth],
                );
                records.clear();
                records.insert(
                    (frameshift, vars.to_string(), germline_vars.to_string()),
                    vec![value_tuple],
                );
                region_sites = current_sites;
            //records.insert(row2, String::from_utf8_lossy(n_peptide).to_string(), String::from_utf8_lossy(w_peptide).to_string()));
            } else {
                debug!(
                    "Adding to record list {}",
                    &String::from_utf8_lossy(n_peptide)
                );
                //if current_depth > 0 {
                depth
                    .entry((frameshift, vars.to_string(), germline_vars.to_string()))
                    .or_insert(vec![current_depth])
                    .push(current_depth);
                frequencies
                    .entry((frameshift, vars.to_string(), germline_vars.to_string()))
                    .or_insert(vec![current_freq * current_depth as f64])
                    .push(current_freq * current_depth as f64);
                records
                    .entry((frameshift, vars.to_string(), germline_vars.to_string()))
                    .or_insert(vec![value_tuple.clone()])
                    .push(value_tuple);
                //}
                //records.push((row2, String::from_utf8_lossy(n_peptide).to_string(), String::from_utf8_lossy(w_peptide).to_string()));
            }
            // check if the somatic peptide is present in the reference normal peptidome
            // for (row2, n_peptide, w_peptide) in &records {
            //     row2.freq = ml;
            //     match ref_set.contains(n_peptide) {
            //         true => {
            //             removed_fasta_writer.write(&format!("{}", row2.id), None, &n_peptide)?;
            //             removed_writer.serialize(row2)?;
            //             debug!("Removed Peptide due to germline similar: {}", &String::from_utf8_lossy(n_peptide));
            //         },
            //         false => {
            //             fasta_writer.write(&format!("{}", row2.id), None, &n_peptide)?;
            //             //if we don't have a matching normal, do not write an empty entry to the output
            //             if w_peptide.len() > 0 {
            //                 normal_writer.write(&format!("{}", row2.id), None, &w_peptide)?;
            //             }
            //             tsv_writer.serialize(row2)?;
            //         }
            //     }
            // }
        }
    }
    debug!("Printing records");
    for (key, entries) in &records {
        let ml = compute_ml(&frequencies.get(key).unwrap(), &depth.get(key).unwrap()).unwrap();
        for (row, np, wp) in entries {
            let n_peptide = np.as_bytes();
            let w_peptide = wp.as_bytes();
            let mut out_row = row.clone();
            out_row.freq = ml;
            debug!("Handling Peptide {}", &String::from_utf8_lossy(n_peptide));
            // check if the somatic peptide is present in the reference normal peptidome
            match ref_set.contains(n_peptide) {
                true => {
                    removed_fasta_writer.write(&format!("{}", out_row.id), None, &n_peptide)?;
                    removed_writer.serialize(out_row)?;
                    debug!(
                        "Removed Peptide due to germline similar: {}",
                        &String::from_utf8_lossy(n_peptide)
                    );
                }
                false => {
                    fasta_writer.write(&format!("{}", out_row.id), None, &n_peptide)?;
                    //if we don't have a matching normal, do not write an empty entry to the output
                    if w_peptide.len() > 0 {
                        normal_writer.write(&format!("{}", out_row.id), None, &w_peptide)?;
                    }
                    tsv_writer.serialize(out_row)?;
                }
            }
        }
    }
    Ok(())
}
