use std::error::Error;
use std::fs;
use std::io;

use std::collections::{HashMap, HashSet};

use alphabets::dna;
use bio::alphabets;
use bio::io::fasta;

extern crate bincode;
use bincode::serialize_into;

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
        let pepseq = to_protein(seq, frame).unwrap();
        ref_set.insert(pepseq);
    }
    // save as binary
    serialize_into(binary_writer, &ref_set)?;
    debug!("Reference is done!");
    Ok(())
}
