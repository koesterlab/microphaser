use std::error::Error;
use std::collections::{VecDeque, BTreeMap};
use std::io;
use std::fs;

use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;

use csv;
use sha1;

use itertools::Itertools;

use vec_map::VecMap;

use bio::io::fasta;

use bio::alphabets;

use alphabets::Alphabet;
use alphabets::dna;

fn make_pairs() -> HashMap<&'static [u8], &'static [u8]> {
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
    let map = make_pairs();
    match map.get(v) {
        Some(aa) => Ok(aa),
        None => Err(())
    }
}

fn to_protein(s: &[u8], mut frame: i32) ->  Result<Vec<u8>, ()> {
    let mut r = s.to_vec();
    if frame < 0 {
        r = dna::revcomp(&s.as_bytes().to_vec());
        frame = frame * (-1)
    }
    let mut p = vec![];
    let mut i = frame as usize - 1;
    while i < r.len() - 2 {
        let sub = &r[i..(i + 3)];
        println!("{:?}", String::from_utf8_lossy(sub));
        let aa = to_aminoacid(sub)?;
        p.extend_from_slice(aa);
        i += 3;
    }
    Ok(p)
}

pub fn filter<F: io::Read + io::Seek, O: io::Write>(
    normal_reader: &mut fasta::Reader<F>,
    tumor_reader: &mut fasta::Reader<F>,
    prot_writer: &mut fasta::Writer<O>
) -> Result<(), Box<Error>> {
    // build hashSet from reference
    let mut ref_set = HashSet::new();
    for record in normal_reader.records {
        let id = record.id;
        let seq = record.seq;
        // check if peptide is in forward or reverse orientation
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1
        };

        let pepseq = to_protein(seq, frame)
        ref_set.insert(pepseq)

    }

    for record in tumor_reader.records {
        let id = record.id;
        let seq = record.seq;
         // check if peptide is in forward or reverse orientation
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1
        };

        let neopeptide = to_protein(seq, frame);

        match ref_set.contains(neopeptide) {
            true => (),
            false => prot_writer.write(&format!("{}", record.id), None, &neopeptide)?
        }
    }
}