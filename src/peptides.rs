use std::error::Error;
use std::io;
use std::fs;

use std::collections::{HashMap, HashSet};

use bio::io::fasta;
use bio::alphabets;
use alphabets::dna;

extern crate bincode;
use bincode::{deserialize_from, serialize_into};

#[derive(Deserialize,Debug,Serialize)]
pub struct IDRecord{
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
    somatic_positions: String,
    somatic_aa_change: String,
    germline_positions: String,
    germline_aa_change: String,
    normal_sequence: String,
    mutant_sequence: String
}

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
        None => Err(())
    }
}

fn to_protein(s: &[u8], mut frame: i32) ->  Result<Vec<u8>, ()> {
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
) -> Result<(), Box<Error>> {
    // build hashSet from reference peptide sequences
    let mut ref_set = HashSet::new();
    for record in (reference_reader).records() {
        let record = record?;
        let id = record.id();
        let seq = record.seq();
        // check if peptide is in forward or reverse orientation
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1
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

pub fn filter<F: io::Read, O: io::Write>(
    reference_reader: fs::File,
    tsv_reader: &mut csv::Reader<F>,
    fasta_writer: &mut fasta::Writer<O>,
    normal_writer: &mut fasta::Writer<fs::File>,
    tsv_writer: &mut csv::Writer<fs::File>
) -> Result<(), Box<Error>> {
    // load refernce HashSet from file
    let ref_set: HashSet<Vec<u8>> = deserialize_from(reference_reader).unwrap();

    // get peptide info from info.tsv table (including sequences)
    for record in tsv_reader.records() {
        let record = record?;
        let row: IDRecord = record.deserialize(None)?;
        let id = &row.id;
        let mt_seq = &row.mutant_sequence.as_bytes();
        let wt_seq = &row.normal_sequence.as_bytes();
        let frame = match id.ends_with("F") {
            true => 1,
            false => -1
        };
        let neopeptide = to_protein(mt_seq, frame).unwrap();
        let wt_peptide = match wt_seq.len() == 0 {
            true => vec![],
            false => to_protein(wt_seq, frame).unwrap()
        };

        // exclude silent mutations
        if neopeptide == wt_peptide {
            continue;
        }

        //check if neopeptide sequence is also found in the reference sequence
        match ref_set.contains(&neopeptide) {
            true => (),
            false => {
                fasta_writer.write(&format!("{}", id), None, &neopeptide)?;
                // if we don't have a matching normal, do not write an empty entry to the output
                if wt_peptide.len() > 0 {
                    normal_writer.write(&format!("{}", id), None, &wt_peptide)?;}
                tsv_writer.serialize(row)?;
            }
        }
    }
    Ok(())
}
