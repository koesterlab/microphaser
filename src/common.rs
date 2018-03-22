use std::str;
use std::error::Error;

use rust_htslib::bcf;


#[derive(Debug)]
pub enum Variant {
    SNV { pos: u32, alt: u8, is_germline: bool },
    Insertion { pos: u32, seq: Vec<u8>, is_germline: bool },
    Deletion { pos: u32, len: u32, is_germline: bool }
}


impl Variant {
    pub fn new(rec: &mut bcf::Record) -> Result<Vec<Self>, Box<Error>> {
        let is_germline = !rec.info(b"SOMATIC").flag().unwrap_or(false);

        let pos = rec.pos();
        let alleles = rec.alleles();
        let refallele = alleles[0];
        let mut _alleles = Vec::with_capacity(alleles.len() - 1);
        for a in &alleles[1..] {
            if a.len() == 1 && refallele.len() > 1 {
                _alleles.push(
                        Variant::Deletion {
                        pos: pos,
                        len: (refallele.len() - 1) as u32,
                        is_germline: is_germline
                    }
                );
            } else if a.len() > 1 && refallele.len() == 1 {
                _alleles.push(Variant::Insertion {
                    pos: pos,
                    seq: a[1..].to_owned(),
                    is_germline: is_germline
                });
            } else if a.len() == 1 && refallele.len() == 1 {
                _alleles.push(Variant::SNV {
                    pos: pos,
                    alt: a[0],
                    is_germline: is_germline
                });
            } else {
                warn!(
                    "Unsupported variant {} -> {}",
                    str::from_utf8(refallele).unwrap(), str::from_utf8(a).unwrap()
                );
            }
        }

        Ok(_alleles)
    }

    pub fn pos(&self) -> u32 {
        match self {
            &Variant::SNV { pos, .. } => pos,
            &Variant::Deletion { pos, .. } => pos,
            &Variant::Insertion { pos, .. } => pos
        }
    }

    pub fn is_germline(&self) -> bool {
        match self {
            &Variant::SNV { is_germline, .. } => is_germline,
            &Variant::Deletion { is_germline, .. } => is_germline,
            &Variant::Insertion { is_germline, .. } => is_germline
        }
    }

    pub fn frameshift(&self) -> u32 {
        match self {
            &Variant::SNV { .. } => 0,
            &Variant::Deletion { len, .. } => len % 3,
            &Variant::Insertion { ref seq, .. } => seq.len() as u32 % 3
        }
    }
}


#[derive(Debug)]
pub struct Gene {
    pub id: String,
    pub transcripts: Vec<Transcript>,
    pub interval: Interval,
    pub chrom: String
}


impl Gene {
    pub fn new(id: &str, chrom: &str, interval: Interval) -> Self {
        Gene {
            id: id.to_owned(),
            transcripts: Vec::new(),
            chrom: chrom.to_owned(),
            interval: interval
        }
    }

    pub fn start(&self) -> u32 {
        self.interval.start
    }

    pub fn end(&self) -> u32 {
        self.interval.end
    }
}


#[derive(Debug)]
pub struct Transcript {
    pub id: String,
    pub exons: Vec<Interval>
}


impl Transcript {
    pub fn new(id: &str) -> Self {
        Transcript {
            id: id.to_owned(),
            exons: Vec::new()
        }
    }
}


#[derive(Debug)]
pub struct Interval {
    pub start: u32,
    pub end: u32
}


impl Interval {
    pub fn new(start: u32, end: u32) -> Self {
        Interval {
            start: start,
            end: end
        }
    }
}
