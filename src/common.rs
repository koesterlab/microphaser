use std::str;
use std::error::Error;
use std::cmp::Ordering;


use bio::utils::Strand;
use rust_htslib::bcf;

use std::borrow::ToOwned;


#[derive(Debug)]
pub struct Annotation {
    pub prot_change: String
}

impl Annotation {
    pub fn new(rec: &mut bcf::Record) -> Self {
        //let info = str::from_utf8(rec.info(b"ANN").string().unwrap().unwrap_or(Vec::new())[0]).unwrap_or("");
        let info = match rec.info(b"ANN").string() {
            Err(_e) => "",
            Ok(v) => str::from_utf8(v.unwrap_or(Vec::new())[0]).unwrap(),
        };
        let pc = match info {
            "" => "".to_string(),
            _ => info.split('|').nth(10).unwrap().to_string()
        };//info.split('|').nth(10).unwrap().to_string();
        Annotation {
            prot_change: pc
        }
    }
}


#[derive(Debug, Clone)]
pub enum Variant {
    SNV { pos: u32, alt: u8, is_germline: bool, prot_change: String },
    Insertion { pos: u32, seq: Vec<u8>, is_germline: bool, prot_change: String },
    Deletion { pos: u32, len: u32, is_germline: bool, prot_change: String }
}


impl Variant {
    pub fn new(rec: &mut bcf::Record) -> Result<Vec<Self>, Box<Error>> {
        let is_germline = !rec.info(b"SOMATIC").flag().unwrap_or(false);

        let ann = Annotation::new(rec);

        let prot_change = ann.prot_change.as_str();

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
                        is_germline: is_germline,
                        prot_change: prot_change.to_owned()
                    }
                );
            } else if a.len() > 1 && refallele.len() == 1 {
                _alleles.push(Variant::Insertion {
                    pos: pos,
                    seq: a[0..].to_owned(),
                    is_germline: is_germline,
                    prot_change: prot_change.to_owned()
                });
            } else if a.len() == 1 && refallele.len() == 1 {
                _alleles.push(Variant::SNV {
                    pos: pos,
                    alt: a[0],
                    is_germline: is_germline,
                    prot_change: prot_change.to_owned()
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

    pub fn end_pos(&self) -> u32 {
        match self {
            &Variant::SNV { pos, .. } => pos,
            &Variant::Deletion { pos, len, .. } => pos + len - 1,
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

    pub fn prot_change(&self) -> String {
        match self {
            &Variant::SNV { ref prot_change, .. } => prot_change.to_owned(),
            &Variant::Deletion { ref prot_change, .. } => prot_change.to_owned(),
            &Variant::Insertion { ref prot_change, .. } => prot_change.to_owned()
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
    pub name: String,
    pub transcripts: Vec<Transcript>,
    pub interval: Interval,
    pub chrom: String,
    pub biotype: String
}


impl Gene {
    pub fn new(id: &str, name: &str, chrom: &str, interval: Interval, biotype: &str) -> Self {
        Gene {
            id: id.to_owned(),
            name: name.to_owned(),
            transcripts: Vec::new(),
            chrom: chrom.to_owned(),
            interval: interval,
            biotype: biotype.to_owned()
        }
    }

    pub fn start(&self) -> u32 {
        self.interval.start
    }

    pub fn end(&self) -> u32 {
        self.interval.end
    }
}

#[derive(Debug, PartialEq)]
pub enum PhasingStrand {
    Forward,
    Reverse,
}

impl From<Strand> for PhasingStrand{
    fn from(strand: Strand) -> Self {
        let s = match strand {
            Strand::Forward => PhasingStrand::Forward,
            Strand::Reverse => PhasingStrand::Reverse,
            _ => panic!(
                "Unsupported Strand orientation! Only Forward (+) and Reverse(-) allowed"
            ),
        };
        s
    }
}

#[derive(Debug)]
pub struct Transcript {
    pub id: String,
    pub strand: PhasingStrand,
    pub exons: Vec<Interval>,
}


impl Transcript {
    pub fn new(id: &str, strand: PhasingStrand) -> Self {
        Transcript {
            id: id.to_owned(),
            strand: strand,
            exons: Vec::new(),
        }
    }
    pub fn is_coding(&self) -> bool{
        return !(self.exons.is_empty());
    }
}


#[derive(Debug)]
pub struct Interval {
    pub start: u32,
    pub end: u32
}

impl Ord for Interval {
    fn cmp(&self, other: &Interval) -> Ordering {
        self.start.cmp(&other.start)
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Interval) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Interval {
    fn eq(&self, other: &Interval) -> bool {
        self.start == other.start
    }
}

impl Eq for Interval {}

impl Clone for Interval {
    fn clone(&self) -> Interval {
    *self
    }
}

impl Copy for Interval {}

impl Interval {
    pub fn new(start: u32, end: u32) -> Self {
        Interval {
            start: start,
            end: end
        }
    }
}




