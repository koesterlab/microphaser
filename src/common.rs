use std::str;
use std::error::Error;

use itertools::any;
use rust_htslib::bcf;



pub enum Variant {
    SNV { pos: u32, alt: u8, is_germline: bool },
    Insertion { pos: u32, seq: Vec<u8>, is_germline: bool },
    Deletion { pos: u32, len: u32, is_germline: bool }
}


impl Variant {
    pub fn new(rec: &bcf::Record) -> Result<Vec<Self>, Box<Error>> {
        let pos = rec.pos();
        let alleles = rec.alleles();
        let refallele = alleles[0];
        let is_germline = any(rec.genotypes()?.get(1).iter().map(|gt_allele| {
            match gt_allele {
                GenotypeAllele::Unphased(i) if i > 0 => true,
                GenotypeAllele::Phased(i) if i > 0 => true,
                _ => false
            }
        }), |is_germline| is_germline);

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
                    seq: (*a).to_owned(),
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
}


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
