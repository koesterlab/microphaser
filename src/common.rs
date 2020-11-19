use std::cmp::Ordering;

use bio_types::strand::Strand;

use std::error::Error;
use std::str;
use std::str::FromStr;

use rust_htslib::bcf;

use std::borrow::ToOwned;

#[derive(Debug)]
pub struct Annotation {
    pub prot_change: String,
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
            _ => match info.split('|').into_iter().position(|e| e.contains("p.")) {
                Some(index) => info.split('|').nth(index).unwrap().to_string(),
                _ => "".to_string(),
            },
        }; //info.split('|').nth(10).unwrap().to_string();
        Annotation { prot_change: pc }
    }
}

#[derive(Debug, Clone)]
pub enum Variant {
    SNV {
        pos: u32,
        alt: u8,
        is_germline: bool,
        prot_change: String,
    },
    Insertion {
        pos: u32,
        seq: Vec<u8>,
        len: u32,
        is_germline: bool,
        prot_change: String,
    },
    Deletion {
        pos: u32,
        len: u32,
        is_germline: bool,
        prot_change: String,
    },
}

impl Variant {
    pub fn new(rec: &mut bcf::Record) -> Result<Vec<Self>, Box<dyn Error>> {
        let is_germline = !rec.info(b"SOMATIC").flag().unwrap_or(false);

        let ann = Annotation::new(rec);

        let prot_change = ann.prot_change.as_str();

        let pos = rec.pos();
        let alleles = rec.alleles();
        let refallele = alleles[0];
        let mut _alleles = Vec::with_capacity(alleles.len() - 1);
        for a in &alleles[1..] {
            if a.len() == 1 && refallele.len() > 1 {
                _alleles.push(Variant::Deletion {
                    pos: pos,
                    len: (refallele.len() - 1) as u32,
                    is_germline: is_germline,
                    prot_change: prot_change.to_owned(),
                });
            } else if a.len() > 1 && refallele.len() == 1 {
                _alleles.push(Variant::Insertion {
                    pos: pos,
                    seq: a[0..].to_owned(),
                    len: (a.len() - 1) as u32,
                    is_germline: is_germline,
                    prot_change: prot_change.to_owned(),
                });
            } else if a.len() == 1 && refallele.len() == 1 {
                _alleles.push(Variant::SNV {
                    pos: pos,
                    alt: a[0],
                    is_germline: is_germline,
                    prot_change: prot_change.to_owned(),
                });
            } else {
                warn!(
                    "Unsupported variant {} -> {}",
                    str::from_utf8(refallele).unwrap(),
                    str::from_utf8(a).unwrap()
                );
            }
        }

        Ok(_alleles)
    }

    pub fn pos(&self) -> u32 {
        match self {
            &Variant::SNV { pos, .. } => pos,
            &Variant::Deletion { pos, .. } => pos,
            &Variant::Insertion { pos, .. } => pos,
        }
    }

    pub fn end_pos(&self) -> u32 {
        match self {
            &Variant::SNV { pos, .. } => pos,
            &Variant::Deletion { pos, len, .. } => pos + len - 1,
            &Variant::Insertion { pos, .. } => pos,
        }
    }

    pub fn is_germline(&self) -> bool {
        match self {
            &Variant::SNV { is_germline, .. } => is_germline,
            &Variant::Deletion { is_germline, .. } => is_germline,
            &Variant::Insertion { is_germline, .. } => is_germline,
        }
    }

    pub fn prot_change(&self) -> String {
        match self {
            &Variant::SNV {
                ref prot_change, ..
            } => prot_change.to_owned(),
            &Variant::Deletion {
                ref prot_change, ..
            } => prot_change.to_owned(),
            &Variant::Insertion {
                ref prot_change, ..
            } => prot_change.to_owned(),
        }
    }

    pub fn frameshift(&self) -> u32 {
        match self {
            &Variant::SNV { .. } => 0,
            &Variant::Deletion { len, .. } => len % 3,
            &Variant::Insertion { ref seq, .. } => (3 - ((seq.len() as u32 - 1) % 3)) % 3,
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
    pub biotype: String,
}

impl Gene {
    pub fn new(id: &str, name: &str, chrom: &str, interval: Interval, biotype: &str) -> Self {
        Gene {
            id: id.to_owned(),
            name: name.to_owned(),
            transcripts: Vec::new(),
            chrom: chrom.to_owned(),
            interval: interval,
            biotype: biotype.to_owned(),
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

impl From<Strand> for PhasingStrand {
    fn from(strand: Strand) -> Self {
        let s = match strand {
            Strand::Forward => PhasingStrand::Forward,
            Strand::Reverse => PhasingStrand::Reverse,
            _ => panic!("Unsupported Strand orientation! Only Forward (+) and Reverse(-) allowed"),
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
    pub fn is_coding(&self) -> bool {
        return !(self.exons.is_empty());
    }
}

#[derive(Debug)]
pub struct Interval {
    pub start: u32,
    pub end: u32,
    pub frame: u32
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
    pub fn new(start: u32, end: u32, frame: &str) -> Self {
        Interval {
            start: start,
            end: end,
            frame: match frame {
                "." => 0 as u32,
                _ => u32::from_str(frame).unwrap(),
            }
        }
    }
/*     pub fn update(&mut self, start: u32, end: u32, frame: &str) -> Result<(), Box<dyn Error>> {
        self.start = start;
        self.end = end;
        self.frame = match frame {
            "." => 0 as u32,
            _ => u32::from_str(frame).unwrap(),
        };
        Ok(())
    } */
}

#[derive(Deserialize, Debug, Serialize, Clone)]
pub struct IDRecord{
    pub id: String,
    pub transcript: String,
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    pub offset: u32,
    pub frame: u32,
    pub freq: f64,
    pub depth: u32,
    pub nvar: u32,
    pub nsomatic: u32,
    pub nvariant_sites: u32,
    pub nsomvariant_sites: u32,
    pub strand: String,
    pub variant_sites: String,
    pub somatic_positions: String,
    pub somatic_aa_change: String,
    pub germline_positions: String,
    pub germline_aa_change: String,
    pub normal_sequence: String,
    pub mutant_sequence: String
}

impl IDRecord {
    pub fn update(&self, rec: &IDRecord, offset: u32, frame: u32, freq: f64, wt_seq: Vec<u8>, mt_seq: Vec<u8>) -> Self {
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

        // let new_freq = match self.freq == rec.freq {
        //     true => self.freq,
        //     false => self.freq * rec.freq,
        // };

        let new_offset = match self.strand == "Forward" {
            true => self.offset + offset,
            false => self.offset - offset,
        };

        debug!("nvars {} {}", self.nvar, rec.nvar);
        IDRecord {
            id: fasta_id,
            transcript: self.transcript.to_owned(),
            gene_id: self.gene_id.to_owned(),
            gene_name: self.gene_name.to_owned(),
            chrom: self.chrom.to_owned(),
            offset: new_offset,
            frame: frame,
            freq: freq,
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
        debug!("Freq to add: {}", freq);
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
            normal_sequence: self.normal_sequence.to_owned(),
            mutant_sequence: self.mutant_sequence.to_owned(),
        }
    }
}