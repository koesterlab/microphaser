pub enum Variant {
    SNV { pos: u32, alt: u8 },
    Insertion { pos: u32, seq: Vec<u8> },
    Deletion { pos: u32, len: u32 }
}


impl Variant {
    pub fn new(rec: &bcf::Record) -> Self {
        // TODO implement conversion
    }
}


pub struct Gene {
    pub name: String,
    pub transcripts: Vec<Transcript>,
    pub interval: Interval
}


impl Gene {
    pub fn start(&self) -> u32 {
        self.interval.start
    }

    pub fn end(&self) -> u32 {
        self.interval.end
    }
}


pub struct Transcript {
    pub exons: Vec<Interval>
}


pub struct Interval {
    pub chrom: String,
    pub start: u32,
    pub end: u32
}
