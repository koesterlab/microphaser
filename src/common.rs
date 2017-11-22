pub enum Variant {
    SNV { pos: u32, alt: u8 },
    Insertion { pos: u32, seq: Vec<u8> },
    Deletion { pos: u32, len: u32 }
}


pub struct Gene {
    name: String,
    transcripts: Vec<Transcript>,
    interval: Interval
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
    exons: Vec<Interval>
}


pub struct Interval {
    chrom: String,
    start: u32,
    end: u32
}
