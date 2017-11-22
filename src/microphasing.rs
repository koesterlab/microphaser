use std::error::Error;

use bio::io::fasta;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::prelude::*;

use common::{Gene, Variant};


pub fn supports_variant(read: &bam::Record, variant: &Variant) -> Result<bool, Box<Error>> {
    match variant {
        Variant::SNV { pos: pos, alt: alt } => {
            let b = read.seq()[read.cigar().read_pos(pos)?.unwrap()];
            Ok(b == alt)
        },
        Variant::Insertion { pos: pos, seq: seq } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar() {
                match c {
                    Cigar::Ins(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        },
        Variant::Deletion { pos: pos, len: varlen } => {
            // TODO compare the two using a pair HMM or use cigar string
            for c in read.cigar() {
                match c {
                    Cigar::Del(_) => return Ok(true),
                    _ => ()
                }
            }
            Ok(false)
        }
    }
}


pub struct ObservationMatrix {
    inner: Vec<u64>,
    reads: BTreeMap<u32, Vec<&bam::Record>>,
    ncols: u32
}


impl ObservationMatrix {
    pub fn shrink_left(&mut self, k: u32) {
        self.ncols -= k;
        for b in inner.iter_mut() {
            *b = *b & (2.pow(self.ncols) - 1);
        }
    }

    pub fn extend_right(&mut self, k: u32) {
        self.ncols += k;
        for b in inner.iter_mut() {
            *b = *b << k;
        }
    }

    /// Remove all reads that do not enclose interval end.
    pub fn cleanup_reads(&mut self, interval_end: u32) {
        self.reads = self.reads.split_off(interval_end);
    }

    /// Add read, while considering given interval end.
    pub fn push_read(&mut self, read: &bam::Record, interval_end:u32, variants: &[Variant]) {
        let end_pos = read.cigar().end_pos();
        if end_pos >= interval_end {
            // only insert if end_pos is larger than the interval end
            self.reads.entry(end_pos).or_insert_with(|| Vec::new()).push(read);
            let mut b = 0;
            for (i, variant) in variants.iter().enumerate() {
                if supports_variant(read, variant) {
                    b |= 1 << i;
                }
            }
            self.inner.push(b);
        }
    }
}


pub fn microphasing(
    gene: Gene,
    fasta_reader: &mut fasta::IndexedReader,
    read_buffer: &mut ReadBuffer,
    variant_buffer: &mut VariantBuffer,
    window_len: u32
) {
    let observations = ObservationMatrix::new();

    let mut offset = gene.start();
    variant_buffer.fetch(offset, offset + window_len);
    for read in read_buffer.fetch(offset, offset + window_len) {
        observations.push_read(read, offset + window_len);
    }

    while offset + window_len < gene.end() {


        // advance window to next position
        let deleted = variant_buffer.shrink_left();
        observations.shrink_left(deleted);
        let added = variant_buffer.extend_right();
        observations.extend_right(added);
        offset += 1;
        observations.cleanup_reads(offset + window_len);
        // add new reads
        for read in read_buffer.starts_at(offset) {
            observations.push_read(read, offset + window_len);
        }
        observations.push_read()
    }
}
