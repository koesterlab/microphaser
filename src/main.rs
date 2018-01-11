#[macro_use]
extern crate log;
#[macro_use]
extern crate clap;
extern crate vec_map;
extern crate bio;
extern crate itertools;
extern crate rust_htslib;

use std::process;
use std::error::Error;
use std::io;

use clap::App;

use rust_htslib::{bam, bcf};
use bio::io::{fasta, gff};

pub mod microphasing;
pub mod common;


pub fn run() -> Result<(), Box<Error>> {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .get_matches();

    let bam_reader = bam::IndexedReader::from_path(matches.value_of("tumor-sample").unwrap())?;
    let bcf_reader = bcf::Reader::from_stdin()?;
    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;
    let mut gtf_reader = gff::Reader::from_file(
        matches.value_of("annotation").unwrap(), gff::GffType::GTF2
    )?;
    let mut fasta_writer = fasta::Writer::new(io::stdout());
    let window_len = value_t!(matches, "window-len", u32)?;
    microphasing::phase(
        &mut fasta_reader, &mut gtf_reader, bcf_reader,
        bam_reader, &mut fasta_writer, window_len
    )
}


pub fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        process::exit(1);
    }
}
