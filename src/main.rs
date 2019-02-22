#[macro_use]
extern crate log;
extern crate fern;
extern crate env_logger;
#[macro_use]
extern crate clap;
extern crate vec_map;
extern crate bio;
extern crate itertools;
extern crate rust_htslib;
extern crate csv;
#[macro_use]
extern crate serde_derive;
extern crate sha1;

use std::process;
use std::error::Error;
use std::io;



use clap::{ArgMatches, App, SubCommand};

use rust_htslib::{bam, bcf};
use bio::io::{fasta, gff};

pub mod microphasing;
pub mod common;
pub mod normal_microphasing;
pub mod filter;




pub fn run() -> Result<(), Box<Error>> {
    println!("Starting");
    let yaml = load_yaml!("cli.yaml");
    let somatic_yaml = load_yaml!("somatic_cli.yaml");
    let germline_yaml = load_yaml!("germline_cli.yaml");
    let filter_yaml = load_yaml!("filter_cli.yaml");
    let matches = App::from_yaml(yaml)
                    .version(env!("CARGO_PKG_VERSION"))
                    .subcommand(SubCommand::from_yaml(somatic_yaml))
                    .subcommand(SubCommand::from_yaml(germline_yaml))
                    .subcommand(SubCommand::from_yaml(filter_yaml))
                    .get_matches();

    match matches.subcommand() {
        ("somatic", Some(m)) => run_somatic(m),
        ("normal", Some(m)) => run_normal(m),
        ("filter", Some(m)) => run_filtering(m),
        _ => Ok(())
    }

}

pub fn run_somatic(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
                   .format(|out, message, _| out.finish(format_args!("{}", message)))
                   .level(
                       if matches.is_present("verbose") {
                           log::LevelFilter::Debug
                       } else {
                           log::LevelFilter::Info
                       }
                   )
                   .chain(std::io::stderr())
                   .apply().unwrap();


    let mut gtf_reader = gff::Reader::new(
        io::stdin(), gff::GffType::GTF2
    );

    let bam_reader = bam::IndexedReader::from_path(matches.value_of("tumor-sample").unwrap())?;

    let bcf_reader = bcf::Reader::from_path(matches.value_of("variants").unwrap())?;

    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    let only_relevant = matches.is_present("relevant");

    let mut prot_writer = fasta::Writer::to_file(matches.value_of("proteome").unwrap())?;

    let mut tsv_writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(matches.value_of("tsv").unwrap())?;

    let window_len = value_t!(matches, "window-len", u32)?;
    microphasing::phase(
        &mut fasta_reader, &mut gtf_reader, bcf_reader,
        bam_reader, &mut fasta_writer, &mut tsv_writer, &mut prot_writer, window_len, only_relevant
    )
}

pub fn run_normal(matches: &ArgMatches) -> Result<(), Box<Error>> {
    println!("Normal Mode");
    fern::Dispatch::new()
                   .format(|out, message, _| out.finish(format_args!("{}", message)))
                   .level(
                       if matches.is_present("verbose") {
                           log::LevelFilter::Debug
                       } else {
                           log::LevelFilter::Info
                       }
                   )
                   .chain(std::io::stderr())
                   .apply().unwrap();


    let mut gtf_reader = gff::Reader::new(
        io::stdin(), gff::GffType::GTF2
    );

    let bam_reader = bam::IndexedReader::from_path(matches.value_of("normal-sample").unwrap())?;

    let bcf_reader = bcf::Reader::from_path(matches.value_of("variants").unwrap())?;

    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    let window_len = value_t!(matches, "window-len", u32)?;
    normal_microphasing::phase(
        &mut fasta_reader, &mut gtf_reader, bcf_reader,
        bam_reader, &mut fasta_writer, window_len
    )
}

pub fn run_filtering(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
                   .format(|out, message, _| out.finish(format_args!("{}", message)))
                   .level(
                       if matches.is_present("verbose") {
                           log::LevelFilter::Debug
                       } else {
                           log::LevelFilter::Info
                       }
                   )
                   .chain(std::io::stderr())
                   .apply().unwrap();

    let normal_reader = fasta::Reader::from_file(&matches.value_of("normal").unwrap())?;

    let tumor_reader = fasta::Reader::from_file(&matches.value_of("tumor").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    filter::filter(normal_reader, tumor_reader, &mut fasta_writer)
}

pub fn main() {

    if let Err(e) = run() {
        error!("{}", e);
        process::exit(1);
    }
}
