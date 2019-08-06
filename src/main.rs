#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate bio;
extern crate csv;
extern crate itertools;
extern crate rust_htslib;
extern crate vec_map;
#[macro_use]
extern crate serde_derive;
extern crate sha1;

use std::error::Error;
use std::fs::File;
use std::io;
use std::process;

use clap::{App, ArgMatches, SubCommand};

use bio::io::{fasta, gff};
use rust_htslib::{bam, bcf};

pub mod build_reference;
pub mod common;
pub mod filter;
pub mod microphasing_wholegenome;
pub mod peptides;
pub mod microphasing;
pub mod normal_microphasing;

pub fn run() -> Result<(), Box<Error>> {
    let yaml = load_yaml!("cli.yaml");
    let somatic_yaml = load_yaml!("somatic_cli.yaml");
    let germline_yaml = load_yaml!("germline_cli.yaml");
    let filter_yaml = load_yaml!("filter_cli.yaml");
    let build_yaml = load_yaml!("build_ref_cli.yaml");
    let wgs_yaml = load_yaml!("wgs.yaml");
    let matches = App::from_yaml(yaml)
                    .version(env!("CARGO_PKG_VERSION"))
                    .subcommand(SubCommand::from_yaml(somatic_yaml))
                    .subcommand(SubCommand::from_yaml(germline_yaml))
                    .subcommand(SubCommand::from_yaml(filter_yaml))
                    .subcommand(SubCommand::from_yaml(build_yaml))
                    .subcommand(SubCommand::from_yaml(wgs_yaml))
                    .get_matches();

    match matches.subcommand() {
        ("somatic", Some(m)) => run_somatic(m),
        ("normal", Some(m)) => run_normal(m),
        ("build_reference", Some(m)) => run_build(m),
        ("filter", Some(m)) => run_filtering(m),
        ("whole_genome", Some(m)) => run_wg(m),
        _ => Ok(())
    }
}

pub fn run_somatic(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    let mut gtf_reader = gff::Reader::new(io::stdin(), gff::GffType::GTF2);

    let bam_reader = bam::IndexedReader::from_path(matches.value_of("tumor-sample").unwrap())?;

    let bcf_reader = bcf::Reader::from_path(matches.value_of("variants").unwrap())?;

    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    let mut normal_writer = fasta::Writer::to_file(matches.value_of("normal").unwrap())?;

    let mut tsv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(matches.value_of("tsv").unwrap())?;

    let window_len = value_t!(matches, "window-len", u32)?;
    microphasing::phase(
        &mut fasta_reader,
        &mut gtf_reader,
        bcf_reader,
        bam_reader,
        &mut fasta_writer,
        &mut tsv_writer,
        &mut normal_writer,
        window_len,
    )
}

pub fn run_normal(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    let mut gtf_reader = gff::Reader::new(io::stdin(), gff::GffType::GTF2);

    let bam_reader = bam::IndexedReader::from_path(matches.value_of("normal-sample").unwrap())?;

    let bcf_reader = bcf::Reader::from_path(matches.value_of("variants").unwrap())?;

    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    let window_len = value_t!(matches, "window-len", u32)?;
    normal_microphasing::phase(
        &mut fasta_reader,
        &mut gtf_reader,
        bcf_reader,
        bam_reader,
        &mut fasta_writer,
        window_len,
    )
}

pub fn run_build(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    let reference_reader = fasta::Reader::from_file(&matches.value_of("reference").unwrap())?;
    let binary_writer = File::create(&matches.value_of("output").unwrap())?;

    peptides::build(reference_reader, binary_writer)
}

pub fn run_filtering(matches: &ArgMatches) -> Result<(), Box<Error>> {
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    let reference_reader = File::open(&matches.value_of("reference").unwrap())?;

    let mut tsv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&matches.value_of("tsv").unwrap())?;

    let mut tsv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(matches.value_of("tsvoutput").unwrap())?;
    let mut fasta_writer = fasta::Writer::new(io::stdout());
    let mut normal_writer = fasta::Writer::to_file(matches.value_of("normaloutput").unwrap())?;

    peptides::filter(
        reference_reader,
        &mut tsv_reader,
        &mut fasta_writer,
        &mut normal_writer,
        &mut tsv_writer,
    )
}

pub fn run_wg(matches: &ArgMatches) -> Result<(), Box<Error>> {
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



    let bam_reader = bam::IndexedReader::from_path(matches.value_of("tumor-sample").unwrap())?;

    let bcf_reader = bcf::Reader::from_path(matches.value_of("variants").unwrap())?;

    let mut fasta_reader = fasta::IndexedReader::from_file(&matches.value_of("ref").unwrap())?;

    let mut fasta_writer = fasta::Writer::new(io::stdout());

    let only_relevant = matches.is_present("relevant");

    let mut normal_writer = fasta::Writer::to_file(matches.value_of("normal").unwrap())?;

    let mut tsv_writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(matches.value_of("tsv").unwrap())?;

    let window_len = value_t!(matches, "window-len", u32)?;
    microphasing_wholegenome::phase(
        &mut fasta_reader, bcf_reader,
        bam_reader, &mut fasta_writer, &mut tsv_writer, &mut normal_writer, window_len, only_relevant
    )
}

pub fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        process::exit(1);
    }
}
