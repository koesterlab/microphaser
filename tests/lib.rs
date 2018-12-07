extern crate hyper;
extern crate flate2;
extern crate http;

use std::process::Command;
use std::fs;
use std::path::Path;
// use std::process;
use std::error::Error;
use std::io;
use http::Uri;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


fn microphaser(cmd: &str) {
    assert!(Command::new("bash")
            .arg("-c")
            .arg(format!("RUST_BACKTRACE=1 target/debug/microphaser {}", cmd))
            .spawn().unwrap().wait().unwrap().success());
}


fn download_reference(chrom: &str) -> String {
    let reference = format!("tests/resources/{}.fa", chrom);
    if !Path::new(&reference).exists() {
        let client = hyper::Client::new();
        let res = client.get(
            &format!("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/{}.fa.gz", chrom)
        ).send().unwrap();
        let mut reference_stream = flate2::read::GzDecoder::new(res).unwrap();
        let mut reference_file = fs::File::create(&reference).unwrap();

        io::copy(&mut reference_stream, &mut reference_file).unwrap();
    }
    assert!(Path::new(&reference).exists());
    if !Path::new(&(reference.clone() + ".fai")).exists() {
        Command::new("samtools").args(&["faidx", &reference])
                                .status()
                                .expect("failed to create fasta index");
    }

    reference
}

#[test]
fn test_empty() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr14");
    println!("{}",reference);
    microphaser(&format!("tests/resources/test_forward/forward_test.bam --variants tests/resources/test_forward/empty_test.vcf --tsv tests/output/empty_test.tsv --proteome tests/output/empty_test.prot.fa --ref {} > tests/output/empty_test.fa < tests/resources/test_forward/forward_test.gtf", reference));
    test_output("tests/output/empty_test.fa", "tests/resources/test_forward/expected_output/empty_test.fasta.out");
    test_output("tests/output/empty_test.prot.fa", "tests/resources/test_forward/expected_output/empty_test.prot.out");
}

#[test]
fn test_forward() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr14");
    microphaser(&format!("tests/resources/test_forward/forward_test.bam --variants tests/resources/test_forward/forward_test.vcf --tsv tests/output/forward_test.csv --proteome tests/output/forward_test.prot.fa --ref {} > tests/output/forward_test.fa < tests/resources/test_forward/forward_test.gtf", reference));
    test_output("tests/output/forward_test.fa", "tests/resources/test_forward/expected_output/forward_test.fasta.out");
    test_output("tests/output/forward_test.prot.fa", "tests/resources/test_forward/expected_output/forward_test.prot.out");
    test_output("tests/output/forward_test.csv", "tests/resources/test_forward/expected_output/forward_test.csv");
}

#[test]
fn splice_test_forward() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr7");
    microphaser(&format!("tests/resources/splice_forward_test/INSIG1.test.bam --variants tests/resources/splice_forward_test/INSIG1.test.vcf --tsv tests/output/splice_forward_test.tsv --proteome tests/output/splice_forward_test.prot.fa --ref {} > tests/output/splice_forward_test.fa < tests/resources/splice_forward_test/INSIG1.test.gtf", reference));
    test_output("tests/output/splice_forward_test.fa", "tests/resources/splice_forward_test/expected_output/INSIG1.test.fa");
    test_output("tests/output/splice_forward_test.prot.fa", "tests/resources/splice_forward_test/expected_output/INSIG1.test.prot.fa");
    test_output("tests/output/splice_forward_test.tsv", "tests/resources/splice_forward_test/expected_output/INSIG1.test.tsv");
}

#[test]
fn test_reverse() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr1");
    microphaser(&format!("tests/resources/test_reverse/reverse_test.bam --variants tests/resources/test_reverse/reverse_test.vcf --tsv tests/output/reverse_test.tsv --proteome tests/output/reverse_test.prot.fa --ref {} > tests/output/reverse_test.fa < tests/resources/test_reverse/reverse_test.gtf", reference));
    test_output("tests/output/reverse_test.fa", "tests/resources/test_reverse/expected_output/reverse_test.fasta.out");
    test_output("tests/output/reverse_test.prot.fa", "tests/resources/test_reverse/expected_output/reverse_test.prot.out");
    test_output("tests/output/reverse_test.tsv", "tests/resources/test_reverse/expected_output/reverse_test.tsv");
}

#[test]
fn splice_test_reverse() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr6");
    microphaser(&format!("tests/resources/splice_reverse_test/MMS22L.test.bam --variants tests/resources/splice_reverse_test/MMS22L.test.vcf --tsv tests/output/splice_reverse_test.tsv --proteome tests/output/splice_reverse_test.prot.fa --ref {} > tests/output/splice_reverse_test.fa < tests/resources/splice_reverse_test/MMS22L.test.gtf", reference));
    test_output("tests/output/splice_reverse_test.fa", "tests/resources/splice_reverse_test/expected_output/MMS22L.test.fa");
    test_output("tests/output/splice_reverse_test.prot.fa", "tests/resources/splice_reverse_test/expected_output/MMS22L.test.prot.fa");
    test_output("tests/output/splice_reverse_test.tsv", "tests/resources/splice_reverse_test/expected_output/MMS22L.test.tsv");
}
