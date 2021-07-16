extern crate flate2;
extern crate http;
extern crate hyper;

use std::fs;
use std::path::Path;
use std::process::Command;
// use std::process;
use std::io;

fn test_output(result: &str, expected: &str) {
    assert!(Command::new("diff")
        .arg(result)
        .arg(expected)
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    fs::remove_file(result).unwrap();
}

fn microphaser_somatic(cmd: &str) {
    assert!(Command::new("bash")
        .arg("-c")
        .arg(format!(
            "RUST_BACKTRACE=1 target/debug/microphaser somatic {}",
            cmd
        ))
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
}

fn microphaser_normal(cmd: &str) {
    assert!(Command::new("bash")
        .arg("-c")
        .arg(format!(
            "RUST_BACKTRACE=1 target/debug/microphaser normal {}",
            cmd
        ))
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
}

fn microphaser_filter(cmd: &str) {
    assert!(Command::new("bash")
        .arg("-c")
        .arg(format!(
            "RUST_BACKTRACE=1 target/debug/microphaser filter {}",
            cmd
        ))
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
}

fn microphaser_build(cmd: &str) {
    assert!(Command::new("bash")
            .arg("-c")
            .arg(format!("RUST_BACKTRACE=1 target/debug/microphaser build_reference {}", cmd))
            .spawn().unwrap().wait().unwrap().success());
}

fn download_reference(chrom: &str) -> String {
    let reference = format!("tests/resources/{}.fa", chrom);
    if !Path::new(&reference).exists() {
        let client = hyper::Client::new();
        let res = client
            .get(&format!(
                "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/{}.fa.gz",
                chrom
            ))
            .send()
            .unwrap();
        let mut reference_stream = flate2::read::GzDecoder::new(res).unwrap();
        let mut reference_file = fs::File::create(&reference).unwrap();

        io::copy(&mut reference_stream, &mut reference_file).unwrap();
    }
    assert!(Path::new(&reference).exists());
    if !Path::new(&(reference.clone() + ".fai")).exists() {
        Command::new("samtools")
            .args(&["faidx", &reference])
            .status()
            .expect("failed to create fasta index");
    }

    reference
}

#[test]
fn test_empty() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr14");
    println!("{}", reference);
    microphaser_somatic(&format!(
        "tests/resources/test_forward/forward_test.bam \
         --variants tests/resources/test_empty/empty_test.vcf --tsv tests/output/empty_test.tsv \
         --normal-output tests/output/empty_test.normal.fa \
         --ref {} > tests/output/empty_test.fa < tests/resources/test_forward/forward_test.gtf",
        reference
    ));
    test_output(
        "tests/output/empty_test.fa",
        "tests/resources/test_empty/expected_output/empty_test.fa",
    );
    test_output(
        "tests/output/empty_test.normal.fa",
        "tests/resources/test_empty/expected_output/empty_test.normal.fa",
    );
    test_output(
        "tests/output/empty_test.tsv",
        "tests/resources/test_empty/expected_output/empty_test.tsv",
    );
}

#[test]
fn test_build_ref() {
    fs::create_dir("tests/output");
    microphaser_build("--reference tests/resources/test_build/reference.fa -l4 \
        --output tests/output/test_build_ref.binary > tests/output/test_build_ref.fasta");
    test_output("tests/output/test_build_ref.fasta", "tests/resources/test_build/expected_output/reference_peptides.fasta");
}

#[test]
fn test_filter() {
    fs::create_dir("tests/output");
    microphaser_filter(
        "--reference tests/resources/test_filter/reference.binary -l 9 \
         --tsv tests/resources/test_filter/info.tsv --tsv-output tests/output/info.filtered.tsv \
         --normal-output tests/output/normal.filtered.fa > tests/output/tumor.filtered.fa",
    );
    test_output(
        "tests/output/tumor.filtered.fa",
        "tests/resources/test_filter/expected_output/tumor.filtered.fa",
    );
    test_output(
        "tests/output/normal.filtered.fa",
        "tests/resources/test_filter/expected_output/normal.filtered.fa",
    );
    test_output(
        "tests/output/info.filtered.tsv",
        "tests/resources/test_filter/expected_output/info.filtered.tsv",
    );
}

#[test]
fn test_filter_long() {
    fs::create_dir("tests/output");
    microphaser_filter(
        "--reference tests/resources/test_filter_long/reference.binary -l 9 \
         --tsv tests/resources/test_filter_long/info.tsv --tsv-output tests/output/info.filtered_long.tsv \
         --normal-output tests/output/normal.filtered_long.fa > tests/output/tumor.filtered_long.fa",
    );
    test_output(
        "tests/output/tumor.filtered_long.fa",
        "tests/resources/test_filter_long/expected_output/tumor.filtered_long.fa",
    );
    test_output(
        "tests/output/normal.filtered_long.fa",
        "tests/resources/test_filter_long/expected_output/normal.filtered_long.fa",
    );
    test_output(
        "tests/output/info.filtered_long.tsv",
        "tests/resources/test_filter_long/expected_output/info.filtered_long.tsv",
    );
}

#[test]
fn test_filter_fs() {
    fs::create_dir("tests/output");
    microphaser_filter(
        "--reference tests/resources/test_filter_fs/reference.binary -l 9 \
         --tsv tests/resources/test_filter_fs/info.tsv --tsv-output tests/output/info.filtered_fs.tsv \
         --normal-output tests/output/normal.filtered_fs.fa > tests/output/tumor.filtered_fs.fa",
    );
    test_output(
        "tests/output/tumor.filtered_fs.fa",
        "tests/resources/test_filter_fs/expected_output/tumor.filtered_fs.fa",
    );
    test_output(
        "tests/output/normal.filtered_fs.fa",
        "tests/resources/test_filter_fs/expected_output/normal.filtered_fs.fa",
    );
    test_output(
        "tests/output/info.filtered_fs.tsv",
        "tests/resources/test_filter_fs/expected_output/info.filtered_fs.tsv",
    );
}

#[test]
fn test_forward_somatic() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr14");
    microphaser_somatic(&format!(
        "tests/resources/test_forward/forward_test.bam \
         --variants tests/resources/test_forward/forward_test.vcf \
         --tsv tests/output/forward_test.tsv --ref {} \
         --normal-output tests/output/forward_test.normal.fa \
         > tests/output/forward_test.fa < tests/resources/test_forward/forward_test.gtf",
        reference
    ));
    test_output(
        "tests/output/forward_test.fa",
        "tests/resources/test_forward/expected_output/forward_test.fa",
    );
    test_output(
        "tests/output/forward_test.tsv",
        "tests/resources/test_forward/expected_output/forward_test.tsv",
    );
    test_output(
        "tests/output/forward_test.normal.fa",
        "tests/resources/test_forward/expected_output/forward_test.normal.fa",
    );
}

#[test]
fn test_forward_germline() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr14");
    microphaser_normal(&format!("tests/resources/test_forward/forward_test.bam \
        --variants tests/resources/test_forward/forward_test.germline.vcf \
        --tsv tests/output/forward_test.germline-origins.tsv \
        --ref {} > tests/output/forward_test.germline.fa < tests/resources/test_forward/forward_test.gtf", reference));
    test_output(
        "tests/output/forward_test.germline.fa",
        "tests/resources/test_forward/expected_output/forward_test.germline.fa",
    );
}

#[test]
fn splice_test_forward() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr7");
    microphaser_somatic(&format!("tests/resources/splice_forward_test/INSIG1.test.bam \
        --variants tests/resources/splice_forward_test/INSIG1.test.vcf \
        --tsv tests/output/splice_forward_test.tsv --normal-output tests/output/splice_forward_test.normal.fa --ref {} \
        > tests/output/splice_forward_test.fa < tests/resources/splice_forward_test/INSIG1.test.gtf", reference));
    test_output(
        "tests/output/splice_forward_test.fa",
        "tests/resources/splice_forward_test/expected_output/splice_forward_test.fa",
    );
    test_output(
        "tests/output/splice_forward_test.normal.fa",
        "tests/resources/splice_forward_test/expected_output/splice_forward_test.normal.fa",
    );
    test_output(
        "tests/output/splice_forward_test.tsv",
        "tests/resources/splice_forward_test/expected_output/splice_forward_test.tsv",
    );
}

#[test]
fn splice_test_forward_germline() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr7");
    microphaser_normal(&format!("tests/resources/splice_forward_test/INSIG1.test.bam \
        --variants tests/resources/splice_forward_test/INSIG1.test.germline.vcf \
        --tsv tests/output/splice_forward_test.germline-origins.tsv \
        --ref {} > tests/output/splice_forward_test.germline.fa < tests/resources/splice_forward_test/INSIG1.test.gtf", reference));
    test_output(
        "tests/output/splice_forward_test.germline.fa",
        "tests/resources/splice_forward_test/expected_output/splice_forward_test.germline.fa",
    );
}

#[test]
fn test_reverse() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr1");
    microphaser_somatic(&format!("tests/resources/test_reverse/reverse_test.bam \
        --variants tests/resources/test_reverse/reverse_test.vcf \
        --tsv tests/output/reverse_test.tsv --normal-output tests/output/reverse_test.normal.fa --ref {} \
        > tests/output/reverse_test.fa < tests/resources/test_reverse/reverse_test.gtf", reference));
    test_output(
        "tests/output/reverse_test.fa",
        "tests/resources/test_reverse/expected_output/reverse_test.fa",
    );
    test_output(
        "tests/output/reverse_test.normal.fa",
        "tests/resources/test_reverse/expected_output/reverse_test.normal.fa",
    );
    test_output(
        "tests/output/reverse_test.tsv",
        "tests/resources/test_reverse/expected_output/reverse_test.tsv",
    );
}

//#[test]
//fn test_reverse_germline() {
//    fs::create_dir("tests/output");
//    let reference = download_reference("chr1");
//    microphaser_normal(&format!("tests/resources/test_reverse/reverse_test.bam \
//        --variants tests/resources/test_reverse/reverse_test.germline.vcf \
//        --ref {} > tests/output/reverse_test.germline.fa < tests/resources/test_reverse/reverse_test_germline.gtf", reference));
//    test_output(
//        "tests/output/reverse_test.germline.fa",
//        "tests/resources/test_reverse/expected_output/reverse_test.germline.fa",
//    );
//}

#[test]
fn splice_test_reverse() {
    fs::create_dir("tests/output");
    let reference = download_reference("chr6");
    microphaser_somatic(&format!("tests/resources/splice_reverse_test/MMS22L.test.bam \
       --variants tests/resources/splice_reverse_test/MMS22L.test.vcf --tsv tests/output/splice_reverse_test.tsv \
       --normal-output tests/output/splice_reverse_test.normal.fa --ref {} \
       > tests/output/splice_reverse_test.fa < tests/resources/splice_reverse_test/MMS22L.test.gtf", reference));
    test_output(
        "tests/output/splice_reverse_test.fa",
        "tests/resources/splice_reverse_test/expected_output/splice_reverse_test.fa",
    );
    test_output(
        "tests/output/splice_reverse_test.normal.fa",
        "tests/resources/splice_reverse_test/expected_output/splice_reverse_test.normal.fa",
    );
    test_output(
        "tests/output/splice_reverse_test.tsv",
        "tests/resources/splice_reverse_test/expected_output/splice_reverse_test.tsv",
    );
}

/* #[test]
fn three_way_splice() {
   fs::create_dir("tests/output");
   let reference = download_reference("chr19");
   microphaser_somatic(&format!("tests/resources/three_way_splice/three_way_splice.bam \
       --variants tests/resources/three_way_splice/three_way_splice.vcf --tsv tests/output/three_way_splice.tsv \
       --normal-output tests/output/three_way_splice.wt.fa --ref {} \
       > tests/output/three_way_splice.mt.fa < tests/resources/three_way_splice/three_way_splice.gtf", reference));
   test_output("tests/output/three_way_splice.mt.fa", "tests/resources/three_way_splice/expected_output/three_way_splice.mt.fa");
   test_output("tests/output/three_way_splice.wt.fa", "tests/resources/three_way_splice/expected_output/three_way_splice.wt.fa");
   test_output("tests/output/three_way_splice.tsv", "tests/resources/three_way_splice/expected_output/three_way_splice.tsv");
}

#[test]
fn frameshift_test() {
   fs::create_dir("tests/output");
   let reference = download_reference("chr11");
   microphaser_somatic(&format!("tests/resources/frameshift_test/frameshift_test.bam \
       --variants tests/resources/frameshift_test/frameshift_test.vcf --tsv tests/output/frameshift_test.tsv \
       --normal-output tests/output/frameshift_test.wt.fa --ref {} \
       > tests/output/frameshift_test.mt.fa < tests/resources/frameshift_test/frameshift_test.gtf", reference));
   test_output("tests/output/frameshift_test.mt.fa", "tests/resources/frameshift_test/expected_output/frameshift_test.mt.fa");
   test_output("tests/output/frameshift_test.wt.fa", "tests/resources/frameshift_test/expected_output/frameshift_test.wt.fa");
   test_output("tests/output/frameshift_test.tsv", "tests/resources/frameshift_test/expected_output/frameshift_test.tsv");
} */
