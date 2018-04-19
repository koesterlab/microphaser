use std::process::Command;
use std::fs;


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
            .arg(format!("target/debug/microphaser {}", cmd))
            .spawn().unwrap().wait().unwrap().success());
}


#[test]
fn test_empty() {
    fs::create_dir("tests/output");
    microphaser("tests/resources/test_empty/empty_test_golden.bam --variants tests/resources/test_empty/empty_test_golden.vcf --ref tests/resources/ensembl_test.fa > tests/output/test_empty.fa < tests/resources/test.gtf")
}


#[test]
fn test_short() {
    fs::create_dir("tests/output");
    microphaser("tests/resources/test_short/shorttest.bam --variants tests/resources/test_short/shorttest.vcf --ref tests/resources/test_short/shorttest.fa > tests/output/test_short.fa < tests/resources/test_short/shorttest.gtf");
    test_output("tests/output/test_short.fa", "tests/resources/test_short/expected_output/test_short.fa");
}

#[test]
fn test_forward() {
    fs::create_dir("tests/output");
    microphaser("tests/resources/ensembl_test.sorted.bam --variants tests/resources/ensembl_test.vcf --ref tests/resources/ensembl_test.fa > tests/output/test_forward.fa < tests/resources/test_forward/test_forward.gtf");
    test_output("tests/output/test_forward.fa", "tests/resources/test_forward/expected_output/test_forward.fa");
}

#[test]
fn test_reverse() {
    fs::create_dir("tests/output");
    microphaser("tests/resources/ensembl_test.sorted.bam --variants tests/resources/ensembl_test.vcf --ref tests/resources/ensembl_test.fa > tests/output/test_reverse.fa < tests/resources/test_reverse/reverse.gtf");
    test_output("tests/output/test_reverse.fa", "tests/resources/test_reverse/expected_output/test_reverse.fa");
}
