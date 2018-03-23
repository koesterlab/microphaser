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
    microphaser("tests/resources/test_empty/empty_test_golden.bam --track tests/resources/test.gtf --ref tests/resources/ensembl_test.fa > tests/output/test_empty.fa < tests/resources/test_empty/empty_test_golden.vcf")
}


#[test]
fn test_short() {
    fs::create_dir("tests/output");
    microphaser("tests/resources/test_short/shorttest.bam --track tests/resources/test_short/shorttest.gtf --ref tests/resources/test_short/shorttest.fa > tests/output/test_short.fa < tests/resources/test_short/shorttest.vcf");
    test_output("tests/output/test_short.fa", "tests/resources/test_short/expected_output/test_short.fa");
}
