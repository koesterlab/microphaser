# microphaser

Microphaser is a tool for phasing small tumor DNA sequences - e.g. coding for small peptides - in linear time.  
It can be used in tumor neoantigen prediction to generate the neo-peptidome.

## Installation

  Microphaser will soon be available for installation via Bioconda.
  
## Usage

  To use microphaser, you need a sorted and indexed bam file containing mapped tumor reads, 
  a reference genome in fasta format and matching gtf annotation. Germline and somatic variants should be in the same file in vcf/bcf format where somatic variants have to be flagged with a `SOMATIC` info tag.
  
  You can run microphaser like this:

  ```microphaser tumor.bam -r reference.fa -b variants.vcf -t haplotypes.info.tsv < reference.gtf > haplotypes.fa```

  `haplotypes.fa` will be a fasta file containing the phased haplotype sequences, with unique hash values as identifier.  
  These IDs can be used to find other information about this haplotype in `haplotypes.info.tsv`.  
  Microphaser returns chromosome, position, gene and transcript IDs, gene name,
  strand orientation, number of variants, number of somatic variants and the frequency for every haplotype.
  
