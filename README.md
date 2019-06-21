# microphaser

Microphaser is a tool for phasing small tumor DNA sequences - e.g. coding for small peptides - in linear time.  
It can be used in tumor neoantigen prediction to generate the neo-peptidome.

## Installation

  Microphaser will soon be available for installation via Bioconda.
  
## Usage

### Input
  To use microphaser, you need the following input files:
  
  - a sorted and indexed bam file containing mapped tumor reads
  - your reference genome in fasta format
  - the matching gene and transcript annotation in gtf format
  - a bcf/vcf file containing germline and somatic variants, where the somatic variants should be flagged with a ```SOMATIC``` INFO tag
  - optional: a bcf/vcf file containing only germline variants
  
### Output
  Microphaser returns three important files
  - two filtered fasta files containing all neo-peptides and their wildtype counterparts for further use with MHC-binding prediction tools
  - an info file in tsv format containing meta-information about every neo-peptide
  
  
  The info table consist of the following fields:
  - id: peptide identifier as found in the fasta files
  - transcript: Ensembl transcript name
  - gene_id: Ensembl gene name
  - gene_name: Gene symbol
  - chrom: Chromosome
  - offset: Position of the neopeptide on the chromosome
  - freq: Frequency of the neopeptide occurring in all reads overlapping the peptide position
  - nvar: number of variants in the neopeptide
  - nsomatic: number of somatic variants in the neopeptide
  - nvariant_sites: number of variant sites in the range of the neopeptide
  - nsomvariant_sites: number of somatic variant sites in the range of the neopeptide
  - strand: Strand orientation of the transcript (forward or reverse)
  - somatic_positions: Positions of the somatic variants in the neopeptide
  - somatic_aa_change: Somatic Amino Acid changes occuring in the neopeptide
  - germline_positions: Positions of germline variants in the neopeptide
  - germline_aa_change: Germline Amino Acid changes occuring in the neopeptide
  - normal_sequence: Nucleotide sequence of the wildtype peptide
  - mutant_sequence: Nucleotide sequence of the neopeptide
  
### Run
  
  Currently, microphaser consists of four different submodules:
  - somatic (returns neo-peptides and their corresponding wildtype peptides)
  - normal (returns all wildtype peptides of the patient)
  - build_reference (returns a binary file representing the patients wildtype peptidome)
  - filter (compares neo-peptides against the wildtype peptidome and removes self-similar candidates)
  
  You can run microphaser like this:
  
  Phasing of the tumor reads and variants:

  ```microphaser somatic tumor.bam -r reference.fa -b all_variants.bcf -t neopeptides.info.tsv -n wildtype_peptides.fa < reference.gtf > neopeptides.fa```
  
  Generation of the patients germline peptidome:
  
  ```microphaser normal healthy.bam -r reference.fa -b germline_variants.bcf -n wildtype_peptides.fa < reference.gtf > germline_peptidome.fa```
  
  Building the reference binary file of the germline peptidome: 
  
  ```microphaser build_reference germline_proteome.fa > germline_proteome.bin```

  Filtering the neopeptide candidates from subcommand ```microphaser somatic```: 

  ```microphaser filter -r germline_proteome.bin -t neopeptides.info.tsv -o neopeptides.filtered.info.tsv -n wildtype_peptides.filtered.fa > neopeptides.filtered.fa```
  
  
