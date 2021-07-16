# Microphaser

Microphaser is a tool for phasing small tumor DNA sequences - e.g. coding for small peptides - in linear time.  
It can be used in tumor neoantigen prediction to generate the neo-peptidome.

## Installation

  Microphaser is available for installation via conda. 
  Use `conda install -c biconda microphaser` to easily install the current version.
  
## Usage

### Input
  To use microphaser, you need the following input files:
  
  * a sorted and indexed bam file containing mapped tumor reads
  * a reference genome in fasta format
  * a matching gene and transcript annotation in gtf format
  * a bcf/vcf file containing germline and somatic variants, where somatic variants should be flagged with a ```SOMATIC``` INFO tag
  * optional: a bcf/vcf file containing only germline variants
  
### Output
  Microphaser returns three important files:

  * two filtered fasta files containing all neopeptides and their wildtype counterparts for further use with MHC-binding prediction tools
  * an info file in tsv format containing meta-information about every neopeptide
  
  
  The info table consist of the following fields:
  * id: peptide identifier as found in the fasta files
  * transcript: Ensembl transcript name
  * gene_id: Ensembl gene name
  * gene_name: Gene symbol
  * chrom: Chromosome
  * offset: Position of the neopeptide on the chromosome
  * freq: Frequency of the neopeptide occurring in all reads overlapping the peptide position
  * depth: Read depth at the peptide position
  * nvar: number of variants in the neopeptide
  * nsomatic: number of somatic variants in the neopeptide
  * nvariant_sites: number of variant sites in the range of the neopeptide
  * nsomvariant_sites: number of somatic variant sites in the range of the neopeptide
  * strand: Strand orientation of the transcript (forward or reverse)
  * somatic_positions: Positions of the somatic variants in the neopeptide
  * somatic_aa_change: Somatic Amino Acid changes occuring in the neopeptide
  * germline_positions: Positions of germline variants in the neopeptide
  * germline_aa_change: Germline Amino Acid changes occuring in the neopeptide
  * normal_sequence: Nucleotide sequence of the wildtype peptide
  * mutant_sequence: Nucleotide sequence of the neopeptide
  
### Run
  
  Currently, microphaser consists of four different submodules:

  - somatic (returns neopeptides and their corresponding normal peptides)
  - normal (returns all normal peptides of the patient)
  - build_reference (returns a binary file representing the patients normal peptidome)
  - filter (compares neopeptides against the normal peptidome and removes self-similar candidates)

  
  You can run microphaser like this:
  
  Phasing of the tumor reads and variants:

  ```
  microphaser somatic tumor.bam -r reference.fasta -b all_variants.bcf -t neopeptides.info.tsv -n peptides.wt.fasta < annotation.gtf > peptides.mt.fasta
  ```
  
  Generation of the patients healthy peptidome:
  
  ```
  microphaser normal normal.bam -r reference.fasta -b germline_variants.bcf < annotation.gtf > healthy_peptides.fasta
  ```
  
  Building the reference binary file of the healthy peptidome: 
  
  ```
  microphaser build_reference -r healthy_peptides.fasta -o peptides.bin > peptides.translated.fasta
  ```

  Filtering the neopeptide candidates from subcommand ```microphaser somatic```: 

  ```
  microphaser filter -r peptides.bin -t neopeptides.info.tsv -o neopeptides.filtered.info.tsv -n normal_peptides.filtered.fasta > neopeptides.filtered.fasta
  ```
  
## Authors

* Jan Forster (https://github.com/jafors)
* Johannes Köster (https://koesterlab.github.io) 
