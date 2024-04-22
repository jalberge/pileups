# pileups

A fast wolF Flow to build an AD/DP matrix with `bcftools` and indexed lists of `bams`. Works on 1GB machines. Tested up to 250 variants per machine with no issue so far. Performs 100,000 mutations forcecalling in 20 genomes in < 10 minutes if `bam` files are localized. Streaming available, performance limited to ~1 SNV per mem GB per worker and < 500 total tasks per flow.

## Dockerfile

Unlike most wolF repos, this one only hosts a generic htslib+samtools+bcftools install with gcp enabled. Is it hosted on `gcr.io/broad-getzlab-workflows/samtools`

## to-do

* Join with original MAF file to rescue funcotator
* Pipe to Sequoia directly
* Manually verify counts
* Option to report more than AD/DP?
