# pileups

A fast wolF Flow to build an AD/DP matrix with `bcftools` and indexed lists of `bams`. Works on 1GB machines. Tested up to 250 variants per machine with no issue so far. Performs 100,000 mutations forcecalling in 20 genomes in < 5 minutes if `bam` files are localized. Streaming available, performance limited to ~1 SNV per mem GB and < 500 total tasks per flow.
