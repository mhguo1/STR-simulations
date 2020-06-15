# STR-simulations

**Simulation of STR reads**

STR_read_simulator.py allows you to generate fake NGS reads around a given STR, which can facilitate downstream optimization of alignment and STR genotyping. This code will make a "pseudo" reference genome with variable STR lengths and then simulate a user specified number of reads.

#Required options:
1) --chr: chromosome of STR. Can be in ucsc format beginning with "chr" or GRCh without the "chr". However, this must match the convention used in the reference genome specified
2) --strstart: The start coordinate of your STR of interest. It is provided in -1 bed coordinates
3) --strend: The end coordinate of your STR of interest.
4) --nucleotide: The repeat nucleotides for your STR. An example is --nucleotide CAG
5) --genotype: The genotypes of the STR you'd like to simulate. If you have a haploid genome (e.g., with X chromosome), then you can just provide a single number (e.g., 15). If you have a diploid genome, you can provide comma separated values (e.g., 15,30). If you have a diploid genome but would like to simulate a homozygous genotype, you can provide a single number or a comma separated number (i.e., 15 and 15,15 are treated equivalently)
6) --ref: This is the hard path to your reference genome. It should be a gzipped fasta file. All chromosomes can be present in the ref file. The current default is /project/jcreminslab/guomic_projects/ref/genomes/hg38.fa.gz
7) --out: This is the prefix for your output. output will be written to the current directory. If paired end is specified, then output files will be ${out}_1.fastq.gz and ${out}_2.fastq.gz. If single end is specified, then output file will be ${out}.fastq.gz. Default is "temp"

#Additional options
1) --stutter: This allows you to simulate stutter in your STR genotypes. The value provided is the standard deviation of your STR genotypes in repeat units. Default is 0.1.
2) --window: A window upstream and downstream of the STR to simulate. For example, if 10000 is provided, the code will simulate reads for 10kb upstream and 10kb downstream of the STR. A value of at least 1000 should be provided. Default is 10000.
3) --coverage: Coverage of the region to simulate. A value of 10 would simulate approximately 10x coverage. Default is 10
4) --numreads: Number of reads to simulate. This is only used if --coverage is not provided
5) --insertmean: Mean distance between read pairs for PE reads. Only used if --pairedend mode is specified. Default is 500.
6) --insertstdev: Mean distance between read pairs for PE reads. Only used if --pairedend mode is specified. Default is 100.
7) --pairedend: Specifies usage of paired end mode
8) --readlen: Specifies length of reads. Default is 150 bp.

**An example**
Here, we will simulate reads around HTT CAG expansion. STR repeat lengths of 19 and 50 will be simulated. Paired end mode with 150 bp PE reads are used. Simulations will be done to 10x coverage. 

python STR_read_simulator.py --chr chr4 --strstart 3074877 --strend 3074933 --nucleotide CAG --genotype 19,50 --stutter 0.2 --out test --window 10000 --coverage 10 --insertmean 500 --insertstdev 120 --pairedend --readlen 150 
