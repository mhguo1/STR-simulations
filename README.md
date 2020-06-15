# STR-simulations

**Simulation of STR reads**
STR_read_simulator.py allows you to generate fake NGS reads around a given STR, which can facilitate downstream optimization of alignment and STR genotyping. This code will make a "pseudo" reference genome with variable STR lengths and then simulate a user specified number of reads.

#Required options:
--chr: chromosome of STR. Can be in ucsc format beginning with "chr" or GRCh without the "chr". However, this must match the convention used in the reference genome specified
--strstart: The start coordinate of your STR of interest. It is provided in -1 bed coordinates
--strend: The end coordinate of your STR of interest.
--nucleotide: The repeat nucleotides for your STR. An example is --nucleotide CAG
--genotype: The genotypes of the STR you'd like to simulate. If you have a haploid genome (e.g., with X chromosome), then you can just provide a single number (e.g., 15). If you have a diploid genome, you can provide comma separated values (e.g., 15,30). If you have a diploid genome but would like to simulate a homozygous genotype, you can provide a single number or a comma separated number (i.e., 15 and 15,15 are treated equivalently)
--ref: This is the hard path to your reference genome. It should be a gzipped fasta file. All chromosomes can be present in the ref file. The current default is /project/jcreminslab/guomic_projects/ref/genomes/hg38.fa.gz
--out: This is the prefix for your output. output will be written to the current directory. If paired end is specified, then output files will be ${out}_

#Additional options
--stutter: This allows you to simulate stutter in your STR genotypes. The value provided is the standard deviation of your STR genotypes in repeat units. Default is 0.1.



parser.add_option("--stutter", action="store", dest="stutter", default=0.1) #stutter, as expressed in SD of # of repeats
parser.add_option("--ref", action="store", dest="ref", default="/project/jcreminslab/guomic_projects/ref/genomes/hg38.fa.gz")
parser.add_option("--out", action="store",dest="outprefix", default="temp") #output file
parser.add_option("--window", action="store", dest="window", default=10000) #Window around STR to simulate
parser.add_option("--coverage", action="store", dest="coverage") #X coverage
parser.add_option("--numreads", action="store", dest="numreads", default=10000) # of reads to simulate. Used only if cov is not provided
parser.add_option("--insertmean", action="store", dest="insertmean", default=500) # insert mean. Used only if PE
parser.add_option("--insertstdev", action="store", dest="insertstdev", default=100) #st dev of insert length. Used only if PE
parser.add_option("--pairedend", action="store_true", dest="pairedend") #Paired end
parser.add_option("--readlen", action="store", dest="readlen", default=150) #length of reads
