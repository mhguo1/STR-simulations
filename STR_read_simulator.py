#!/usr/bin/python

import optparse
import sys
import gzip
import numpy as np
import random

#Parse options
parser = optparse.OptionParser()

#Required options
parser.add_option("--chr", action="store",dest="chr") #chromosome
parser.add_option("--strstart", action="store",dest="strstart") #start coordinate of STR in bed format
parser.add_option("--strend", action="store",dest="strend") #end coordinate of STR in bed format
parser.add_option("--nucleotide", action="store", dest="nucleotide") #repeat letters
parser.add_option("--genotype", action="store", dest="genotype") #genotypes
parser.add_option("--ref", action="store", dest="ref", default="/project/jcreminslab/guomic_projects/ref/genomes/hg38.fa.gz")
parser.add_option("--out", action="store",dest="outprefix", default="temp") #output file

#Optional
parser.add_option("--stutter", action="store", dest="stutter", default=0.1) #stutter, as expressed in SD of # of repeats
parser.add_option("--window", action="store", dest="window", default=10000) #Window around STR to simulate
parser.add_option("--coverage", action="store", dest="coverage") #X coverage
parser.add_option("--numreads", action="store", dest="numreads", default=10000) # of reads to simulate. Used only if cov is not provided
parser.add_option("--insertmean", action="store", dest="insertmean", default=500) # insert mean. Used only if PE
parser.add_option("--insertstdev", action="store", dest="insertstdev", default=100) #st dev of insert length. Used only if PE
parser.add_option("--pairedend", action="store_true", dest="pairedend") #Paired end
parser.add_option("--readlen", action="store", dest="readlen", default=150) #length of reads


options, args = parser.parse_args()

window=int(options.window)
chr=options.chr
str_start=int(options.strstart)
str_end=int(options.strend)
stutter=float(options.stutter)
readlen=int(options.readlen)

#Boundaries of region to simulate
region_start=str_start-window
region_end=str_end+window

pad=readlen*5

#Generate error messages
if options.chr is None:
   sys.stdout.write("Error: No chromosome for STR provided. Killing\n")
   sys.exit()
   
if options.strstart is None:
   sys.stdout.write("Error: No start position for STR provided. Killing\n")
   sys.exit()
   
if options.strend is None:
   sys.stdout.write("Error: No end position for STR provided. Killing\n")
   sys.exit()
   
if options.nucleotide is None:
   sys.stdout.write("Error: No nucleotides provided for STR\n")
   sys.stdout.write("An example is --nucleotide CAG\n")
   sys.stdout.write("Killing\n")
   sys.exit()

if options.nucleotide is None:
   sys.stdout.write("Error: No genotypes provided for STR\n")
   sys.stdout.write("An example if you want a haploid with 15 repeats: --genotype 15\n")
   sys.stdout.write("An example if you want a diploid with 15,30 repeats: --genotype 15,30\n")
   sys.stdout.write("Killing\n")
   sys.exit()
   
#Generate notes
if options.pairedend is None:
   sys.stdout.write("Note: insert mean and standard deviation are ignored as single end mode\n")

if options.coverage is not None:
   sys.stdout.write("Note: using provided coverage to calculate number of reads to simulate\n")

   
   
   
###Calculate number of reads to simulate
if options.coverage is not None:
   numreads=round(int(options.coverage)*(region_end-region_start)/(readlen))
   if options.pairedend:
      numreads=round(numreads/2)
else:
   numreads=int(options.numreads)


###generate simulated expansion lengths
genotypes=str(options.genotype).split(",")
if len(set(genotypes))==1:
   simrepeatlengths = list(np.random.normal(float(genotypes[0]), stutter, int(numreads)))
if len(set(genotypes))>1:
   simrepeatlengths = list(np.random.normal(float(genotypes[0]), stutter, int(numreads/2)))+list(np.random.normal(float(genotypes[1]), stutter, int(numreads/2)))

#In case repeat lenghts are < 0, set to 0
simrepeatlengths=map(lambda x: 0 if x<0 else x, simrepeatlengths)

##Read in reference genome and region of interest
refgenomefile=gzip.open(options.ref, "rb")
refgenomechr=""
if chr.isdigit():
   chrtester=">"+str(chr)
else:
   chrtester=">"+chr
   
for line_g1 in refgenomefile:
  if line_g1.rstrip().split(" ")[0]==chrtester: #find fragment corresponding to chromosome of interest
    while True: #Read in all of ref genome corresponding to that chromosome and stop once finished reading that chromosome
      temp=refgenomefile.readline().rstrip()
      refgenomechr=refgenomechr+temp
      if temp[0]==">":
        break
    break
refgenomefile.close()
genomefrag1=refgenomechr[(region_start-pad):(str_start-1)] #Region upstream of STR
genomefrag2=refgenomechr[str_end:(region_end+pad)] #Region downstream of STR
   
  
   
qualstring="I"*readlen #generate fake quality string

if options.pairedend: ##Generate reads if PE
   filename1=options.outprefix+"_1.fastq.gz"
   filename2=options.outprefix+"_2.fastq.gz"
   outfile1=gzip.open(filename1, "wb")
   outfile2=gzip.open(filename2, "wb") 
   for i in range(0,len(simrepeatlengths)):
      simgenome=genomefrag1+options.nucleotide*int(simrepeatlengths[i])+genomefrag2

      #Simulate reads
      read1_start=random.choice(range(pad, len(simgenome)-pad))
      read1=simgenome[read1_start:(read1_start+readlen)]
   
      insertlen=int(np.random.normal(float(options.insertmean), float(options.insertstdev), 1))  
      read2=simgenome[(read1_start+readlen+insertlen):(read1_start+2*readlen+insertlen)]
      
      #Write out reads
      outfile1.write("@read"+str(i+1)+"/1\n")
      outfile1.write(read1+"\n")
      outfile1.write("+\n")
      outfile1.write(qualstring+"\n")
      
      outfile2.write("@read"+str(i+1)+"/2\n")
      outfile2.write(read2+"\n")
      outfile2.write("+\n")
      outfile2.write(qualstring+"\n")
      
   outfile1.close()
   outfile2.close()
else: #Generate reads if SE
   filename=options.outprefix+".fastq.gz"
   outfile=gzip.open(filename, "wb")
   for i in range(0,len(simrepeatlengths)):
      simgenome=genomefrag1+options.nucleotide*int(simrepeatlengths[i])+genomefrag2
      
      #Simulate SE reads
      read1_start=random.choice(range(pad, len(simgenome)-pad))
      read1=simgenome[read1_start:(read1_start+readlen)]
      
      #Write out reads
      outfile.write("@read"+str(i+1)+"\n")
      outfile.write(read1+"\n")
      outfile.write("+\n")
      outfile.write(qualstring+"\n")
   outfile.close()

   
##completion message##
sys.stdout.write("STR position: "+str(chr)+":"+str(str_start)+"-"+str(str_end)+"\n")
sys.stdout.write("STR genotypes: " + str(genotypes)+"\n")
sys.stdout.write("STR nucleotide: "+str(options.nucleotide) + "\n")
sys.stdout.write("STR stutter model: "+str(stutter) + "\n")
sys.stdout.write("Simulated +/-"+str(window)+"kb around STR" + "\n")
sys.stdout.write("STR stutter model: "+str(stutter) + "\n")
sys.stdout.write("Number of reads were generated"+str(numreads) +"\n")
sys.stdout.write("Read lengths were "+str(readlen)+"\n")
if options.pairedend is not None:
   sys.stdout.write("Paired end mode was used\n")
   sys.stdout.write("Mean insert length was " + str(options.insertmean)+"\n")
   sys.stdout.write("Standard deviation of insert length was " + str(options.insertstdev)+"\n")
else:
   sys.stdout.write("Single end mode was used\n")
sys.stdout.write("Congratulations, you have generated some fake STR repeats!\n")
sys.stdout.write("If you have found this program of use, you can buy Mike a beer at some point\n")
      
