#! /usr/bin/python3
     
from Bio import SeqIO
import glob
import sys

# arguments
arguments = len(sys.argv) - 1
if(arguments!=2): 
    quit()


inputFiles = sys.argv[1]+'/*.fasta'
outputFile = sys.argv[2]

#######################################
# 1. reads sequences and concatenates them
#######################################
#inputFiles = '/home/pturjanski/kapow/data/omics/dferreiro/SRR1570769/*.fasta'
# inputFiles = 'input/*.fasta'
#outputFile = 'input/input.txt' 

T = ""
for filename in glob.glob(inputFiles):
  fasta_sequences = SeqIO.parse(open(filename),'fasta')
  for fasta in fasta_sequences:
    T+=str(fasta.seq)+"+"

T+="|" # > Alphabet (used by l_intervals algorithm)

with open(outputFile, "w") as text_file:
    print(T, file=text_file)

