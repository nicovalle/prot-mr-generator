#! /usr/bin/python3
     
from Bio import SeqIO
import glob
import sys
import os

# arguments
arguments = len(sys.argv) - 1
if(arguments!=2): 
    quit()


inputFiles = os.path.join(sys.argv[1],'**/*.fasta')
outputFile = sys.argv[2]

#######################################
# 1. reads sequences and concatenates them
####################################### 

T = ""
for filename in glob.iglob(inputFiles,recursive=True):
  fasta_sequences = SeqIO.parse(open(filename),'fasta')
  for fasta in fasta_sequences:
    T+=str(fasta.seq)+"+"

T+="|" # > Alphabet (used by l_intervals algorithm)

with open(outputFile, "w") as text_file:
    print(T, file=text_file)

