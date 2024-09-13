#! /usr/bin/python3
     
from Bio import SeqIO
import pandas as pd
import glob
import sys
import os

# arguments
arguments = len(sys.argv) - 1
if(arguments!=3): 
    quit()

data_group_name = os.path.basename(os.path.dirname(sys.argv[1]))  
inputFiles = os.path.join(sys.argv[1],'**/*.fasta')
outputFile = sys.argv[2]

#######################################
# 1. reads sequences and concatenates them
# while also numerating and creating uified protein dataset
####################################### 

# unified concatenated sequence
T = ""

# sequence dataset
sequence_table = []

# counter for sequence id assignment
seq_id = 0

for filename in glob.iglob(inputFiles,recursive=True):
  fasta_sequences = SeqIO.parse(open(filename),'fasta')
  family_name = os.path.basename(os.path.dirname(filename))
  if seq_id % 50000 == 0:
     print(seq_id)
  for fasta in fasta_sequences:
    sequence_name = str(fasta.name)
    sequence = str(fasta.seq)
    T+=sequence+"+"
    sequence_table.append([seq_id, family_name, sequence_name, sequence])
    seq_id +=1

T+="|" # > Alphabet (used by l_intervals algorithm)

# create unified concatenated sequence output
with open(outputFile, "w") as text_file:
    print(T, file=text_file)

# create sequence dataset output
sequence_dataset_output_path = os.path.join(sys.argv[3], data_group_name+"_sequence_dataset.csv")
sequence_df = pd.DataFrame(sequence_table, columns=["id", "family_name", "name", "sequence"])
sequence_df.to_csv(sequence_dataset_output_path, index=False)

