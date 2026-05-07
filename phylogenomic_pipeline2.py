# Import needed python modules
# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo


# Create a string variable that stores the folder path to the input files
file_name = "At_4G09350.fasta"

input_folder_path = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"

full_path_to_file = input_folder_path + file_name

print(full_path_to_file)

# Perform mulitple sequence alignment using MAFFt


# Perform phylogenetic tree construction using IQtree


# Perform tree topology test with biopython


# store output
