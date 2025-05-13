# Python script for my pipeline
# This code should perform several phylogenomic analysis steps 

# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

# Store the path to the directory containing input files
in_dir="/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
print(in_dir)


# Store the path to a new folder where we want to write the output
out_dir="/scratch/forsythe/BB485/Week06/phy_pipeline_output/"
print(out_dir)


# Get a list of all fasta files in the input folder
all_in_files=glob.glob(in_dir+"*fasta")

# start a for-loop to loop through each file name, Inside of the loop:
    # create a mafft command for the file of interest
    # Call that mafft command from the command line using system call

for file in all_in_files[0:10]: # remember to remove this for the full run
    print(file)
    new_file_path = file.replace(in_dir, out_dir)
    print(new_file_path)

    # Create a command string (this is what get called using the 'system call'.
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    print(aln_cmd)

    #Run the command
    os.system(aln_cmd) #Uncomment this once you've double-checked that it's looking good.

#Create a for loop to run iqtree
# Get a list of all fasta files in the input folder
all_aln_files=glob.glob(out_dir+"*fasta")

# Loop through aln files and run iqtree on each (see example code and alter for my needs)
























































sys.exit()
#Directory where we're doing our work
#/scratch/forsythe/BB485_examples/Week06/phylogenomic_pipeline

# Store the paths to our input and output directories
indir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
outdir = "/scratch/forsythe/BB485_examples/Week06/phylogenomic_pipeline/alns_n_trees/"

# Create a list of all the fasta files in the indir
in_file_list=glob.glob(indir+"*fasta")

## TESTING: make the list shorter for now
in_file_list=in_file_list[0:3]

# Loop through items in the list and create alignments
for file in in_file_list:
    new_file_path = file.replace(indir, outdir)
    print(new_file_path)
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    print(aln_cmd)
    os.system(aln_cmd)


# Get a list of all the alignments
aln_list=glob.glob(outdir+"*fasta")

# Loop through all the alignments
for aln in aln_list:
    #tree_command = ["iqtree"+"-s"+aln+"-m", "TEST", "-nt 2"]
    tree_command = f"iqtree -s {aln} -m TEST -nt 2"
    print(tree_command)
    os.system(tree_command)

# Get a list of all the alignments
tree_list=glob.glob(outdir+"*.treefile")

#Create a list to store results
topo_results=[]

#loop through the tree files 
for tree in tree_list:
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")

    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #Stope the loop once we found the correct tip
            break
    
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
    
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
        
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
    

    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"
    
    #Add the results to your blank list
    topo_results.append(topo_str)

#Print the list
print(topo_results)








