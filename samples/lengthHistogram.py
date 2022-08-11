from Bio import SeqIO

import sys

import seaborn as sns
import matplotlib.pyplot as plt


sns.set_style('darkgrid')

files_to_process = sys.argv[1:]

print(files_to_process)

all_seqs = []
all_lengths = []


import os


cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd)  # Get all the files in that directory
print("Files in %r: %s" % (cwd, files))

directory = os.fsencode(files_to_process[0])
    
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     with open(files_to_process[0] + filename, 'r') as f:
        fasta_sequences = SeqIO.parse(f, 'fasta' if filename.endswith('.fa') or filename.endswith('.fasta') else 'fastq')

        read_list = []
        length_list = []

        for fasta in fasta_sequences:
            read_list.append(fasta)
            length_list.append(len(fasta.seq))

        file_plot = sns.distplot(length_list, bins=10)
        file_plot.set_title(filename)

        all_seqs += read_list
        all_lengths += length_list

total_plot = sns.distplot(all_lengths)
total_plot.set_title('All Reads in sample.')
plt.show()

# for i in files_to_process:

    



    #     for read in reads:
    #         all_seqs.append(seq)

# print(all_seqs)
    
    # fasta_sequences = SeqIO.parse(open(i, 'r'),'fasta')

    # all_seqs += fasta_sequences
    
    # length_list = []

    # for fasta in fasta_sequences:
    #     name, sequence = fasta.id, str(fasta.seq)
    #     length_list.append(len(sequence))


    # print(length_list)

    # print([len(read.seq) for read in fasta_sequences])

    # sns.displot([len(read.seq) for read in fasta_sequences])






