# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 11:59:48 2021

@author: aorti
"""

import statistics as stats

import matplotlib.pyplot as plt


from Bio import SeqIO
import numpy as np 
np.random.seed(444)


FILE="../data/controlseqs/uniprot-steroid+receptor-filtered-homo.fasta"
FILE="../data/controlseqs/uniprot-dehydrogenase-filtered-homo.fasta"
# =============================================================================
# FILE="../data/controlseqs/uniprot-acyltransferase-filtered-homo.fasta"
# =============================================================================
hist={}
lens=[]
filtered=[]

for record in SeqIO.parse(FILE, "fasta"):
    lens.append(len(record.seq))
med=stats.median(lens)
upper = med+0.35*med#stats.stdev(lens)
lower = med-0.35*med#stats.stdev(lens)/3
print(lower,upper)
print(stats.median(lens))
print(stats.median(lens)-stats.stdev(lens)/3)

n, bins, patches = plt.hist(x=lens, bins=50, color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.vlines(lower,0,50)
plt.vlines(upper,0,50)
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

for record in SeqIO.parse(FILE, "fasta"):
    if len(record.seq) > lower:
        if len(record.seq) < upper:
            record.description = record.description.replace(' ','.')
            filtered.append(record)
            
print('../data/controlseqs/filtered/2021-06-22_'+FILE[:-6]+'-formuscle.fasta')
with open('../data/controlseqs/filtered/2021-06-22_'+FILE[20:-6]+'-formuscle.fasta', "w") as output_handle:
    SeqIO.write(filtered, output_handle, "fasta")
