# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 11:59:48 2021

@author: aorti
"""

import statistics as stats

import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np 
import statistics



np.random.seed(444)

start = 20
saveloci = '../data/controlaln/null/2021-06-28_'
typ = 'fasta'
end = 6
FILE="../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
# =============================================================================
# FILE="../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta"
# =============================================================================
# =============================================================================
# FILE="../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"
# =============================================================================

# =============================================================================
# 
# start = 26
# end = -3
# saveloci = '../data/BSHaln/null/2021-06-23_'
# typ = 'nexus'
# FILE="../data/BSHaln/2021-03-21_BSH_filtered_muscle.nex"
# 
typ = 'fasta'
FILE= "../data/viralaln/2021-06-18_sarscov2sgly_muscle.fasta"
# =============================================================================
# FILE= "../data/viralaln/2021-06-21_ebolasgly_muscle_trimmed.fasta"
# =============================================================================
# 
# typ = 'fasta'
# FILE= "../data/pfamaln/PF00067_seed.txt"
# FILE= "../data/pfamaln/PF00005_seed.txt"
# FILE= "../data/pfamaln/PF00069_seed.txt"
# FILE= "../data/pfamaln/PF00072_seed.txt"
# FILE= "../data/pfamaln/PF00271_seed.txt"
# FILE= "../data/pfamaln/PF00528_seed.txt"
# FILE= "../data/pfamaln/PF00501_seed.txt"
# FILE= "../data/pfamaln/PF02518_seed.txt"
# FILE= "../data/pfamaln/PF07690_seed.txt"
# 
# =============================================================================
# =============================================================================
# start = 26
# end = -3
# saveloci = '../data/OMPaln/null/2021-06-23_'
# typ = 'fasta'
# FILE="../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
# =============================================================================

# =============================================================================
# start = len('../data/viralaln/')
# end = -5
# saveloci = '../data/viralaln/null/2021-06-23_'
# typ = 'fasta'
# FILE="../data/viralaln/sarscov2gen_muscle.fasta"
# FILE="../data/viralaln/2021-06-21_ebolasgly_muscle_trimmed.fasta"
# 
# =============================================================================


mean = 0
X = []
N = 0
L = 0
parsed = SeqIO.parse(FILE, typ)


for record in parsed:
    for rec in parsed:
        recmatch = 0
        N+=1
        L=0
        length = 0
        for i in range(0,len(record.seq)):
            length +=1
            if (record.seq[i] != '-'):
                if (rec.seq[i] != '-'):
                    L+=1
                    if record.seq[i] == rec.seq[i]: 
                        #print(record.seq[i],rec.seq[i])
                        recmatch += 1
        #print(recmatch/L,L)
        mean += recmatch/L
        X.append(recmatch/L)
        
print(statistics.mean(X))
print(statistics.stdev(X))
print(N)
print('Length: ',L)
print('Length: ',length)
                    
                    
            
            



# =============================================================================
#     print(record)
#     length = len(record.seq)
#     record.seq = Seq('N' * length)
#     nulled.append(record)
#     
# print(nulled)
#     
# print(saveloci+FILE[start:end]+'null.'+typ)
# with open(saveloci+FILE[start:end]+'null.'+typ, "w") as output_handle:
#     SeqIO.write(nulled, output_handle, typ)
# 
# =============================================================================
