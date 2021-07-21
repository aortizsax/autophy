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
np.random.seed(444)

start = 20
saveloci = '../data/controlaln/null/2021-06-28_'
typ = 'fasta'
end = 6
FILE="../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
FILE="../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta"
FILE="../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"


start = 26
end = -3
saveloci = '../data/BSHaln/null/2021-06-23_'
typ = 'nexus'
FILE="../data/BSHaln/2021-03-21_BSH_filtered_muscle.nex"

start = 26
end = -3
saveloci = '../data/OMPaln/null/2021-06-23_'
typ = 'fasta'
FILE="../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"

start = len('../data/viralaln/')
end = -5
saveloci = '../data/viralaln/null/2021-06-23_'
typ = 'fasta'
FILE="../data/viralaln/sarscov2gen_muscle.fasta"
FILE="../data/viralaln/2021-06-21_ebolasgly_muscle_trimmed.fasta"



hist={}
lens=[]
nulled=[]

for record in SeqIO.parse(FILE, typ):
    print(record)
    length = len(record.seq)
    record.seq = Seq('N' * length)
    nulled.append(record)
    
print(nulled)
    
print(saveloci+FILE[start:end]+'null.'+typ)
with open(saveloci+FILE[start:end]+'null.'+typ, "w") as output_handle:
    SeqIO.write(nulled, output_handle, typ)
