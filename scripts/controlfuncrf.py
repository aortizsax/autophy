# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 11:59:48 2021

@author: aorti
"""

import pprint

import statistics as stats

import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

np.random.seed(444)

start = 20
saveloci = "../data/controlaln/null/2021-07-15_"
typ = "fasta"
end = 6
FILE = "../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
INF = "../data/controlcolored/2021-07-15controlreceptor_colored.csv"
INRF = "../data/controlcolored/2021-07-15controlreceptor_RF.csv"


FILE = "../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"
INF = "../data/controlcolored/2021-07-15controlacyltrans_colored.csv"
INRF = "../data/controlcolored/2021-07-15controlacyltrans_RF.csv"

# =============================================================================
# FILE="../data/controlaln/2021-06-22_dehydrogenase-filtered-homo-muscle.fasta"
# INF="../data/controlcolored/2021-07-15controldehydrogenase_colored.csv"
# INRF='../data/controlcolored/2021-07-15controldehydrogenase_RF.csv'
#
# =============================================================================

# =============================================================================
#
# start = 26
# end = -3
# saveloci = '../data/BSHaln/null/2021-06-23_'
# typ = 'nexus'
# FILE="../data/BSHaln/2021-03-21_BSH_filtered_muscle.nex"
#
# start = 26
# end = -3
# saveloci = '../data/OMPaln/null/2021-06-23_'
# typ = 'fasta'
# FILE="../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
#
# start = len('../data/viralaln/')
# end = -5
# saveloci = '../data/viralaln/null/2021-06-23_'
# typ = 'fasta'
# FILE="../data/viralaln/sarscov2gen_muscle.fasta"
# FILE="../data/viralaln/2021-06-21_ebolasgly_muscle_trimmed.fasta"
#
# =============================================================================


hist = {}
lens = []
idfunc = {}

for record in SeqIO.parse(FILE, typ):
    length = len(record.seq)
    ID = record.id
    des = record.description.split("OS")[0]
    print(des.split(" ")[-1])
    idfunc[ID] = des.split(" ")[-1]


pprint.pprint(idfunc)

import csv

inrf = {}
with open(INRF, newline="\n") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for row in reader:
        print(row)
        inrf["|".join(row[0].split("|")[:-1])] = row[1:]
print(inrf)
with open(INF, newline="\n") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for row in reader:
        print(row, idfunc[row[1].replace(".", "_")])
        # print(inrf['>'+row[1].replace('_', ' ')][:5])

print(saveloci + FILE[start:end] + "null." + typ)
# =============================================================================
# with open(saveloci+FILE[start:end]+'null.'+typ, "w") as output_handle:
#     SeqIO.write(nulled, output_handle, typ)
# =============================================================================

annarray = {}
for record in SeqIO.parse(FILE, typ):
    anno = str(record.description.split("|")[-1][:3])
    annarray[anno] = 0


for record in SeqIO.parse(FILE, typ):
    anno = record.description.split("|")[-1][:3]
    annarray[anno] += 1
print(annarray)
print(len(annarray))

N = 0
for k, v in annarray.items():
    if v >= 2:
        N += 1
print(N)
