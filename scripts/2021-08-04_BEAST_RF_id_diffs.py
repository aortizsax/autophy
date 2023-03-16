# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:37:51 2021

@author: aorti
"""

import dendropy

import statistics

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from collections import OrderedDict

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn import metrics

from prettytable import PrettyTable

from sklearn.linear_model import LogisticRegression
from yellowbrick.model_selection import FeatureImportances

# make matrix
# row taxons
# columns are locations of gaps


ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"

# TREEPATH = "../data/viraltrees/sarscov2sgly_muscle_x_cah.tree"
TREEPATH = "../data/viraltrees/clustered/2021-07-23_sarscov2sgly_umapemm.nwk"
N = 1300
FAM = "viral"
TREEID = "sarscov2sgly"
REP = "YP "

ALNPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle.fasta"
TREEPATH = "../data/BSHtree/clustered/2021-08-04_BSH_filterd_BEAST_umapemm.nwk"
N = 417
TREEID = "mouse"
FAM = "BSH"

# =============================================================================
# ALNPATH = '../data/pantheraln/2021-07-28_PTHR10000_aln_new.fasta'
# TREEPATH = "../data/panthertrees/2021-08-04_PTHR10000_aln.umapemm.tree"
# TREEID='PTHR1000'
# FAM = 'panther'
# N=700
# HEIGHT = 22
# =============================================================================


ALNPATH = "../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
TREEPATH = "../data/OMPtree/2021-08-04_OMPuniprot_muscle_umapemm.tree"
N = 1230
TREEID = "mb"
FAM = "OMP"
# =============================================================================
#
# ALNPATH = "../data/pfamaln/PF00005_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# FAM = 'pfam'
# HEIGHT = 22
#
# =============================================================================
# =============================================================================
# ALNPATH = "../data/pfamaln/PF00271_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00271_seed_umapemm.tree"
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# =============================================================================
# =============================================================================
# # # ALNPATH = "../data/pfamaln/PF00067_seed.txt"
# # # TREEPATH ="../data/pfamtrees/clustered/2021-07-23_PF00067_seed_umapemm.tree"
# # # TREEID='PF00067'
# # # FAM = 'pfam'
# # # HEIGHT = 22
# =============================================================================
# =============================================================================
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00501_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00501_seed_umapemm.tree"
# TREEID='PF00501'
# FAM = 'pfam'
# N=835
# HEIGHT = 22
#
# =============================================================================
# =============================================================================
#         ALNPATH = "../data/pfamaln/PF00528_seed.txt"
#         TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00528_seed_umapemm.tree"
#         TREEID='PF00528'
#         REP=' '
#         FAM = 'pfam'
#         HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00069_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00069_seed_umapemm.tree"
# TREEID='PF00069'
# REP=' '
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00072_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00072_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF07690_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"
# TREEPATH ="../data/controltrees/clustered/2021-07-23_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
# TREEID='acyltrans'
# N=1000
# FAM = 'control'
# HEIGHT = 22
# REP='.'
# =============================================================================

# =============================================================================
# ALNPATH = "../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
# TREEPATH = "../data/controltrees/clustered/2021-07-23_receptor-filtered-homo-muscle_umapemm.mean.tree"
# TREEID='receptor'
# N=1289
# FAM = 'control'
# HEIGHT = 22
# REP='.'
# =============================================================================

# =============================================================================
# ALNPATH= '../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta'
# TREEPATH = "../data/controltrees/clustered/2021-07-23_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.median.tree"
# TREEID='dehydrogenase'
# REP=' '
# FAM = 'control'
# HEIGHT = 22
# =============================================================================


treeOR = dendropy.Tree.get(
    path=TREEPATH,
    schema="newick",
    label=None,
    taxon_namespace=None,
    collection_offset=None,
    tree_offset=None,
    rooting="default-unrooted",
    edge_length_type=float,
    suppress_edge_lengths=False,
    extract_comment_metadata=True,
    store_tree_weights=False,
    finish_node_fn=None,
    case_sensitive_taxon_labels=False,
    preserve_underscores=False,
    suppress_internal_node_taxa=True,
    suppress_leaf_node_taxa=False,
    terminating_semicolon_required=True,
    ignore_unrecognized_keyword_arguments=False,
)
treeOR.ladderize()

fileType = "a"
i = 0
taxons = []
catDict = {}
for taxon in treeOR.taxon_namespace:
    j = 0
    print(taxon)
    print("|".join(str(taxon)[1:-1].split("|")[:]))
    print()
    taxons.append(str(taxon)[1:-1])

    if FAM == "control":
        catDict["|".join(str(taxon)[1:-1].split("|")[:-1])] = str(taxon)[1:-1]
    elif TREEID == "PF00528":
        catDict["|".join(str(taxon)[1:-1].split("|")[:])] = str(taxon)[1:-1]
    else:
        catDict["|".join(str(taxon)[1:-1].split("|")[:-1])] = str(taxon)[1:-1]


#'0': ['A', 'P', 'T'], '1': ['C', 'F', 'W', 'Y'], '2': ['D', 'E'], \
#'3': ['I', 'L', 'M', 'V'], '4': ['H', 'K'], '5': ['R', 'Q'],
#'6': ['N', 'G', 'S']

#'0': ['A', 'P', 'T'],
dict0 = {}
#'1': ['C', 'F', 'W', 'Y'],
dict1 = {}
#'2': ['D', 'E'],
dict2 = {}
#'3': ['I', 'L', 'M', 'V'],
dict3 = {}
#'4': ['H', 'K'],
dict4 = {}
#'5': ['R', 'Q'],
dict5 = {}
#'6': ['N', 'G', 'S']
dict6 = {}

gapsDict = {}
seqCount = -1
taxa = []
alnDict = {}

gapsDict["clust"] = []
for i in range(N + 1):
    gapsDict[i] = [0] * (len(taxons))
    dict0[str(i) + "APT"] = [0] * (len(taxons))
    dict1[str(i) + "CFWY"] = [0] * (len(taxons))
    dict2[str(i) + "DE"] = [0] * (len(taxons))
    dict3[str(i) + "ILMV"] = [0] * (len(taxons))
    dict4[str(i) + "HK"] = [0] * (len(taxons))
    dict5[str(i) + "RQ"] = [0] * (len(taxons))
    dict6[str(i) + "NGS"] = [0] * (len(taxons))


if (ALNPATH[-5:] == "fasta") | (ALNPATH[-3:] == "txt"):
    with open(ALNPATH, "r") as reader:
        for l in reader:
            l = l.strip()
            print(l)

            if l[0] == ">":
                seqCount += 1
                print(catDict)
                print(l)
                if (
                    (FAM == "control")
                    | (TREEID == "sarscov2sgly")
                    | (TREEID == "PF00528")
                ):
                    try:
                        taxon = catDict[l.strip()[1:].split(" ")[0].replace("_", REP)]
                    except:
                        taxon = catDict[l.strip()[1:].split(" ")[0].replace("YP", REP)]
                else:
                    print(l.strip()[1:].split(" ")[0])
                    taxon = catDict[l.strip()[1:].split(" ")[0]]
                gapsDict["clust"].append(taxon.split("|")[-1])
                print(taxon, seqCount)
                taxa.append(taxon)

                counter = 0
                line = ""
            else:
                print(l)
                line = l

            gapCounter = 0

            for char in line:
                if char == "-":
                    # if (char == 'H') | (char == 'K'):
                    gapCounter += 1
                    gapsDict[counter][seqCount] = gapCounter
                elif (char == "A") | (char == "P") | (char == "T"):
                    gapCounter = 0
                    dict0[str(counter) + "APT"][seqCount] = 1

                elif (char == "C") | (char == "F") | (char == "W") | (char == "Y"):
                    gapCounter = 0
                    dict1[str(counter) + "CFWY"][seqCount] = 1

                elif (char == "D") | (char == "E"):
                    gapCounter = 0
                    dict2[str(counter) + "DE"][seqCount] = 1

                elif (char == "I") | (char == "L") | (char == "V") | (char == "M"):
                    gapCounter = 0
                    dict3[str(counter) + "ILMV"][seqCount] = 1

                elif (char == "H") | (char == "K"):
                    gapCounter = 0
                    dict4[str(counter) + "HK"][seqCount] = 1

                elif (char == "R") | (char == "Q"):
                    gapCounter = 0
                    dict5[str(counter) + "RQ"][seqCount] = 1

                elif (char == "N") | (char == "G") | (char == "S"):
                    gapCounter = 0
                    dict6[str(counter) + "NGS"][seqCount] = 1

                else:
                    gapCounter = 0

                counter += 1


d = {}
d["clust"] = gapsDict["clust"]
l = np.asarray(d["clust"], dtype=np.float32, order="C")
print(max(l))
print(d["clust"])


maxclust = int(max(l)) + 1
print(maxclust)


found = {}
c = int(maxclust)
for i in range(len(gapsDict["clust"])):
    clusti = gapsDict["clust"][i]
    if int(float(gapsDict["clust"][i])) == round(float(gapsDict["clust"][i]), 1):
        d["clust"][i] = int(float(d["clust"][i]))
    else:
        if clusti in found.keys():
            print(clusti)
        else:
            found[clusti] = c
            c += 1
print()
for i in range(len(d["clust"])):
    clus = d["clust"][i]
    print(type(clus))
    if isinstance(clus, str):
        print(clus)
        d["clust"][i] = found[clus]

print(d["clust"])
print(found)


for k in gapsDict.keys():
    if k != "clust":
        if sum(gapsDict[k]) > 3:
            print(k, len(gapsDict[k]))
            medianGap = statistics.median(gapsDict[k])
            d[str(k)] = gapsDict[k]

            for i in range(len(gapsDict[k])):
                d[str(k)][i] = int(d[str(k)][i] > 0)  # medianGap)
            # d[str(k)] = gapsDict[k]
            if len(gapsDict[k]) != len(taxa):
                print("woah")
                print(k, len(gapsDict[k]))

d0 = {}
for k in dict0.keys():
    if k != "clust":
        if sum(dict0[k]) > 3:
            d0[str(k)] = dict0[k]

d1 = {}
for k in dict1.keys():
    if k != "clust":
        if sum(dict1[k]) > 3:
            d1[str(k)] = dict1[k]
d2 = {}
for k in dict2.keys():
    if k != "clust":
        if sum(dict2[k]) > 3:
            d2[str(k)] = dict2[k]

d3 = {}
for k in dict3.keys():
    if k != "clust":
        if sum(dict3[k]) > 3:
            d3[str(k)] = dict3[k]

d4 = {}
for k in dict4.keys():
    if k != "clust":
        if sum(dict4[k]) > 3:
            d4[str(k)] = dict4[k]

d5 = {}
for k in dict5.keys():
    if k != "clust":
        if sum(dict5[k]) > 3:
            d5[str(k)] = dict5[k]

d6 = {}
for k in dict6.keys():
    if k != "clust":
        if sum(dict6[k]) > 0:
            d6[str(k)] = dict6[k]

df = pd.DataFrame(data=d, index=taxa)
df0 = pd.DataFrame(data=d0, index=taxa)
df1 = pd.DataFrame(data=d1, index=taxa)
df2 = pd.DataFrame(data=d2, index=taxa)
df3 = pd.DataFrame(data=d3, index=taxa)
df4 = pd.DataFrame(data=d4, index=taxa)
df5 = pd.DataFrame(data=d5, index=taxa)
df6 = pd.DataFrame(data=d6, index=taxa)
df = df.merge(df0, left_index=True, right_index=True)
df = df.merge(df1, left_index=True, right_index=True)
df = df.merge(df2, left_index=True, right_index=True)
df = df.merge(df3, left_index=True, right_index=True)
df = df.merge(df4, left_index=True, right_index=True)
df = df.merge(df5, left_index=True, right_index=True)
df = df.merge(df6, left_index=True, right_index=True)


# Split the Groups from the dataset where y is category and x is data with species
y = df.iloc[:, 0]
x = df.iloc[:, 1:]

print(y)
print(x)

#######################################FEATURE IMPORTANCE
model = LogisticRegression(multi_class="auto", solver="liblinear")
fig = plt.figure(figsize=(8, 6))
viz = FeatureImportances(model, stack=True, relative=False, topn=20)
viz.fit(x, y)
viz.show()
fig.savefig(
    "../data/" + FAM + "colored/2021-08-04_" + FAM + TREEID + "_logimportance.tiff",
    dpi=300,
)

model = LogisticRegression(
    multi_class="auto", solver="liblinear", random_state=1111
).fit(x, y)

# print(model.coef_)
for i in range(len(model.coef_)):
    for j in range(len(model.coef_[i])):
        if abs(model.coef_[i][j]) > 0.06:
            print(i, df.columns.values[j], model.coef_[i][j])


modelB = RandomForestClassifier(n_estimators=10).fit(x, y)

impdict = []
index = []

for i in range(len(modelB.feature_importances_)):
    if modelB.feature_importances_[i] > 0.01:
        print(modelB.feature_importances_[i], df.columns.values[i])
        impdict.append([df.columns.values[i], modelB.feature_importances_[i]])
        index.append(df.columns.values[i])

impdf = pd.DataFrame(
    impdict, columns=["MSA Feature", "Relative Importance"], index=index
)

imphead = impdf.sort_values("Relative Importance")
print(imphead)

fig = plt.figure(figsize=(8, 6))
imphead.plot.barh()
fig.legend(loc="lower right")
fig.savefig(
    "../data/" + FAM + "colored/2021-08-04_" + FAM + TREEID + "_RFClassimportance.tiff",
    dpi=300,
)
fig.show()
