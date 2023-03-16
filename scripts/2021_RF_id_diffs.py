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
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc

from prettytable import PrettyTable


# open aln do random forest ################get pdist for location
# investigate gaps

# make matrix
# row taxons
# columns are locations of gaps

ALNPATH = "../data/SPIKEaln/PF19214_full.fasta"

ALNPATH = "../data/SPIKEaln/PF19209_full.fasta"
TREEPATH = "../data/SPIKEtree/clustered/2021-02-02_SPIKE09_MP.nwk"
N = 140
FAM = "SPIKE"

# =============================================================================
# ALNPATH = "../data/BSHaln/PF19209_full.fasta"
# TREEPATH = "../data/BSHtree/clustered/2021-02-02_SPIKE09_MP.nwk"
# N = 2100
# FAM = 'SPIKE'
# =============================================================================

ALNPATH = "../data/BSHaln/2021-01-15_BSH_filtered.meg"
TREEPATH = "../data/2021-02-02_BSH_filterd_ML.nwk"
FAM = "BSH"
N = 2100
HEIGHT = 22

# =============================================================================
# ALNPATH =  '../data/OMPaln/uniprot-gene_omp+reviewed_yes-8924.meg'
# TREEPATH = "../data/2021-02-02_ML_OMP.nwk"
# FAM = 'OMP'
# N=1650
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
    #    encode_splits=False,
    finish_node_fn=None,
    case_sensitive_taxon_labels=False,
    preserve_underscores=False,
    suppress_internal_node_taxa=True,
    suppress_leaf_node_taxa=False,
    terminating_semicolon_required=True,
    ignore_unrecognized_keyword_arguments=False,
)
treeOR.ladderize()
# treeOR.reroot_at_midpoint(update_bipartitions=False)
# print(pdm.as_data_table())

fileType = "a"
i = 0
taxons = []
catDict = {}
for taxon in treeOR.taxon_namespace:
    j = 0
    print(taxon)
    print("|".join(str(taxon)[1:-1].replace("...", ".").split("|")[:-1]))
    taxons.append(str(taxon)[1:-1].replace("...", "."))

    catDict["|".join(str(taxon)[1:-1].replace("...", ".").split("|")[:-1])] = str(
        taxon
    )[1:-1].replace("...", ".")

    if fileType == "fasta":
        catDict[str(taxon)[1:-1].replace(".", "_").split("|")[0]] = str(taxon)[
            1:-1
        ].replace(".", "_")
print(len(taxons))
print((catDict))


gapsDict = {}
seqCount = 0
taxa = []
alnDict = {}

gapsDict["clust"] = []
for i in range(N):
    gapsDict[i] = [0] * len(taxons)

if ALNPATH[-4:] == ".meg":
    with open(ALNPATH, "r") as reader:
        line = ""
        counter = 0
        txcount = -1
        for l in reader:
            l = l.strip()
            #  if seqCount != 0 & l[0] == >:
            #      gapsDict[]
            if len(l) == 0:
                gapCounter = 0
                cc = 0
                for char in line:
                    # print(cc,l)
                    cc += 1
                    # if len(gapsDict[counter]) <68:
                    # print(len(gapsDict[counter]))
                    if char == "-":
                        gapCounter += 1
                        gapsDict[counter][txcount] = gapCounter
                        # print(taxon,counter,gapCounter)
                    else:
                        gapCounter = 0
                        gapsDict[counter][txcount] = gapCounter
                    # else:
                    #   exit()
                    counter += 1
                pass

            elif l[0] == "#":
                print(l)
                if l.find("mega") == -1:
                    seqCount += 1
                    print(
                        l.strip()[1:]
                        .replace("_", ".")
                        .replace("(", "")
                        .replace(")", "")
                    )
                    taxon = catDict[
                        l.strip()[1:]
                        .replace("_", ".")
                        .replace("(", "")
                        .replace(")", "")
                        .replace("|org|", "|")
                        .replace("|org-type|Bacteria", "")
                    ]
                    gapsDict["clust"].append(taxon.split("|")[-1])
                    print(taxon, seqCount)
                    taxa.append(taxon)
                    alnDict[taxon] = ""
                    txcount += 1
                    counter = 0
                    line = ""
            else:
                line += l


seqCount = -1

if ALNPATH[-5:] == "fasta":
    with open(ALNPATH, "r") as reader:
        for l in reader:
            l = l.strip()
            #  if seqCount != 0 & l[0] == >:
            #      gapsDict[]

            if l[0] == ">":
                seqCount += 1
                taxon = catDict[l.strip()[1:].replace("_", ".")]
                gapsDict["clust"].append(taxon.split("|")[-1])
                print(taxon, seqCount)
                taxa.append(taxon)
                counter = 0
                line = ""
            else:
                print(l)
                line += l

            gapCounter = 0
            print(line)
            for char in line:
                # print(cc,l)
                # cc+=1
                # if len(gapsDict[counter]) <68:
                # print(len(gapsDict[counter]))
                if char == "-":
                    gapCounter += 1
                    gapsDict[counter][seqCount] = gapCounter
                    # print(taxon,counter,gapCounter)
                else:
                    gapCounter = 0
                    gapsDict[counter][seqCount] = gapCounter
                # else:
                #   exit()
                counter += 1
# =============================================================================
#             for char in line:
#                 if char == '-':
#                     gapCounter+=1
#                     gapsDict[counter].append(gapCounter)
#                    # print(counter,gapCounter)
#                 else:
#                     gapCounter = 0
#                     gapsDict[counter].append(gapCounter)
#                 counter += 1
# =============================================================================


d = {}
d["clust"] = gapsDict["clust"]
for k in gapsDict.keys():
    if k != "clust":
        if sum(gapsDict[k]) != 0:
            # print(k,len(gapsDict[k]))
            medianGap = statistics.median(gapsDict[k])
            d[str(k)] = gapsDict[k]
            for i in range(len(gapsDict[k])):
                d[str(k)][i] = int(d[str(k)][i])  # > medianGap)
            # d[str(k)] = gapsDict[k]
            if len(gapsDict[k]) != 184:
                print("woah")
                print(k, len(gapsDict[k]))

#   print(k,len(gapsDict[k]))#,gapsDict[k])

print(len(d["0"]), d)
print(len(taxa), taxa)
df = pd.DataFrame(data=d, index=taxa)
print(df)


# Split the Groups from the dataset where y is category and x is data with species
y = df.iloc[:, 0]
x = df.iloc[:, 1:]

print(y)
print(x)


# Split the data into training and test data for the categories(y) and dataset(x)
# Here we are spliting it 65% training and 35% test
X_train, X_test, y_train, y_test = train_test_split(
    x, y, test_size=0.35, random_state=42
)


ensemble_clfs = [
    (
        "RandomForestClassifier, max_features='sqrt'",
        RandomForestClassifier(
            warm_start=True, oob_score=True, max_features="sqrt", random_state=42
        ),
    ),
    (
        "RandomForestClassifier, max_features='log2'",
        RandomForestClassifier(
            warm_start=True, max_features="log2", oob_score=True, random_state=42
        ),
    ),
    (
        "RandomForestClassifier, max_features=None",
        RandomForestClassifier(
            warm_start=True, max_features=None, oob_score=True, random_state=42
        ),
    ),
]

# Map a classifier name to a list of (<n_estimators>, <error rate>) pairs.
error_rate = OrderedDict((label, []) for label, _ in ensemble_clfs)

# Range of `n_estimators` values to explore.
min_estimators = 50
max_estimators = 100

for label, clf in ensemble_clfs:
    for i in range(min_estimators, max_estimators + 1):
        clf.set_params(n_estimators=i)
        clf.fit(X_train, y_train)

        # Record the OOB error for each `n_estimators=i` setting.
        oob_error = 1 - clf.oob_score_
        error_rate[label].append((i, oob_error))

# Generate the "OOB error rate" vs. "n_estimators" plot.
for label, clf_err in error_rate.items():
    xs, ys = zip(*clf_err)
    plt.plot(xs, ys, label=label)

plt.xlim(min_estimators, max_estimators)
plt.xlabel("n_estimators")
plt.ylabel("OOB error rate")
plt.legend(loc="upper right")
plt.show()

clf = RandomForestClassifier(n_estimators=70, max_features=None, random_state=42)
all_accuracies = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=5)
print("Mean Validation Scores: ", end="")
print(np.mean(all_accuracies))

clf_final = RandomForestClassifier(
    n_estimators=70, bootstrap=True, max_features=None, oob_score=True, random_state=42
)
clf_final.fit(X_train, y_train)
y_pred = clf_final.predict(X_test)
print("Test Set Accuracy:", metrics.accuracy_score(y_test, y_pred))

# rf_probs = clf_final.predict_proba(X_test)[:, 1]
# roc_value = roc_auc_score(y_test, rf_probs)
# roc_value

print(clf_final.oob_score_)

feats = {}  # a dict to hold feature_name: feature_importance
for feature, importance in zip(df.columns, clf_final.feature_importances_):
    feats[feature] = importance  # add the name/value pair

importances = pd.DataFrame.from_dict(feats, orient="index").rename(
    columns={0: "Gini-importance"}
)
imp = importances.sort_values(by="Gini-importance", ascending=False)
print(imp.head(20))
print(list(imp.index)[:20])
alnImp = list(imp.index)[:5]
imp.head(6).to_csv("../data/RF/2021-03-01_RF_" + FAM + ".csv")


alnDict = {}

# gapsDict['clust'] = []
for i in range(len(taxons)):
    alnDict[i] = [0] * N


with open(ALNPATH, "r") as reader:
    for l in reader:
        # print(l.strip())
        line = l.strip()
        if len(line) == 0:
            pass
        elif line[0] == "#" or line[0] == ">":
            taxa = line[1:].replace("_", ".")
            #           print(taxa)
            alnDict[taxa] = ""
        else:
            alnDict[taxa] += line

c = 0
for imp in alnImp:
    alnImp[c] = int(imp)
    c += 1

alnImp.sort()
print(alnImp)
l = []
lDict = {}
for taxon, aln in alnDict.items():
    if type(taxon) == int or taxon == "mega":
        pass
    else:
        # print(len(aln))
        info = [
            ">"
            + catDict[taxon.strip().replace("_", ".").replace("(", "").replace(")", "")]
            .replace(".|org|", "|")
            .replace("|org-type|Bacteria", "")
            .replace("Streptococcus.faecium.", ".")
            .replace(".strain.DSM.17509./.CIP.109821./.100-23.", ".")
            .replace(".subsp..saprophyticus.", "")
            .replace(".strain.ATCC.BAA-472./.TX0016./.DO.", ""),
            aln[int(alnImp[0]) - 5 : int(alnImp[0]) + 5],
            aln[int(alnImp[1]) - 5 : int(alnImp[1]) + 5],
            aln[int(alnImp[2]) - 5 : int(alnImp[2]) + 5],
            aln[int(alnImp[3]) - 5 : int(alnImp[3]) + 5],
            aln[int(alnImp[4]) - 5 : int(alnImp[4]) + 5],
        ]
        l.append(info)
        print(
            ">"
            + catDict[
                taxon.strip().replace("_", ".").replace("(", "").replace(")", "")
            ],
            aln[int(alnImp[0]) - 5 : int(alnImp[0]) + 5],
            aln[int(alnImp[1]) - 5 : int(alnImp[1]) + 5],
            aln[int(alnImp[2]) - 5 : int(alnImp[2]) + 5],
            aln[int(alnImp[3]) - 5 : int(alnImp[3]) + 5],
            aln[int(alnImp[4]) - 5 : int(alnImp[4]) + 5],
        )
        # fig gini imp
        # string order by cluster number
# =============================================================================
#         ,
#               aln[int(alnImp[11])-5:int(alnImp[11])+5],
#               aln[int(alnImp[12])-5:int(alnImp[12])+5],
#               aln[int(alnImp[13])-5:int(alnImp[13])+5],
#               aln[int(alnImp[14])-5:int(alnImp[14])+5],
#               aln[int(alnImp[15])-5:int(alnImp[15])+5])
# =============================================================================
# =============================================================================
#         ,
#               aln[int(alnImp[16])-5:int(alnImp[16])+5],
#               aln[int(alnImp[17])-5:int(alnImp[17])+5],
#               aln[int(alnImp[18])-5:int(alnImp[18])+5],
#               aln[int(alnImp[19])-5:int(alnImp[19])+5],
#               )
# ============================================================================


# l = [["Hassan", 21, "LUMS"], ["Ali", 22, "FAST"], ["Ahmed", 23, "UET"]]

table = PrettyTable(["Name", alnImp[0], alnImp[1], alnImp[2], alnImp[3], alnImp[4]])

for rec in l:
    table.add_row(rec)

print(table)


from sklearn import svm

X = [[0, 0], [1, 1]]
y = [0, 1]
clf = svm.SVC()
clf.fit(X, y)
print(clf.predict([[2.0, 2.0]]))
