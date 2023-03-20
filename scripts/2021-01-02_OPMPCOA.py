# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 15:15:32 2020

@author: aorti
"""
import dendropy
import numpy as np
import scipy.spatial.distance as ssd
from sklearn.decomposition import PCA
from sklearn import preprocessing
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

TREEPATH = "../../PhyloClust/data/2020-04-06_MaxParsimony_BSH_idFull.nwk"


tree = dendropy.Tree.get(
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
pdm = tree.phylogenetic_distance_matrix()
print(pdm.as_data_table())


pdmA = np.zeros((len(tree.taxon_namespace), len(tree.taxon_namespace)))
counter = 0
print(pdmA.shape)
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon))
    for taxon2 in tree.taxon_namespace:
        pdmA[i, j] = float(pdm.distance(taxon, taxon2))
        j += 1
    i += 1

pdist = ssd.pdist(pdmA)

##setting data to final dataframe
df = pd.DataFrame(taxons, columns=["target"])

##center and scale PDM AND fit
pdmScaled = preprocessing.scale(pdmA)

# make PCA object
pca = PCA(0.975)
pca.fit(pdmA)
principalComponents = pca.transform(pdmA)

# scree plot
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ["PC" + str(x) for x in range(1, len(per_var) + 1)]

plt.bar(x=range(1, len(per_var) + 1), height=per_var, tick_label=labels)
plt.ylabel("Percntage of explained varience")
plt.xlabel("Principal Component")
plt.title("Scree Plot")
plt.show()


principalDf = pd.DataFrame(data=principalComponents, columns=labels)

finalDf = pd.concat([principalDf, df[["target"]]], axis=1)


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1)

ax.set_xlabel(
    "Principal Component 1: " + str(round(pca.explained_variance_ratio_[0], 2)),
    fontsize=15,
)
ax.set_ylabel(
    "Principal Component 2: " + str(round(pca.explained_variance_ratio_[1], 2)),
    fontsize=15,
)
ax.set_title("2 Component PCA", fontsize=20)

targets = ["OMP", "PMP", "PAL"]
colors = ["r", "g", "b"]
marks = ["o", "v", "*"]

for target, color, mark in zip(targets, colors, marks):
    indicesToKeep = finalDf["target"].str.contains(target)
    ax.scatter(
        finalDf.loc[indicesToKeep, "PC1"],
        finalDf.loc[indicesToKeep, "PC2"],
        c=color,
        s=50,
        marker=mark,
    )
ax.legend(targets)
ax.grid()
plt.show()


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel(
    "Principal Component 1: " + str(round(pca.explained_variance_ratio_[0], 2)),
    fontsize=15,
)
ax.set_ylabel(
    "Principal Component 2: " + str(round(pca.explained_variance_ratio_[1], 2)),
    fontsize=15,
)
ax.set_title("3 Component PCA", fontsize=20)

targets = ["OMP", "PMP", "PAL"]
colors = ["r", "g", "b"]
marks = ["o", "v", "*"]

for target, color, mark in zip(targets, colors, marks):
    indicesToKeep = finalDf["target"].str.contains(target)
    ax.scatter(
        finalDf.loc[indicesToKeep, "PC1"],
        finalDf.loc[indicesToKeep, "PC2"],
        finalDf.loc[indicesToKeep, "PC3"],
        c=color,
        s=50,
        marker=mark,
    )
ax.legend(targets)
ax.grid()
plt.show()
