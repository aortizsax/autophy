# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:05:17 2021

@author: aorti
"""

# this code imports Jones JTT substitution matrix
# has a litle code to HeirCluster for proof of concept
# START REAL CLUSTERING
# then UMAP dimensionality reduction
# then GMM-EM clustering while scoring each incrementing model
# GMM-EM is a k-means like with maxlikelihood scoring
# pick GMM-EM model based on lowest BIC score
# BIC is AIC like with number of observations in the fudge factor
# BIC: LN(obs) * number of parameters
# AIC:    2    * number of parameters
# plots desicion making progess
# prints labels for observations(row labels of orginal matrix)

import numpy as np
import itertools

from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture

from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
import sys
import dendropy

import numpy as np
import scipy.spatial.distance as ssd
from sklearn.decomposition import PCA
from sklearn import preprocessing
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import mixture
import itertools
from scipy import linalg
import matplotlib as mpl
from Bio import Phylo
import matplotlib
import re
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
from Bio.Align import substitution_matrices
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc


def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=25)
        faces.add_face_to_node(N, node, 0, position="branch-right")


N = 3

# umap distance metrics
metric = "manhattan"  #'euclidean'
metrics = [
    "manhattan",
    "minkowski",
    "canberra",
    "braycurtis",  #'haversine',
    "mahalanobis",
    "wminkowski",
    "cosine",
]

# fill in array with different number to use for loop structure below to test
# those values
n_Hs = [100]  # ,100,125,150,200]# [20,75,100,222]#[20,100,150,180,220,1000,2000,3000]


print(substitution_matrices.load())


matrix = substitution_matrices.load("JONES")
print(matrix)

print(matrix.alphabet)
blosum62 = np.zeros((len(matrix.alphabet), len(matrix.alphabet)))
for i in range(0, len(matrix.alphabet)):
    print(matrix.alphabet[i])
    ia = matrix.alphabet[i]
    for j in range(0, len(matrix.alphabet)):
        # print(matrix.alphabet[j])
        # ja=matrix.alphabet[j]
        # print(matrix[ia][ja])
        blosum62[i][j] = matrix[i][j]

print(blosum62)

# REPLACE blosum62 with your martix of observations(row) and features(col)
# or with your distance matrix
penguins = blosum62


## Start HeirClustering
linked = linkage(penguins, "single")
labelList = matrix.alphabet
plt.figure(figsize=(10, 7))
dendrogram(
    linked,
    orientation="top",
    labels=labelList,
    distance_sort="descending",
    show_leaf_counts=True,
)
plt.show()

plt.figure(figsize=(10, 7))
plt.title("Hierarchical Clustering Amino Acid Dendogram")
dend = shc.dendrogram(shc.linkage(penguins, method="weighted"), labels=labelList)
plt.show()
##END HeirClustering


# Start Real clustering procedure
for n_H in n_Hs:
    print("metric:", n_H)
    reducer = umap.UMAP(
        n_components=N, n_neighbors=n_H, min_dist=0.0, random_state=42, metric=metric
    )

    penguin_data = penguins
    scaled_penguin_data = StandardScaler().fit_transform(penguin_data)

    embedding = reducer.fit_transform(scaled_penguin_data)

    ###########PLOT ORIGNAL UMAP REDUCTION ON 2 DIMENSIONS
    plt.scatter(embedding[:, 0], embedding[:, 1])  # ,
    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")

    for i, txt in enumerate(list(matrix.alphabet)):
        print(txt, embedding[i, 0], embedding[i, 1])
        plt.annotate(txt, (float(embedding[i, 0]), float(embedding[i, 1])), size=15)
    plt.title("UMAP projection of the PDM", fontsize=24)
    plt.show()

    # AXIS LABELS
    labels = ["UMAP Coordinate" + str(x) for x in range(1, N + 1)]

    # ROWNAMES
    df = pd.DataFrame(list(matrix.alphabet), columns=["target"])

    # DF of UMAP reduction
    principalDf = pd.DataFrame(data=embedding[:, :], columns=labels)

    # MERGE DFs
    finalDf = pd.concat([principalDf, df[["target"]]], axis=1)

    #######################CLUTSER

    X = embedding[:, :]

    print(matrix.alphabet)
    lowest_bic = np.infty
    bic = []
    N_start = 4

    # Incrementing "N-means"
    n_components_range = range(N_start, 11)

    # different models
    cv_types = ["full"]  #'diag', 'tied',
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(
                n_components=n_components,
                random_state=111,
                covariance_type=cv_type,
                reg_covar=30e-06,
            )
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
    bic = np.array(bic)

    # color lists: first is colorblind friendly
    BPcolors = [
        "#999999",
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#f46d43",
        "#abd9e9",
    ]  # COLORBLIND FRIENDLY

    BPcolors = [
        "#e6194B",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#f58231",
        "#ba3fdf",
        "#42d4f4",
        "#bfef45",
        "#fabed4",
        "#469990",
        "#dcbeff",
        "#9A6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#0072B2",
    ]
    color_iter = itertools.cycle(BPcolors)

    # set model to use
    clf = best_gmm
    bars = []

    ########## Plot the BIC scores
    plt.figure(figsize=(6, 6))
    spl = plt.subplot(2, 1, 1)
    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        xpos = np.array(n_components_range) + 0.2 * (i - 2)
        bars.append(
            plt.bar(
                xpos,
                bic[i * len(n_components_range) : (i + 1) * len(n_components_range)],
                width=0.2,
                color=color,
            )
        )
    plt.xticks(n_components_range)
    plt.ylim([bic.min() * 1.01 - 0.01 * bic.max(), bic.max()])
    plt.title("BIC score per model")
    xpos = (
        np.mod(bic.argmin(), len(n_components_range))
        + 0.65
        + 0.2 * np.floor(bic.argmin() / len(n_components_range))
    )
    plt.text(xpos + N_start - 1, bic.min() * 0.97 + 0.03 * bic.max(), "*", fontsize=14)
    spl.set_xlabel("Number of Gaussian components")
    spl.legend([b[0] for b in bars], cv_types)
    plt.ylabel("BIC score")

    ######### Plot the winner
    color_iter = itertools.cycle(BPcolors)

    splot = plt.subplot(2, 1, 2)
    Y_ = clf.predict(X)
    print(clf.covariances_.shape, clf.covariances_, clf.means_)
    for i, (mean, cov, color) in enumerate(
        zip(clf.means_, clf.covariances_, color_iter)
    ):
        if not np.any(Y_ == i):
            continue

        print(i, mean, color)
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color, marker="o")

    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")

    plt.xlim(min(embedding[:, 0]) - 0.2, max(embedding[:, 0]) + 0.2)
    plt.ylim(min(embedding[:, 1]) - 0.2, max(embedding[:, 1]) + 0.2)
    for i, txt in enumerate(list(matrix.alphabet)):
        print(txt, embedding[i, 0], embedding[i, 1])
        plt.annotate(txt, (float(embedding[i, 0]), float(embedding[i, 1])), size=15)
    plt.title(
        "Selected GMM: full model, " + str(clf.n_components) + " Gaussian components"
    )
    plt.subplots_adjust(hspace=0.35, bottom=0.02)
    plt.savefig(
        "../data/EMMplots/2021-06-30_" + str(n_H) + "_GMMsweep_" + metric + ".png"
    )
    plt.show()

    dictpred = {
        "0": [],
        "1": [],
        "2": [],
        "3": [],
        "4": [],
        "5": [],
        "6": [],
        "7": [],
        "8": [],
        "9": [],
        "10": [],
        "11": [],
        "12": [],
        "13": [],
        "14": [],
    }
    for i in range(len(matrix.alphabet)):
        # print(matrix.alphabet[i],clf.predict(X)[i])
        dictpred[str(clf.predict(X)[i])].append(matrix.alphabet[i])

    print(dictpred)
