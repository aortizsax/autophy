"""
Created on Thu Jan  7 16:13:27 2021

@author: aorti

Description

Usage

Arguments


"""
from matplotlib.ticker import StrMethodFormatter


import numpy as np
import itertools
import re

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
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
import pprint
import time


from datetime import date
import datetime

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc

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

# make matrix
# row taxons
# columns are locations of gaps


ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"

TREEPATH = "../data/viraltrees/sarscov2sgly_muscle_x_cah.tree"
TREEPATH = "../data/viraltrees/clustered/2021-07-23_sarscov2sgly_umapemm.nwk"
N = 1600
FAM = "viral"
TREEID = "sarscov2sgly"
REP = "YP "

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


ALNPATH = "../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
TREEPATH = (
    "../data/controltrees/2021-06-22_receptor-filtered-homo-musclecorr.median.tree"
)
OUTTREE = "../data/controltrees/clustered/2021-07-23_receptor-filtered-homo-muscle_umapemm.mean.tree"
TREEID = "receptor"
N = 1289
FAM = "control"
HEIGHT = 22
REP = "."


def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=25)
        faces.add_face_to_node(N, node, 0, position="branch-right")


# =============================================================================
# ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"
# ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
# TREEPATH = "../data/viraltrees/2021-06-23_sarscov2sgly_muscle_x.median.tree"
# OUTTREE = "../data/viraltrees/clustered/2021-10-29_sarscov2sgly_umapemm.nwk"
# TREEID='sarscov2sgly'
# Naln=2000
# FAM = 'viral'
# HEIGHT = 22
#
#
# ALNPATH = "../data/viralaln/2021-07-27b_ebolaspikeglycoprotein_muscle.fasta"
# TREEPATH = "../data/viraltrees/2021-08-09b_ebolaspikeglycoprotein_muscle_new.11112222.mean.tree"
# OUTTREE = "../data/viraltrees/clustered/2021-10-29_ebolaspikeglycoprotein_muscle_new.umap.tree"
# TREEID='ebola'
# FAM = 'viral'
# Naln=1000
# HEIGHT = 22
# =============================================================================

# =============================================================================
# =============================================================================
# # TREEPATH = "../data/SPIKEtree/2021-01-24_SPIKE09_btsp50_ML.nwk"
# # OUTTREE = "../data/SPIKEtree/clustered/2021-10-29_SPIKE09_ML.nwk"
# # FAM = 'SPIKE'
# # HEIGHT = 22
# # =============================================================================
#
#
# # =============================================================================
# # TREEPATH = "../../Taxon_Enzyme_Phylogeny/BSH_data/2020-03-05_MaxParsimony_BSH_id.nwk"
# =============================================================================
# =============================================================================

# =============================================================================
# ALNPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle.fasta"
# TREEPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle_relaxed.tree"
# OUTTREE = "../data/BSHtree/clustered/2021-10-29_BSH_filterd_BEAST_umapemm.nwk"
# Naln=417
# TREEID='mouse'
# FAM = 'BSH'
# HEIGHT = 22
#
# ALNPATH = "../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
# TREEPATH = "../data/OMPtree/2021-06-16_OMPuniprot_muscle.tree"
# OUTTREE = "../data/OMPtree/2021-10-29_OMPuniprot_muscle_umapemm.tree"
# Naln=1230
# TREEID='uniprot'
# FAM = 'OMP'
# HEIGHT = 22
#
# ALNPATH = '../data/pantheraln/2021-07-28_PTHR10000_aln_new.fasta'
# TREEPATH = "../data/panthertrees/2021-10-29_PTHR10000_aln.umapemm.tree"
# TREEPATH = "../data/panthertrees/2021-08-02_PTHR10000_aln.2222.mean.tree"
# OUTTREE = "../data/panthertrees/2021-10-29_PTHR10000_aln.umapemm.tree"
# TREEID='PTHR10000'
# FAM = 'panther'
# Naln=700
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = '../data/pantheraln/2021-07-28_PTHR10010_aln.fasta'
# TREEPATH = "../data/panthertrees/2021-10-29_PTHR10010_aln.umapemm.tree"
# TREEPATH = "../data/panthertrees/2021-08-09_PTHR10010_aln.mean.2222.tree"
# OUTTREE = "../data/panthertrees/2021-10-29_PTHR10010_aln.umapemm.tree"
# TREEID='PTHR10010'
# FAM = 'panther'
# Naln=1116
# HEIGHT = 22
# =============================================================================
# =============================================================================
# # check
# ALNPATH = '../data/pantheraln/2021-07-28_PTHR10005_aln.fasta'
# TREEPATH = "../data/panthertrees/2021-10-29_PTHR10005_aln.umapemm.tree"
# TREEPATH = "../data/panthertrees/2021-08-09_PTHR10005_aln.mean.2222.tree"
# OUTTREE = "../data/panthertrees/2021-10-29_PTHR10005_aln.umapemm.tree"
# TREEID='PTHR10005'
# FAM = 'panther'
# Naln=2000
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = '../data/pantheraln/2021-08-02_PTHR100037_aln.fasta'
# TREEPATH = "../data/panthertrees/2021-10-29_PTHR10037_aln.umapemm.tree"
# TREEPATH = "../data/panthertrees/2021-08-09_PTHR10037_aln.mean.2222.tree"
# OUTTREE = "../data/panthertrees/2021-10-29_PTHR10037_aln.umapemm.tree"
# TREEID='PTHR10037'
# FAM = 'panther'
# Naln=2000
# HEIGHT = 22
# =============================================================================


# =============================================================================
# ALNPATH = "../data/pfamaln/PF00005_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# TREEPATH = "../data/pfamtrees/2021-06-14_PF00005_seed_strctclk.tree"
# OUTTREE = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# TREEID='PF00005'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00271_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-06-18_PF00271_seed.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00271_seed_umapemm.tree"
# TREEID='PF00271'
# FAM = 'pfam'
# HEIGHT = 22
#
#
# ALNPATH = "../data/pfamaln/PF07690_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-06-18_PF07690_seed.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# ###seccatTax
# ALNPATH = "../data/pfamaln/PF00067_seed.txt"
# TREEPATH ="../data/pfamtrees/clustered/2021-07-23_PF00067_seed_umapemm.tree"
# TREEPATH = "../data/pfamtrees/2021-07-02_PF00067_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00067_seed_umapemm.tree"
# TREEID='PF00067'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
#
# ALNPATH = "../data/pfamaln/PF00501_seed.txt"
# Naln=835
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00501_seed_umapemm.tree"
# TREEID='PF00501'
# FAM = 'pfam'
# HEIGHT = 22
#
#
# =============================================================================
# =============================================================================
# ALNPATH = "../data/pfamaln/PF00528_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00528_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00528_seed_umapemm.tree"
# REP=' '
# TREEID='PF00528'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ##fix arrary seccattax
# ALNPATH = "../data/pfamaln/PF00069_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-23_PF00069_seed_umapemm.tree"
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00069_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00069_seed_umapemm.tree"
# TREEID='PF00069'
# FAM = 'pfam'
# REP=' '
# HEIGHT = 22
# =============================================================================


# =============================================================================
# ALNPATH = "../data/pfamaln/PF00072_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00072_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF00501_seed_umapemm.tree"
# TREEID='PF00501'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF02518_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-07-15_PF02518_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF02518_seed_umapemm.tree"
# TREEID='PF02518'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# ALNPATH = "../data/pfamaln/PF07690_seed.txt"
# TREEPATH = "../data/pfamtrees/2021-06-18_PF07690_seed.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-10-29_PF07690_seedumapgmmhc.median.tree"
# TREEID='PF02518'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
ALNPATH = "../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"
TREEPATH = "../data/controltrees/2021-06-22_acyltransferase-filtered-homo-musclecorr.1111.median.tree"
OUTTREE = "../data/controltrees/clustered/2021-10-29_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
TREEID = "acyltrans"
FAM = "control"
HEIGHT = 22
Naln = 1000
REP = "."
# =============================================================================


ALNPATH = "../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
TREEPATH = "../data/controltrees/2021-06-22_receptor-filtered-homo-musclecorr.mean.tree"
OUTTREE = "../data/controltrees/clustered/2021-10-29_receptor-filtered-homo-muscle_umapemm.mean.tree"
TREEID = "receptor"
FAM = "control"
HEIGHT = 22
Naln = 1289
REP = "."

ALNPATH = (
    "../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta"
)
TREEPATH = "../data/controltrees/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.1111.mean.tree"
OUTTREE = "../data/controltrees/clustered/2021-10-29_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.mean.tree"
TREEID = "dehydrogenase"
FAM = "control"
REP = " "
HEIGHT = 22


# =============================================================================
# TREEPATH = "../data/nextstraintrees/nextstrain_flu_seasonal_h3n2_ha_2y_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_flu_seasonal_h3n2_ha_2y_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# TREEID='FLU'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/nextstraintrees/nextstrain_mumps_na_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_mumps_na_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# HEIGHT = 22
# TREEID="mumps"
# =============================================================================

# =============================================================================
# ALNPATH= '../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta'
# TREEPATH = "../data/controltrees/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-10-29_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.mean.tree"
# TREEID='dehydrogenase'
# REP=' '
# FAM = 'control'
# HEIGHT = 22
# =============================================================================
# =============================================================================
# TREEPATH = "../data/nextstraintrees/nextstrain_ncov_global_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_ncov_global_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# TREEID='ncov'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../data/nextstraintrees/2021-08-11_nextstrain_ncov_gisaid_global_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/2021-10-29_ncov_.nexus"
# FAM = 'nextstrain'
# TREEID='ncov'
# HEIGHT = 22
# =============================================================================
ALNPATH = "../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
TREEPATH = "../data/uniprotbeast/uniprotOMP_muscle_1647964589147_cah.tree"
OUTTREE = "../data/uniprottrees/uniprotOMP_muscle_cah_umapemm.tree"
Naln = 1230
TREEID = "OMP"  #'uniprot'
FAM = "uniprot"  #'OMP'
HEIGHT = 22


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
]

if (
    (FAM == "OMP")
    | (FAM == "nextstrain")
    | (FAM == "pfam")
    | (FAM == "control")
    | (FAM == "viral")
):
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
        "#b3a400",
        "#800000",
        "#82baff",
        "#808000",
        "#ffd8b1",
        "#0072B2",
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
    ]

if (TREEPATH[-4:] == ".nhx") | (TREEPATH[-4:] == ".nwk"):
    tree = dendropy.Tree.get(
        path=TREEPATH,
        schema="newick",
        label=None,
        taxon_namespace=None,
        collection_offset=None,
        tree_offset=None,
        rooting="default-rooted",
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

else:
    tree = dendropy.Tree.get(
        path=TREEPATH,
        schema="nexus",
        rooting="default-rooted",
        edge_length_type=float,
        case_sensitive_taxon_labels=False,
        # =============================================================================
        #         suppress_internal_node_taxa=True,
        # =============================================================================
        terminating_semicolon_required=True,
    )

tree.ladderize()
pdm = tree.phylogenetic_distance_matrix()

start_time = time.time()

pdmA = np.zeros((len(tree.taxon_namespace), len(tree.taxon_namespace)))
counter = 0
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon)[1:-1])
    for taxon2 in tree.taxon_namespace:
        if j > i:
            pdmA[i, j] = float(pdm.distance(taxon, taxon2))
            pdmA[j, i] = float(pdm.distance(taxon, taxon2))
            # if negative branch legths
        # =============================================================================
        #             if pdmA[i,j] < 0:
        #                 pdmA[i,j] = abs(pdmA[i,j])
        #                 pdmA[j,i] = abs(pdmA[j,i])
        # =============================================================================

        j += 1
    i += 1


print("--- %s seconds ---" % (time.time() - start_time))

# dimensions of reduciton
N = 2

if FAM == "PRACT":
    N = 5


metric = "precomputed"

# n neighbors like tsneq
n_Hs = [12]  # int(len(taxons)*3/4),


for n_H in n_Hs:
    print("metric:", metric)
    reducer = umap.UMAP(
        n_components=N,
        min_dist=0.0,
        n_neighbors=n_H,
        random_state=222,
        metric="precomputed",
    )

    embedding = reducer.fit_transform(pdmA)

    labels = ["UMAP-C" + str(x) for x in range(0, N)]

    df = pd.DataFrame(taxons, columns=["target"])

    principalDf = pd.DataFrame(data=embedding[:, :], columns=labels)

    finalDf = pd.concat([principalDf, df[["target"]]], axis=1)

    plt.figure(figsize=(7, 3.9))
    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.title("UMAP projection of the phylogenetic dataset")
    plt.savefig(
        "../data/EMMplots/2021-10-29_"
        + str(n_H)
        + "_ogUMAPprog"
        + FAM
        + TREEID
        + "_"
        + metric
        + ".svg",
        dpi=600,
        format="svg",
    )

    plt.show()
    #######################CLUTSER

    X = embedding[:, :]

    start = 1

    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 30)
    if FAM == "control":
        n_components_range = range(1, 50)
    if FAM == "OMP":
        start = 3
        n_components_range = range(start, 25)
    if FAM == "uniprot":
        start = 1
        n_components_range = range(start, 20)
    if (FAM == "viral") | (FAM == "nextstrain") | (FAM == "viral"):
        start = 8
        n_components_range = range(start, 35)
    if FAM == "PRACT":
        n_components_range = range(1, 6)
    cv_types = ["full"]  # ,'tied']

    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(
                n_components=n_components,
                random_state=111,
                covariance_type=cv_type,
                reg_covar=10e-06,
            )
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic + 10:
                lowest_bic = bic[-1]
                best_gmm = gmm
    print("--- %s seconds ---" % (time.time() - start_time))
    clf = best_gmm
    bars = []

    color_iter = itertools.cycle(["navy", "turquoise", "cornflowerblue", "darkorange"])
    bic = np.array(bic)

    # Plot the BIC scores
    plt.figure(figsize=(8, 6))
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
    plt.text(xpos + 0.5, bic.min() * 0.97 + 0.03 * bic.max(), "*", fontsize=14)
    spl.set_xlabel("Number of components")
    spl.legend([b[0] for b in bars], cv_types)

    # Plot the winner
    splot = plt.subplot(2, 1, 2)
    Y_ = clf.predict(X)
    for i, (mean, cov, color) in enumerate(
        zip(clf.means_, clf.covariances_, color_iter)
    ):
        v, w = linalg.eigh(cov)
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 0.8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan2(w[0][1], w[0][0])
        angle = 180.0 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180.0 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

    plt.xticks(())
    plt.yticks(())
    plt.title(
        f"Selected GMM: {best_gmm.covariance_type} model, "
        f"{best_gmm.n_components} components"
    )
    plt.subplots_adjust(hspace=0.35, bottom=0.02)
    plt.show()

    ##############################################################
    taxonsLabels = []
    # for tree label replacement
    repDict = {}
    # for coloring tree
    catTax = []

    from scipy.stats import rankdata

    clfpredictX = clf.predict(X)
    # rankdata(list(clf.predict(X)))

    print(clf.predict(X), clfpredictX)
    a_set = set(clfpredictX)
    number_of_unique_values = len(a_set)
    print(number_of_unique_values)

    for i in range(number_of_unique_values):
        catTax.append([])
    for i in range(0, len(taxons)):
        taxons[i] = taxons[i].replace(" ", "_")
        taxonsLabels.append(taxons[i] + "|" + str(clfpredictX[i]))

        # for coloring tree
        catTax[int(clfpredictX[i])].append(taxons[i] + "|" + str(clfpredictX[i]))

        # for tree label replacement
        repDict[taxons[i]] = taxons[i] + "|" + str(clfpredictX[i])

    # tree as string
    strTree = tree.as_string(schema="newick")

    # string replacement
    strTree = re.sub("Inner[0-9][0-9][0-9]", "", strTree)
    strTree = re.sub("Inner[0-9][0-9]", "", strTree)
    strTree = re.sub("Inner[0-9]", "", strTree)

    # replace tips
    for key in repDict.items():
        strTree = strTree.replace(key[0], key[1])

    # make copy to edit one to add colon
    finalStrTree = strTree
    for iter in re.finditer("\D(?<=:)", strTree):
        a = iter.span()
        b = a[1]
        a = a[0]
        if strTree[a - 7] == ")":
            finalStrTree = finalStrTree.replace(strTree[a - 7 : b], "):")
    # print(finalStrTree)

    # save tree
    EMtree_file = open(OUTTREE, "w+")
    n = EMtree_file.write(finalStrTree)
    EMtree_file.close()

    #########################################################################
    ###########CHECK MONOPHYLY
    from Bio.Phylo.PhyloXML import Phylogeny

    BPtree = Phylo.read(OUTTREE, "newick")

    terminals = BPtree.get_terminals()
    termdict = {}
    flag = []
    flagpp = []

    for terminal in terminals:
        termdict[str(terminal)] = terminal
        print(str(terminal))
    print(termdict)
    for v in catTax:
        v = v.copy()
        i = 0
        # set length ofvector
        taxa = [0] * len(v)

        for tax in v:
            # iflen taxa label > 40
            if len(tax) > 40:
                taxa[i] = termdict[tax[:37] + "..."]
            else:
                taxa[i] = termdict[tax]
            i += 1

        # MPhy test
        mp = BPtree.is_monophyletic(taxa)

        if (mp) == False:
            flag.append(v)
            flagpp.append(taxa)
    print("--- %s START to Check monophylo seconds ---" % (time.time() - start_time))
    print(flag)
    print()
    print(catTax)

    tree = dendropy.Tree.get(
        path=OUTTREE,
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
    tree.encode_bipartitions()
    print(flag)
    print()
    groups = []

    mppassdict = {}

    tmpflag = flag
    reclustCount = 0
    reclustlist = {}
    n_comp = clf.n_components

    for f in flag:
        reclustCount += 1

        tmpf = f.copy()
        subclust = 0
        y_mp = [0] * len(tmpf)
        clu = f[0].split("|")[-1]

        while len(f) > 0:
            maxdescend = 0
            poptax = []
            for edge in tree.levelorder_edge_iter():
                descendants = edge.bipartition.leafset_as_newick_string(
                    tree.taxon_namespace
                ).split(", (")
                for descend in [descendants[0]]:
                    descend = descend.replace("(", "")
                    descend = descend.replace(")", "")
                    descend = descend.replace(";", "")
                    descend = descend.replace("'", "")

                    descend = descend.split(",")
                    for i in range(len(descend)):
                        descend[i] = descend[i].replace(" ", "")

                    res = len([w for w in descend if w in f])
                    if res == len(descend):
                        tempdescend = len(descend)

                        if maxdescend < tempdescend:
                            print(res, f, descend)

                            maxdescend = tempdescend
                            monophy = descend
            # print(maxdescend,monophy,f)
            # print('working',f,tmpf)
            for x in monophy:
                ogindex = tmpf.index(x)
                # print('siubclust',subclust,len(y_mp))
                # print(f.index(x),tmpf.index(x))
                print(f.pop(f.index(x)))
                y_mp[ogindex] = subclust
            subclust = subclust + 1
        print(tmpf, y_mp)
        secondpass = {}
        maxcluster = 0

        for i in range(len(y_mp)):
            tempcol = y_mp[i]
            # print(tempcol,f[i])
            if tempcol > maxcluster:
                maxcluster = tempcol
        n_comp += maxcluster

        for i in range(len(y_mp)):
            col = str(y_mp[i])
            tax = tmpf[i]
            lab = tax + "." + col.replace("C", "")

            secondpass[tax] = lab

        for k, v in secondpass.items():
            finalStrTree = finalStrTree.replace(k, v)
        reclustlist[str(clu)] = secondpass

    print("--- %s START to Second pass seconds ---" % (time.time() - start_time))
    # if res in descend

    # =============================================================================
    #     tree = dendropy.Tree.get(
    #         path=OUTTREE,
    #         schema="newick",
    #         label=None,
    #         taxon_namespace=None,
    #         collection_offset=None,
    #         tree_offset=None,
    #         rooting="default-unrooted",
    #         edge_length_type=float,
    #         suppress_edge_lengths=False,
    #         extract_comment_metadata=True,
    #         store_tree_weights=False,
    #         finish_node_fn=None,
    #         case_sensitive_taxon_labels=False,
    #         preserve_underscores=False,
    #         suppress_internal_node_taxa=True,
    #         suppress_leaf_node_taxa=False,
    #         terminating_semicolon_required=True,
    #         ignore_unrecognized_keyword_arguments=False,
    #         )
    #     tree.encode_bipartitions()
    #     print(flag)
    #     print()
    #     groups=[]
    #
    #
    #     mppassdict = {}
    #
    #     tmpflag = flag
    #     reclustCount = 0
    #     reclustlist = {}
    #     n_comp = clf.n_components
    #
    #
    #     for f in flag:
    #         print(f)
    #         reclustCount+=1
    #
    #         tmpf = f.copy()
    #         subclust = 0
    #         y_mp = [0]*len(tmpf)
    #         clu = f[0].split('|')[-1]
    #
    #         while(len(f)>0):
    #             for edge in tree.leaf_edges():
    #                 tmpleaf = edge.bipartition.leafset_as_newick_string(tree.taxon_namespace).split(', (')[0]
    #                 if tmpleaf[2:-1] == f[0]:
    #                     leaf = edge
    #             maxres=0
    #             mp=True
    #             parentedge = leaf
    #             print("parentedge", leaf)
    #
    #             while(mp):
    #                     for adj_edg in parentedge.get_adjacent_edges():
    #                         #print(adj_edg.bipartition.leafset_as_newick_string(tree.taxon_namespace).split(', (')[0])
    #                         descend = adj_edg.bipartition.leafset_as_newick_string(tree.taxon_namespace).split(', (')[0]
    #                         descend = descend.replace('(','')
    #                         descend = descend.replace(')','')
    #                         descend = descend.split(',')
    #                         for i in range(len(descend)):
    #                             descend[i] = descend[i].replace(' ','')
    #                         res =len([w for w in descend if w in f])
    #                         if res > maxres:
    #                             parentedge = adj_edg
    #                             maxres=res
    #                             print(res,descend)
    #                             temppartition = descend
    #                     check =  all(item in f for item in descend)
    #                     mp=check
    #                     print()
    #                     if check is True:
    #                         print("This edge is still monophyletic {} contains only elements of the list {}".format(descend, f))
    #                         partition = temppartition
    #                     elif  res >2:
    #                         for x in f:
    #                             ogindex = tmpf.index(x)
    #                             #print('siubclust',subclust,len(y_mp))
    #                             print(tmpf.index(x))
    #                             print(f.pop(f.index(x)))
    #                             y_mp[ogindex] = subclust
    #                         subclust=subclust+1
    #
    #                     else :
    #                         print("No, this edge went to far partition at the previous edge.")
    #                         for x in partition:
    #                             print(tmpf)
    #                             print('f:',f)
    #                             ogindex = tmpf.index(x)
    #                             #print('siubclust',subclust,len(y_mp))
    #                             print(tmpf.index(x))
    #                             print(f.pop(f.index(x)))
    #                             y_mp[ogindex] = subclust
    #                         subclust=subclust+1
    #
    #
    #         print(tmpf,y_mp)
    #         secondpass = {}
    #         maxcluster = 0
    #
    #         for i in range(len(y_mp)):
    #             tempcol = y_mp[i]
    #            # print(tempcol,f[i])
    #             if tempcol > maxcluster:
    #                 maxcluster = tempcol
    #         n_comp += maxcluster
    #
    #         for i in range(len(y_mp)):
    #             col = str(y_mp[i])
    #             tax = tmpf[i]
    #             lab = tax+'.'+col.replace('C','')
    #
    #             secondpass[tax] = lab
    #
    #         for k,v in secondpass.items():
    #             finalStrTree = finalStrTree.replace(k,v)
    #         reclustlist[str(clu)] = secondpass
    #
    #     print("--- %s START to Second pass seconds ---" % (time.time() - start_time))
    #             #if res in descend
    # =============================================================================

    ##second pass string rplace
    sectaxons = []
    seccatTax = []
    secrepDict = {}
    for i in range(n_comp):
        seccatTax.append([])
    print(catTax)

    sec_n = n_comp
    n_list = []
    count = 0
    for (
        k
    ) in (
        catTax
    ):  # ['sp|A6ZRW3|TOM70.YEAS7|7', 'sp|B0B815|OMCB.CHLT2|7', 'sp|P07213|TOM70.YEAST|7', 'sp|P0A3U8|BP26.BRUME|7', 'sp|P0A3U9|BP26.BRUSU|7', 'sp|P0CC04|OMCBD.CHLTR|7', 'sp|P0DJI2|OMCB.CHLTH|7', 'sp|P23603|OMCBE.CHLTH|7', 'sp|P23700|OMCB.CHLPN|7', 'sp|P26758|OMCBC.CHLTH|7', 'sp|P94664|OMCB.CHLCV|7', 'sp|Q253E5|OMCB.CHLFF|7']
        print("n:", k[0].split("|")[-1])
        n = k[0].split("|")[-1]  # .split('|')[-1]
        n_list.append(n)
        for v in k:  # for tom70 in arrary of tom70
            try:
                vtemp = reclustlist[n][v]
                print(vtemp)
                n2_temp = vtemp.split(".")[-1]
                n_temp = vtemp.split("|")[-1]
                if n2_temp == "0":
                    seccatTax[int(n)].append(vtemp)
                else:
                    seccatTax[int(sec_n - int(n2_temp))].append(vtemp)

            except KeyError:
                vtemp = v
                n_temp = v.split("|")[-1]
                seccatTax[int(n)].append(vtemp)

        sec_n -= 1
    print(seccatTax)

    #
    clusthash = {}
    for i in seccatTax:
        if i != []:
            for j in i:
                try:
                    key = j.split("|")[-1]
                    clusthash[key].append(j)
                except KeyError:
                    key = j.split("|")[-1]
                    clusthash[key] = []
                    clusthash[key].append(j)

    reptax = {}
    for i in seccatTax:
        if i != []:
            for j in i:
                reptax[".".join(j.split(".")[:-1])] = j

    for k, v in reptax.items():
        finalStrTree.replace(k, v)
    print(finalStrTree)

    #########   WRITE TREE
    EMtree_file = open(OUTTREE, "w+")
    n = EMtree_file.write(finalStrTree)
    EMtree_file.close()

    dictpred = {}
    for clustlist in seccatTax:
        if len(clustlist) > 0:
            # print(clustlist)
            key = str(clustlist[0].split("|")[-1])
            dictpred[key] = []
            for val in clustlist:
                key = str(val.split("|")[-1])

                try:
                    key = str(val.split("|")[-1])
                    dictpred[key].append("|".join(val.split("|")[:-1]))

                except KeyError:
                    dictpred[key] = []
                    dictpred[key].append("|".join(val.split("|")[:-1]))

    print((dictpred))
    print(time.time() - start_time)

    treeOR = dendropy.Tree.get(
        path=OUTTREE,
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
        taxons.append(str(taxon)[1:-1])
        if FAM == "control":
            catDict["|".join(str(taxon)[1:-1].split("|")[:-1])] = str(taxon)[1:-1]
        elif TREEID == "PF00528":
            catDict["|".join(str(taxon)[1:-1].split("|")[:-1])] = str(taxon)[1:-1]
        else:
            catDict["|".join(str(taxon)[1:-1].split("|")[:-1])] = str(taxon)[1:-1]

    # dict for tree viewer coloring table
    treeViewColor = {}

    # make bic score array a numpy arrary
    bic = np.array(bic)
    plotColor = {}

    i = 0
    for k, v in dictpred.items():
        for val in v:
            kog = str(k)
            k = float(k)
            if k == round(k, 0):  #########################
                k = int(k)
                treeViewColor[k] = BPcolors[i % len(BPcolors)]
            plotColor[kog] = BPcolors[i % len(BPcolors)]
        i += 1

    color_iter = itertools.cycle(plotColor.values())

    # Plot the BIC scores
    plt.figure(figsize=(7, 8))
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
    plt.text(xpos - 1 + start, bic.min() * 0.97 + 0.03 * bic.max(), "*", fontsize=14)
    spl.set_xlabel("Number of Gaussian components")
    spl.legend([b[0] for b in bars], cv_types)
    plt.ylabel("BIC score")

    # Plot the winner
    color_iter = itertools.cycle(treeViewColor.values())
    cluster_iter = itertools.cycle(treeViewColor.keys())

    splot = plt.subplot(2, 1, 2)
    Y_ = clf.predict(X)

    for i, (mean, color, clust) in enumerate(zip(clf.means_, color_iter, cluster_iter)):
        color = treeViewColor[i]
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, marker="o", color=color)
        plt.annotate(clust, mean[:2], color=color)

    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")

    plt.xlim(min(embedding[:, 0]) - 0.2, max(embedding[:, 0]) + 0.2)
    plt.ylim(min(embedding[:, 1]) - 0.2, max(embedding[:, 1]) + 0.2)
    plt.title(
        "Selected GMM: full model, " + str(clf.n_components) + " Gaussian components"
    )
    plt.subplots_adjust(hspace=0.35, bottom=0.02)
    plt.savefig(
        "../data/EMMplots/2022-01-11_"
        + str(n_H)
        + "_GMMsweep_"
        + FAM
        + TREEID
        + "_"
        + metric
        + ".svg",
        dpi=600,
        format="svg",
    )
    plt.show()

    OUTF = "../data/" + FAM + "colored/2022-01-11" + FAM + TREEID + "_colored_RF.csv"
    print(OUTF)
    import csv

    with open(OUTF, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["Label", "Cluster", "Color", "PassOne", "PassOneCol"])
        i = 0
        for k, v in dictpred.items():
            for val in v:
                k1 = int(float(k))
                if TREEID == "dehydrogenase":
                    val = val.replace("_", " ")

                print([val, k, BPcolors[i % len(BPcolors)], k1, treeViewColor[k1]])
                writer.writerow(
                    [val, k, BPcolors[i % len(BPcolors)], k1, treeViewColor[k1]]
                )
            i += 1

    OUTF = (
        "../data/" + FAM + "colored/2022-01-11" + FAM + TREEID + "_coloredPASSONE.csv"
    )
    import csv

    with open(OUTF, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["Label", "Cluster", "Color"])
        i = 0
        for k, v in dictpred.items():
            for val in v:
                k1 = int(float(k))
                writer.writerow([val, k1, treeViewColor[k1]])
            i += 1

    print(OUTF)

    print(time.time() - start_time)

    # OPEN reclustBPtree
    from io import StringIO

    # make tree color-able
    reclustBPtree = Phylogeny.from_tree(BPtree)
    reclustBPtree.rooted = True

    handle = StringIO(finalStrTree)
    reclustBPtree = Phylo.read(OUTTREE, "newick")
    reclustBPtree.ladderize()

    # color tree
    reclustBPtree.root.color = (128, 128, 128)

    mps = []
    k = 0
    for v in clusthash.values():
        print(v)
        if len(v) == 1:
            try:
                mrca = reclustBPtree.common_ancestor(v)
                indexcol = v[0].split("|")[-1]
                mrca.color = plotColor[indexcol]
            except:
                pass
                mrca = reclustBPtree.common_ancestor(".".join(v[0].split(".")[:-1]))
                indexcol = v[0].split("|")[-1]
                mrca.color = plotColor[indexcol]
        for i in range(len(v) - 1):
            try:
                mrca = reclustBPtree.common_ancestor(v[i], v[i + 1])
                indexcol = v[i].split("|")[-1]
                mrca.color = plotColor[indexcol]
            except ValueError:
                # Phylo.draw_ascii(reclustBPtree)

                mrca = reclustBPtree.common_ancestor(
                    ".".join(v[i].split(".")[:-1]), ".".join(v[i + 1].split(".")[:-1])
                )
                indexcol = v[i].split("|")[-1]
                mrca.color = plotColor[indexcol]
                pass
        k += 1

    with plt.rc_context({"lines.linewidth": 4}):
        matplotlib.rc("font", size=0.0)
        fig = plt.figure(figsize=(8, HEIGHT), dpi=200)
        axes = fig.add_subplot(1, 1, 1)
        plt.rc(
            "font", size=6
        )  # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
        plt.rc("axes", titlesize=14)  # fontsize of the axes title
        plt.rc("xtick", labelsize=18)  # fontsize of the tick labels
        plt.rc("ytick", labelsize=18)  # fontsize of the tick labels
        plt.rc("figure", titlesize=18)  # fontsize of the figure title
        Phylo.draw(reclustBPtree, axes=axes, do_show=False, branch_labels=None)
        plt.savefig(
            "../data/"
            + FAM
            + "colored/2022-01-11_"
            + str(n_H)
            + "_EMClust_"
            + FAM
            + TREEID
            + "_"
            + metric
            + "_coloredtree.svg",
            format="svg",
            dpi=600,
        )
        plt.show()


# =============================================================================
#             dn = dendro
#
#     # =============================================================================
#     #             ax.annotate(labels[ind], (x,y), va='top', ha='center')
#     #         plt.tight_layout()
#     #         plt.savefig('./tmp.png')
#     #         plt.close(fig)
#     # =============================================================================
#             # print(y_hc)
#             # print(f)
#
#             #######################GMM UMAP
#     # =============================================================================
#     #         N_2 = 3
#     #         reducer = umap.UMAP(n_components=N_2,
#     #                             n_neighbors=n_H,
#     #                             random_state=222,
#     #                             metric = 'precomputed')
#     #
#     #
#     #         embedding = reducer.fit_transform(pdm2)
#     #
#     #
#     #         labels = ['PC' + str(x) for x in range(0,N_2) ]
#     #
#     #         df = pd.DataFrame(taxons,columns = ['target'])
#     #
#     #         principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
#     #
#     #         finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
#     #
#     #         #######################CLUTSER
#     #
#     #         X = embedding[:, :]
#     #
#     #         start = 1
#     #
#     #         lowest_bic = np.infty
#     #         bic = []
#     #         n_components_range = range(1,int(0.75 * len(X)))
#     #
#     #         cv_types = ['spherical','diag','tied','full']#'tied',
#     #
#     #         for cv_type in cv_types:
#     #             for n_components in n_components_range:
#     #                 # Fit a Gaussian mixture with EM
#     #                 gmm = mixture.GaussianMixture(n_components=n_components,
#     #                                               random_state = 111,
#     #                                               covariance_type=cv_type,
#     #                                               reg_covar=10e-06)
#     #                 gmm.fit(X)
#     #                 bic.append(gmm.bic(X))
#     #                 if bic[-1] < lowest_bic:
#     #                     lowest_bic = bic[-1]
#     #                     best_gmm = gmm
#     #
#     #
#     #         bic = np.array(bic)
#     #         BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
#     #                      "#D55E00", "#CC79A7",
#     #                      '#f46d43',"#abd9e9"]
#     #
#     #         if (FAM == 'OMP')|(FAM == 'nextstrain')|(FAM == 'pfam')|(FAM == 'control')|(FAM == 'viral'):
#     #             BPcolors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
#     #                     '#ba3fdf', '#42d4f4', '#bfef45', '#fabed4',
#     #                     '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
#     #                     '#aaffc3', '#808000', '#ffd8b1', '#0072B2']
#     #         color_iter = itertools.cycle(BPcolors)
#     #         clf = best_gmm
#     #         bars = []
#     #
#     #         # Plot the BIC scores
#     #         plt.figure(figsize=(6,6))
#     #         spl = plt.subplot(2, 1, 1)
#     #         for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
#     #             xpos = np.array(n_components_range) + .2 * (i - 2)
#     #             bars.append(plt.bar(xpos, bic[i * len(n_components_range):
#     #                                           (i + 1) * len(n_components_range)],
#     #                                 width=.2, color=color))
#     #         plt.xticks(n_components_range)
#     #         plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
#     #         plt.title('BIC score per model')
#     #         xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
#     #             .2 * np.floor(bic.argmin() / len(n_components_range))
#     #         plt.text(xpos-1+start, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
#     #         spl.set_xlabel('Number of Gaussian components')
#     #         spl.legend([b[0] for b in bars], cv_types)
#     #         plt.ylabel("BIC score")
#     #
#     #         # Plot the winner
#     #         color_iter = itertools.cycle(BPcolors)
#     #
#     #         splot = plt.subplot(2, 1, 2)
#     #         Y_ = clf.predict(X)
#     #
#     #         for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
#     #                                                    color_iter)):
#     #             if not np.any(Y_ == i):
#     #                 continue
#     #
#     #             plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
#     #
#     #
#     #         plt.ylabel("UMAP Dimension 2")
#     #         plt.xlabel("UMAP Dimension 1")
#     #
#     #         plt.xlim(min(embedding[:,0])-0.2,max(embedding[:,0])+0.2)
#     #         plt.ylim(min(embedding[:,1])-0.2,max(embedding[:,1])+0.2)
#     #         plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
#     #         plt.subplots_adjust(hspace=.35, bottom=.02)
#     #         plt.savefig('../data/EMMplots/2021-10-29_PASS2'+f2.split('|')[-1]+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.svg',dpi=600,format='svg')
#     #         plt.show()
#     #
#     # =============================================================================
#
#             secondpass = {}
#             maxcluster = 0
#
#             #for i in range(len(f)):
#                 #print(y_hc[i],clf.predict(X)[i], f[i])
#
#
#             for i in range(len(y_hc)):
#                 tempcol = y_hc[i]
#                # print(tempcol,f[i])
#                 if tempcol > maxcluster:
#                     maxcluster = tempcol
#             n_comp += maxcluster
#
#             for i in range(len(y_hc)):
#                 col = str(y_hc[i])
#                 tax = f[i]
#                 lab = tax+'.'+col.replace('C','')
#
#                 secondpass[tax] = lab
#
#
#             for k,v in secondpass.items():
#                 finalStrTree = finalStrTree.replace(k,v)
#             reclustlist[str(clu)] = secondpass
#
#     print(finalStrTree)
#
#     i=0
#     for k in catTax:
#         print(k)
#         try:
#             print(reclustlist[str(i)]) #'sp|A6ZRW3|TOM70.YEAS7|7': 'sp|A6ZRW3|TOM70.YEAS7|7.2',
#                                        #'sp|P07213|TOM70.YEAST|7': 'sp|P07213|TOM70.YEAST|7.2',
#                                        #'sp|P0A3U8|BP26.BRUME|7': 'sp|P0A3U8|BP26.BRUME|7.0'
#
#         except:
#              pass
#         i+=1
#
#
#     sectaxons = []
#     seccatTax = []
#     secrepDict = {}
#     for i in range(n_comp):
#         seccatTax.append([])
#
#     for i in range(len(taxons)):
#         print(taxons[i])        #sp|Q9Z3D6|PMP12.CHLPN
#                                 #sp|Q9Z3D6|PMP12.CHLPN|4
#         print(taxonsLabels[i])
#
#     sec_n = n_comp
#     print('n_comp:',n_comp)
#     n_list = []
#     for k in catTax:  # ['sp|A6ZRW3|TOM70.YEAS7|7', 'sp|B0B815|OMCB.CHLT2|7', 'sp|P07213|TOM70.YEAST|7', 'sp|P0A3U8|BP26.BRUME|7', 'sp|P0A3U9|BP26.BRUSU|7', 'sp|P0CC04|OMCBD.CHLTR|7', 'sp|P0DJI2|OMCB.CHLTH|7', 'sp|P23603|OMCBE.CHLTH|7', 'sp|P23700|OMCB.CHLPN|7', 'sp|P26758|OMCBC.CHLTH|7', 'sp|P94664|OMCB.CHLCV|7', 'sp|Q253E5|OMCB.CHLFF|7']
#         n = k[0].split('|')[-1] #.split('|')[-1]
#         print(n)
#         n_list.append(n)
#         for v in k: #for tom70 in arrary of tom70
#             print(v)
#             try:
#                 vtemp = reclustlist[n][v]
#                 n2_temp = vtemp.split('.')[-1]
#                 n_temp = vtemp.split("|")[-1]
#                 print(vtemp)
#                 print(n_temp)
#                 if n2_temp == '0':
#                     seccatTax[int(n)].append(vtemp)
#                     # print(seccatTax[int(n)])
#                 else:
#                     print(int(sec_n-int(n2_temp))+1)
#                     seccatTax[int(sec_n-int(n2_temp))].append(vtemp)
#
#
#
#             except KeyError:
#                 vtemp = v
#                 n_temp = v.split("|")[-1]
#                 print(vtemp)
#                 print(n_temp)
#                 seccatTax[int(n)].append(vtemp)
#                 # print(seccatTax[int(n)])
#
#         sec_n -= 1
#
#     clusthash = {}
#     for i in seccatTax:
#         if i != []:
#             for j in i :
#                 try:
#                     key = j.split('|')[-1]
#                     clusthash[key].append(j)
#                 except KeyError:
#                     key = j.split('|')[-1]
#                     clusthash[key] = []
#                     clusthash[key].append(j)
#
#
#     reptax = {}
#     for i in seccatTax:
#         if i != []:
#             for j in i:
#                 reptax['.'.join(j.split('.')[:-1])] = j
#
#     for k,v in reptax.items():
#         finalStrTree.replace(k,v)
#     print(finalStrTree)
#
# #########   WRITE TREE
#     EMtree_file = open(OUTTREE,'w+')
#     n = EMtree_file.write(finalStrTree)
#     EMtree_file.close()
#
#
#
# #############################################################
#
#
#
#
#
#     ### Change for secondpass
#
#
# # =============================================================================
#     dictpred = {}
#     for clustlist in seccatTax:
#         if len(clustlist)>0:
#             print(clustlist)
#             key = str(clustlist[0].split('|')[-1])
#             dictpred[key] = []
#             for val in clustlist:
#                 key = str(val.split('|')[-1])
#
#                 try:
#                     key = str(val.split('|')[-1])
#                     dictpred[key].append('|'.join(val.split("|")[:-1]))
#
#                 except KeyError:
#                     dictpred[key] = []
#                     dictpred[key].append('|'.join(val.split("|")[:-1]))
#
#
# print(len(dictpred))
# print(time.time() - start_time)
#
#
#
#
#
#
#
#
# ###OPEN reclustBPtree
# from io import StringIO
#
# handle = StringIO(finalStrTree)
# reclustBPtree = Phylo.read(OUTTREE,"newick")
# reclustBPtree.ladderize()
#
# #make tree color-able
# reclustBPtree.rooted = True
# reclustBPtree = Phylogeny.from_tree(BPtree)
#
#
#
# bic = np.array(bic)
# BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
#              "#D55E00", "#CC79A7",
#              '#f46d43',"#abd9e9"]
#
# if (FAM == 'OMP')|(FAM == 'nextstrain')|(FAM == 'pfam')|(FAM == 'control')|(FAM == 'viral'):
#     BPcolors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
#             '#ba3fdf', '#42d4f4', '#bfef45', '#fabed4',
#             '#469990', '#dcbeff', '#9A6324', '#b3a400', '#800000',
#             '#82baff', '#808000', '#ffd8b1', '#0072B2',"#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
#              "#D55E00", "#CC79A7",
#              '#f46d43',"#abd9e9"]
#
#
#
#
#
#
#
#
#
# treeOR = dendropy.Tree.get(
#     path=OUTTREE,
#     schema="newick",
#     label=None,
#     taxon_namespace=None,
#     collection_offset=None,
#     tree_offset=None,
#     rooting="default-unrooted",
#     edge_length_type=float,
#     suppress_edge_lengths=False,
#     extract_comment_metadata=True,
#     store_tree_weights=False,
#     finish_node_fn=None,
#     case_sensitive_taxon_labels=False,
#     preserve_underscores=False,
#     suppress_internal_node_taxa=True,
#     suppress_leaf_node_taxa=False,
#     terminating_semicolon_required=True,
#     ignore_unrecognized_keyword_arguments=False,
#     )
# treeOR.ladderize()
#
# fileType = 'a'
# i = 0
# taxons = []
# catDict = {}
# for taxon in treeOR.taxon_namespace:
#     j = 0
#     taxons.append(str(taxon)[1:-1])
#     print('|'.join(str(taxon)[1:-1].split('|')[:-1]))
#     if FAM == 'control':
#         catDict['|'.join(str(taxon)[1:-1].split('|')[:-1])] = str(taxon)[1:-1]
#     elif TREEID=='PF00528':
#         catDict['|'.join(str(taxon)[1:-1].split('|')[:-1])] = str(taxon)[1:-1]
#     else:
#         catDict['|'.join(str(taxon)[1:-1].split('|')[:-1])] = str(taxon)[1:-1]
#
#
#
#
# #'0': ['A', 'P', 'T'], '1': ['C', 'F', 'W', 'Y'], '2': ['D', 'E'], \
# #'3': ['I', 'L', 'M', 'V'], '4': ['H', 'K'], '5': ['R', 'Q'],
# #'6': ['N', 'G', 'S']
#
# #'0': ['A', 'P', 'T'],
# dict0 = {}
# #'1': ['C', 'F', 'W', 'Y'],
# dict1 = {}
# #'2': ['D', 'E'],
# dict2 = {}
# #'3': ['I', 'L', 'M', 'V'],
# dict3 = {}
# #'4': ['H', 'K'],
# dict4 = {}
# #'5': ['R', 'Q'],
# dict5 = {}
# #'6': ['N', 'G', 'S']
# dict6 = {}
#
# gapsDict = {}
# seqCount = -1
# taxa = []
# alnDict = {}
#
# gapsDict['clust'] = []
# for i in range(Naln+1):
#     gapsDict[i]=[0]*(len(taxons))
#     dict0[str(i)+'APT']=[0]*(len(taxons))
#     dict1[str(i)+'CFWY']=[0]*(len(taxons))
#     dict2[str(i)+'DE']   =[0]*(len(taxons))
#     dict3[str(i)+'ILMV']=[0]*(len(taxons))
#     dict4[str(i)+'HK']=[0]*(len(taxons))
#     dict5[str(i)+'RQ']=[0]*(len(taxons))
#     dict6[str(i)+'NGS']=[0]*(len(taxons))
#
#
#
# if (ALNPATH[-5:] == "fasta")|(ALNPATH[-3:] == "txt"):
#
#     with open(ALNPATH,'r') as reader:
#         for l in reader:
#             l=l.strip()
# # =============================================================================
# #             print(l)
# # =============================================================================
#
#             if l[0] == '>':
#                 seqCount += 1
#                 print(catDict)
#                 print(l)
#
#                 if (FAM =='control') | (TREEID == 'sarscov2sgly') | (TREEID == 'PF00528'):
#                     try:
#                         taxon = catDict[l.strip()[1:].split(" ")[0].replace('_',REP)]
#                     except:
#                         taxon = catDict[l.strip()[1:].split(" ")[0].replace('YP',REP)]
#
#                 else:
#                     print(l.strip()[1:].split(" ")[0])
#                     taxon = catDict[l.strip()[1:].split(" ")[0]]
#                 gapsDict['clust'].append(taxon.split("|")[-1])
#                 taxa.append(taxon)
#
#                 counter = 0
#                 line = ''
#             else:
#                 line = l
#
#             gapCounter = 0
#
#             for char in line:
#                 if char == '-':
#                 #if (char == 'H') | (char == 'K'):
#                     gapCounter+=1
#                     gapsDict[counter][seqCount] = gapCounter
#                 elif (char == 'A') | (char == 'P') | (char == 'T'):
#                     gapCounter = 0
#                     dict0[str(counter)+'APT'][seqCount] = 1
#
#                 elif (char == 'C') | (char == 'F') | (char == 'W') | (char == 'Y'):
#                     gapCounter = 0
#                     dict1[str(counter)+'CFWY'][seqCount] = 1
#
#                 elif (char == 'D') | (char == 'E'):
#                     gapCounter = 0
#                     dict2[str(counter)+'DE'][seqCount] = 1
#
#                 elif (char == 'I') | (char == 'L') | (char == 'V') | (char == 'M'):
#                     gapCounter = 0
#                     dict3[str(counter)+'ILMV'][seqCount] = 1
#
#                 elif (char == 'H') | (char == 'K'):
#                     gapCounter = 0
#                     dict4[str(counter)+'HK'][seqCount] = 1
#
#                 elif (char == 'R') | (char == 'Q'):
#                     gapCounter = 0
#                     dict5[str(counter)+'RQ'][seqCount] = 1
#
#                 elif (char == 'N') | (char == 'G') | (char == 'S'):
#                     gapCounter = 0
#                     dict6[str(counter)+'NGS'][seqCount] = 1
#
#                 else:
#                     gapCounter = 0
#
#                 counter += 1
#
#
#
# d = {}
# d['clust'] = gapsDict['clust']
# l = np.asarray(d['clust'], dtype=np.float32,order='C')
# print(d['clust'])
# print(gapsDict['clust'])
#
#
# maxclust = len(dictpred)#int(max(l)) +1
#
# found = {}
# c = int(maxclust)
# for i in range(len(gapsDict['clust'])):
#     clusti = gapsDict['clust'][i]
#     if int(float(gapsDict['clust'][i])) == round(float(gapsDict['clust'][i]),1):
#         d['clust'][i] = int(float(d['clust'][i]))
#     else:
#         if clusti in found.keys():
#             print(clusti)
#         else:
#             found[clusti] = c
#             c+=1
# print()
# for i in range(len(d['clust'])):
#     clus = d['clust'][i]
#     # print(type(clus))
#     if isinstance(clus, str):
#         # print(clus)
#         d['clust'][i] = found[clus]
#
#
#
#
# for k in gapsDict.keys():
#     if k != 'clust':
#         if sum(gapsDict[k]) > 3 :
#             medianGap = statistics.median(gapsDict[k])
#             d[str(k)] = gapsDict[k]
#
#             for i in range(len(gapsDict[k])):
#                d[str(k)][i] = int(d[str(k)][i] > 0)#medianGap)
#             #d[str(k)] = gapsDict[k]
#             if len(gapsDict[k]) != len(taxa):
#                 # print('woah')
#                 print(k,len(gapsDict[k]))
#
# d0={}
# for k in dict0.keys():
#     if k != 'clust':
#         if sum(dict0[k]) > 3:
#             d0[str(k)] = dict0[k]
#
# d1={}
# for k in dict1.keys():
#     if k != 'clust':
#         if sum(dict1[k]) > 3:
#             d1[str(k)] = dict1[k]
# d2={}
# for k in dict2.keys():
#     if k != 'clust':
#         if sum(dict2[k]) > 3:
#             d2[str(k)] = dict2[k]
#
# d3={}
# for k in dict3.keys():
#     if k != 'clust':
#         if sum(dict3[k]) > 3:
#             d3[str(k)] = dict3[k]
#
# d4={}
# for k in dict4.keys():
#     if k != 'clust':
#         if sum(dict4[k]) > 3:
#             d4[str(k)] = dict4[k]
#
# d5={}
# for k in dict5.keys():
#     if k != 'clust':
#         if sum(dict5[k]) > 3:
#             d5[str(k)] = dict5[k]
#
# d6={}
# for k in dict6.keys():
#     if k != 'clust':
#         if sum(dict6[k]) > 3:
#             d6[str(k)] = dict6[k]
#
# df = pd.DataFrame(data=d,index=taxa)
# df0 = pd.DataFrame(data=d0,index=taxa)
# df1 = pd.DataFrame(data=d1,index=taxa)
# df2 = pd.DataFrame(data=d2,index=taxa)
# df3 = pd.DataFrame(data=d3,index=taxa)
# df4 = pd.DataFrame(data=d4,index=taxa)
# df5 = pd.DataFrame(data=d5,index=taxa)
# df6 = pd.DataFrame(data=d6,index=taxa)
# df = df.merge(df0,left_index=True,right_index=True)
# df = df.merge(df1,left_index=True,right_index=True)
# df = df.merge(df2,left_index=True,right_index=True)
# df = df.merge(df3,left_index=True,right_index=True)
# df = df.merge(df4,left_index=True,right_index=True)
# df = df.merge(df5,left_index=True,right_index=True)
# df = df.merge(df6,left_index=True,right_index=True)
#
#
#
# # Split the Groups from the dataset where y is category and x is data with species
# y = df.iloc[:,0]
# x = df.iloc[:,1:]
#
# # =============================================================================
# print(y)
# print(x.transpose())
# modelB = RandomForestClassifier(n_estimators=100,criterion='gini',random_state=111).fit(x, y)
#
# impdict = []
# index = []
# numfeatperclass2=[]
# numfeatperclass2
#
# impdf = pd.DataFrame(modelB.feature_importances_, columns = ['Relative Importance'],index=list(x.columns.values))
#
# imphead = impdf.sort_values('Relative Importance', ascending=False).head(10)
# print(imphead)
# fig = plt.figure(figsize=(8,6))
# ax = imphead.plot.barh()
# fig = ax.get_figure()
#
# fig.legend(loc='lower right')
# fig.savefig('../data/'+FAM+'colored/2021-10-29_'+FAM+TREEID+'_RFClassimportance.svg',dpi=600,format='svg')
# fig.show()
#
#
# featdict = {}
# for i in imphead.index:
#     featdict[i] = x[i]
#
#
# for i in imphead.index:
#     print(i,x[i])
#     featdict[i] = {}
#     for j in range(len(x[i])):
#         print(i,x[i][j],x.index[j])
#         featdict[i]['|'.join(x.index[j].split('|')[:-1])] = x[i][j]
#
#
# print(featdict)
#
#
#
#
#
# treeViewColor = {}
# plotColor = {}
#
# i=0
# for k,v in dictpred.items():
#     for val in v:
#         kog=str(k)
#         k=float(k)
#         if k == round(k,0):#########################
#             k=int(k)
#             treeViewColor[k] = BPcolors[i%len(BPcolors)]
#         plotColor[kog] = BPcolors[i%len(BPcolors)]
#     i+=1
#
#
#
# pprint.pprint(dictpred)
#
# OUTF="../data/"+FAM+"colored/2021-10-29"+FAM+TREEID+"_colored_RF.csv"
# print(OUTF)
# import csv
# with open(OUTF,'w') as f :
#     writer=csv.writer(f)
#     feats = list(featdict.keys())
#     writer.writerow(['Label','Cluster','Color','PassOne','PassOneCol',
#                      feats[-1],feats[-2], feats[-3], feats[-4], feats[-5],
#                      feats[-6],feats[-7], feats[-8], feats[-9], feats[-10]])
#     i = 0
#     for k,v in dictpred.items():
#         for val in v:
#             k1 = int(float(k))
#             if TREEID == 'dehydrogenase':
#                 val = val.replace('_',' ')
#             print(feats)
#             print(featdict)
#             print(k,val,BPcolors[i%len(BPcolors)],k1,treeViewColor[k1],
#                   feats[-1],featdict[feats[-1]][val],
#                   feats[-2],featdict[feats[-2]][val],
#                   feats[-3],featdict[feats[-3]][val],
#                   feats[-4],featdict[feats[-4]][val],
#                   feats[-5],featdict[feats[-5]][val],
#                   feats[-6],featdict[feats[-6]][val],
#                   feats[-7],featdict[feats[-7]][val],
#                   feats[-8],featdict[feats[-8]][val],
#                   feats[-9],featdict[feats[-9]][val],
#                   feats[-10],featdict[feats[-10]][val])
#
#             writer.writerow([val,k,BPcolors[i%len(BPcolors)],k1,treeViewColor[k1],
#                              featdict[feats[-1]][val],
#                              featdict[feats[-2]][val],
#                              featdict[feats[-3]][val],
#                              featdict[feats[-4]][val],
#                              featdict[feats[-5]][val],
#                              featdict[feats[-6]][val],
#                              featdict[feats[-7]][val],
#                              featdict[feats[-8]][val],
#                              featdict[feats[-9]][val],
#                              featdict[feats[-10]][val]])
#         i+=1
#
#
# # =============================================================================
# # print(clusthash)
# # =============================================================================
#
#
#
#
# # =============================================================================
# # print(numfeatperclass)
# # =============================================================================
# print(numfeatperclass2)
#
#
#
# OUTF="../data/"+FAM+"colored/2021-10-29"+FAM+TREEID+"_coloredPASSONE.csv"
# print(OUTF)
# import csv
# with open(OUTF,'w') as f :
#     writer=csv.writer(f)
#     feats = list(featdict.keys())
#     writer.writerow(['Label','Cluster','Color'])
#     i = 0
#     for k,v in dictpred.items():
#         for val in v:
#             k1 = int(float(k))
#             print(val,k1,treeViewColor[k1])
#
#             writer.writerow([val,k1,treeViewColor[k1]])
#         i+=1
#
#
# print(OUTF)
#
# tot = np.array(numfeatperclass2)>0.005
# total = tot.sum()
#
# print(len(dictpred))
# print(total)
# print(time.time() - start_time)
# tot = np.array(numfeatperclass2)>0.001
# total = tot.sum()
#
#
# print(total)
#
#
# color_iter = itertools.cycle(plotColor.values())
#
#
# # Plot the BIC scores
# plt.figure(figsize=(7,8))
# spl = plt.subplot(2, 1, 1)
# for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
#     xpos = np.array(n_components_range) + .2 * (i - 2)
#     bars.append(plt.bar(xpos, bic[i * len(n_components_range):
#                                   (i + 1) * len(n_components_range)],
#                         width=.2, color=color))
# plt.xticks(n_components_range)
# plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
# plt.title('BIC score per model')
# xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
#     .2 * np.floor(bic.argmin() / len(n_components_range))
# plt.text(xpos-1+start, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
# spl.set_xlabel('Number of Gaussian components')
# spl.legend([b[0] for b in bars], cv_types)
# plt.ylabel("BIC score")
#
# # Plot the winner
# color_iter = itertools.cycle(treeViewColor.values())
# cluster_iter = itertools.cycle(treeViewColor.keys())
#
# splot = plt.subplot(2, 1, 2)
# Y_ = clf.predict(X)
#
# for i, (mean, color,clust) in enumerate(zip(clf.means_, color_iter, cluster_iter)):
#     color = treeViewColor[i]
#     if not np.any(Y_ == i):
#         continue
#     plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22,marker='o', color=color)
#     plt.annotate(clust,mean,color=color)
#
# plt.ylabel("UMAP Dimension 2")
# plt.xlabel("UMAP Dimension 1")
#
# plt.xlim(min(embedding[:,0])-0.2,max(embedding[:,0])+0.2)
# plt.ylim(min(embedding[:,1])-0.2,max(embedding[:,1])+0.2)
# plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
# plt.subplots_adjust(hspace=.35, bottom=.02)
# plt.savefig('../data/EMMplots/2021-10-29_'+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.svg',dpi=600,format='svg')
# plt.show()
#
# plt.figure(figsize=(7,3.9))
# plt.scatter(
#     embedding[:, 0],
#     embedding[:, 1])
# plt.title('UMAP projection of the phylogenetic dataset')
# plt.savefig('../data/EMMplots/2021-10-29_'+str(n_H)+'_ogUMAPprog'+FAM+TREEID+'_'+metric+'.svg',dpi=600,format='svg')
#
# plt.show()
#
#
# #color tree
# reclustBPtree.root.color = (128, 128, 128)
#
#
# mps=[]
# k=0
# for v in clusthash.values():
#     print("V:",v)
#     if len(v)==1:
#         try:
#             mrca = reclustBPtree.common_ancestor(v)
#             indexcol = v[0].split("|")[-1]
#             mrca.color = plotColor[indexcol]
#         except:
#             pass
#             mrca = reclustBPtree.common_ancestor('.'.join(v[0].split('.')[:-1]))
#             indexcol = v[0].split("|")[-1]
#             mrca.color = plotColor[indexcol]
#     for i in range(len(v)-1):
#
#         print(v[i])
#         try:
#             mrca = reclustBPtree.common_ancestor(v[i], v[i+1])
#             indexcol = v[i].split("|")[-1]
#             mrca.color = plotColor[indexcol]
#         except ValueError:
#             # Phylo.draw_ascii(reclustBPtree)
#             print((reclustBPtree))
#             print(finalStrTree)
#             print(v[i], v[i+1])
#
#             mrca = reclustBPtree.common_ancestor('.'.join(v[i].split('.')[:-1]),
#                                                  '.'.join(v[i+1].split('.')[:-1]))
#             indexcol = v[i].split("|")[-1]
#             mrca.color = plotColor[indexcol]
#             pass
#     k+=1
#
#
# with plt.rc_context({'lines.linewidth': 4}):
#     matplotlib.rc('font', size=0.0)
#     fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
#     axes = fig.add_subplot(1, 1, 1)
#     plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
#     plt.rc('axes', titlesize=14)     # fontsize of the axes title
#     plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
#     plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
#     plt.rc('figure', titlesize=18)   # fontsize of the figure title
#     Phylo.draw(reclustBPtree,axes=axes,do_show=False,branch_labels=None)
#     plt.savefig("../data/"+FAM+"colored/2021-10-29_"+str(n_H)+"_EMClust_"+FAM+TREEID+"_"+metric+"_coloredtree.svg",format='svg',dpi=600)
#     plt.show()
#
# #make target vector
# labs = []
# for ii in finalDf['target']:
#     ii=ii.replace(' ','_')
#     labs.append(BPcolors[int(repDict[ii].split('|')[-1])%len(BPcolors)])#print(ii,repDict[ii].split('|')[-1])
# finalDf['labels'] = labs
#
# eteColors = ["","","","","","",""]
#
# eteTree = Tree(finalStrTree[4:].replace("'",""))
# print(eteTree)
#
# ts = TreeStyle()
# ts.show_leaf_name = False
#
# ts.mode = "c"
# ts.arc_start = -180 # 0 degrees = 3 o'clock
# ts.arc_span = 180
#
# ts.root_opening_factor = 1
#
# ts.min_leaf_separation = 1
# ts.optimal_scale_level = 10
# ts.layout_fn = layout
# ts.show_leaf_name = False
# print(dir(ts))
# cc=0
# for key in clusthash.values():
#     print(key)
#     print(plotColor)
#     nst = NodeStyle()
#     color = plotColor[key[0].split('|')[-1]] #BPcolors[cc%len(BPcolors)]
#     nst["vt_line_width"] = 12
#     nst["hz_line_width"] = 12
#     nst["bgcolor"] = color
#     print(nst)
#
#     for k in key:
#         print("KEY:",k.replace(' ','_'))
#
#         ke='|'.join(k.replace(' ','_').split('|')[:-1])
#         leaf = eteTree.search_nodes(name=ke)
#         leafb = eteTree.search_nodes(name="'"+ke+"'")
#         leafc = eteTree.search_nodes(name=k.replace(' ','_'))
#         if leaf != []:
#             leaf[0].set_style(nst)
#
#             for v in key:
#                 try:
#                     val='|'.join(v.replace(' ','_').split('|')[:-1])
#                     node = eteTree.get_common_ancestor(ke,val)
#                 except:
#                     v=v
#                     node = eteTree.get_common_ancestor(ke,v)
#                 node.set_style(nst)
#         elif leafb != []:
#             leafb[0].set_style(nst)
#             for v in key:
#                 #print(v)
#                 node = eteTree.get_common_ancestor("'"+ke+"'",v)
#                 node.set_style(nst)
#
#         elif leafc != []:
#             leafc[0].set_style(nst)
#             for v in key:
#                 try:
#                     node = eteTree.get_common_ancestor(k,v)
#
#                 except:
#                     break
#                     val='|'.join(v.replace(' ','_').split('|')[:-1])
#                     node = eteTree.get_common_ancestor(k,v)
#                 node.set_style(nst)
#         else:
#             pass
#     cc += 1
#
#
#
#
# eteTree.show(tree_style=ts)
# eteTree.render("../data/"+FAM+"colored/2021-10-29"+FAM+TREEID+"_colored.svg", tree_style=ts,dpi=600)
#
#
#
#
#
#
# =============================================================================
