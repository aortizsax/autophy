# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 15:15:32 2020

@author: aorti
"""
import dendropy 
import numpy as np
import scipy.spatial.distance as ssd
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt

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


pdmA = np.zeros((len(tree.taxon_namespace),len(tree.taxon_namespace)))
counter = 0
print(pdmA.shape)
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon))
    for taxon2 in tree.taxon_namespace:
        pdmA[i,j] = float(pdm.distance(taxon,taxon2))
        j+=1
    i+=1

pdist = ssd.pdist(pdmA)


df = pd.DataFrame(taxons,columns = ['target'])



pca = PCA(n_components=2)
pca.fit_transform(pdmA)
print(pca.explained_variance_ratio_)
principalComponents = pca.fit_transform(pdmA)
principalDf = pd.DataFrame(data = principalComponents,
              columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf,df[['target']]], axis = 1)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 

ax.set_xlabel('Principal Component 1: '+str(round(pca.explained_variance_ratio_[0], 2)), fontsize = 15)
ax.set_ylabel('Principal Component 2: '+str(round(pca.explained_variance_ratio_[1], 2)), fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

targets = ['OMP', 'PMP', 'PAL']
colors = ['r', 'g', 'b']
marks = ['o','v','*']

for target, color, mark in zip(targets,colors,marks):
    indicesToKeep = finalDf['target'].str.contains(target)
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50
               , marker = mark)
ax.legend(targets)
ax.grid()

