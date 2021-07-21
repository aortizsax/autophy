# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 12:03:16 2021

@author: aorti
"""

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

TREEPATH = "../../PhyloClust/data/2020-04-06_MaxParsimony_BSH_idFull.nwk" 
OUTTREE = "../data/2021_MP_OMP.nwk"

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
    taxons.append(str(taxon)[1:-1])
    for taxon2 in tree.taxon_namespace:
        pdmA[i,j] = float(pdm.distance(taxon,taxon2))
        j+=1
    i+=1

pdist = ssd.pdist(pdmA)

##setting data to final dataframe
df = pd.DataFrame(taxons,columns = ['target'])

##center and scale PDM AND fit
pdmScaled = preprocessing.scale(pdmA)

#make PCA object
pca = PCA(0.98,random_state = 99)
try:
    pca.fit(pdmA)
except:
    print("-=-=-=-=-=-=-=-=-ERROR: PCA failed. Rerun.-=-=-=-=-=-=-=-=-=")
    sys.exit()
   
principalComponents = pca.transform(pdmA)

#scree plot
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals =1)
labels = ['PC' + str(x) for x in range(1,len(per_var)+1) ]

plt.bar(x=range(1,len(per_var)+1), height = per_var, tick_label=labels)
plt.ylabel('Percntage of explained varience')
plt.xlabel('Principal Component')
plt.title('Scree Plot')
plt.show()


principalDf = pd.DataFrame(data = principalComponents,
              columns = labels)

finalDf = pd.concat([principalDf,df[['target']]], axis = 1)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 

ax.set_xlabel('Principal Component 1: '+str(round(pca.explained_variance_ratio_[0], 2)), fontsize = 15)
ax.set_ylabel('Principal Component 2: '+str(round(pca.explained_variance_ratio_[1], 2)), fontsize = 15)
ax.set_title('2 Component PCA', fontsize = 20)

targets = ['OMP', 'PMP', 'PAL']
colors = ['r', 'g', 'b']
marks = ['o','v','*']

for target, color, mark in zip(targets,colors,marks):
    indicesToKeep = finalDf['target'].str.contains(target)
    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
               , finalDf.loc[indicesToKeep, 'PC2']
               , c = color
               , s = 50
               , marker = mark)
ax.legend(targets)
ax.grid()
plt.show()





fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection = '3d')
ax.set_xlabel('Principal Component 1: '+str(round(pca.explained_variance_ratio_[0], 2)), fontsize = 15)
ax.set_ylabel('Principal Component 2: '+str(round(pca.explained_variance_ratio_[1], 2)), fontsize = 15)
ax.set_title('3 Component PCA', fontsize = 20)

targets = ['OMP', 'PMP', 'PAL']
colors = ['r', 'g', 'b']
marks = ['o','v','*']

for target, color, mark in zip(targets,colors,marks):
    indicesToKeep = finalDf['target'].str.contains(target)
    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
               , finalDf.loc[indicesToKeep, 'PC2']
               , finalDf.loc[indicesToKeep, 'PC3']
               , c = color
               , s = 60
               , marker = mark)
ax.legend(targets)
ax.grid()
plt.show()








#######################CLUTSER

X = principalComponents

#print(X)
lowest_bic = np.infty
bic = []
n_components_range = range(1, 22)
cv_types = ['spherical', 'tied', 'diag', 'full']
for cv_type in cv_types:
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components,
                                      random_state = 111,
                                      covariance_type=cv_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

bic = np.array(bic)
BPcolors = ['navy',#'cornflowerblue', 'turquoise', 'darkorange',
            'olive','teal','purple','gold',
                              'orange','salmon','red','cyan','tan']
color_iter = itertools.cycle(BPcolors)
clf = best_gmm
bars = []

# Plot the BIC scores
plt.figure(figsize=(8, 8))
spl = plt.subplot(2, 1, 1)
for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
    xpos = np.array(n_components_range) + .2 * (i - 2)
    bars.append(plt.bar(xpos, bic[i * len(n_components_range):
                                  (i + 1) * len(n_components_range)],
                        width=.2, color=color))
plt.xticks(n_components_range)
plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
plt.title('BIC score per model')
xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
    .2 * np.floor(bic.argmin() / len(n_components_range))
plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
spl.set_xlabel('Number of components')
spl.legend([b[0] for b in bars], cv_types)

# Plot the winner
splot = plt.subplot(2, 1, 2)
Y_ = clf.predict(X)
print(clf.covariances_.shape,clf.covariances_)
for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
                                           color_iter)):
    #v, w = linalg.eigh(cov)
    if not np.any(Y_ == i):
        continue
    plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

    # Plot an ellipse to show the Gaussian component
    #angle = np.arctan2(w[0][1], w[0][0])
    #angle = 180. * angle / np.pi  # convert to degrees
    #v = 2. * np.sqrt(2.) * np.sqrt(v)
    #ell = mpl.patches.Ellipse(mean, v[0]*40, v[1]*40, 180. + angle, color=color)
    #ell.set_clip_box(splot.bbox)
    #ell.set_alpha(.45)
    #splot.add_artist(ell)

plt.xticks(())
plt.yticks(())
plt.title('Selected GMM: full model, '+ str(clf.n_components) +' components')
plt.subplots_adjust(hspace=.35, bottom=.02)
plt.savefig('../data/2021-01-06_GMMsweepOMP.png')
plt.show()

print(clf.n_components)

# =============================================================================
# print(clf.predict(X))
# print(taxons)
# =============================================================================
taxonsLabels = []
repDict={}
catTax = []

for i in range(clf.n_components):
    catTax.append([])#{0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[]}
for i in range(0,len(taxons)):
    taxonsLabels.append(taxons[i]+"|"+str(clf.predict(X)[i]))
    print(clf.predict(X)[i],taxons[i],BPcolors[clf.predict(X)[i]%len(BPcolors)])
    #for coloring tree
    catTax[int(clf.predict(X)[i])].append(taxons[i]+"|"+str(clf.predict(X)[i]))
    #for tree label replacement
    repDict[taxons[i]] = taxons[i]+"|"+str(clf.predict(X)[i])

strTree = tree.as_string(schema="newick")
# =============================================================================
# print(strTree)
# =============================================================================
strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
strTree = re.sub('Inner[0-9][0-9]','', strTree)
strTree = re.sub('Inner[0-9]','', strTree)

for key in repDict.items():
    #print(key)
    strTree = strTree.replace(key[0],key[1])
    
# =============================================================================
# print(strTree)
# =============================================================================

EMtree_file = open(OUTTREE,'w+')
n = EMtree_file.write(strTree)
EMtree_file.close()





from Bio.Phylo.PhyloXML import Phylogeny

BPtree = Phylo.read(OUTTREE, "newick")
BPtree.ladderize()
#print(BPtree)
#Phylo.draw_ascii(BPtree)
#Phylo.draw(BPtree)
 
#make tree color-able
BPtree.rooted = True
BPtree = Phylogeny.from_tree(BPtree)
#BPtree = BPtree.as_phyloxml()

#color tree
BPtree.root.color = (128, 128, 128)
Phylo.draw(BPtree)

#BPcolors = list(color_iter)#["blue",'salmon','red','yellow','black','orange','purple','green','brown']
k=0
for v in catTax:
    print(k,v,BPcolors[k%len(BPcolors)])
    if len(v)==1:
        mrca = BPtree.common_ancestor(v)#[i], v[i+1])
        mrca.color = BPcolors[k%len(BPcolors)]
    for i in range(len(v)-1):
        mp = BPtree.is_monophyletic(v)
        mrca = BPtree.common_ancestor(v[i], v[i+1])
# =============================================================================
#         for clade in mrca:
#             print(clade)
#             clade.color = BPcolors[k%len(BPcolors)]
# =============================================================================
        #mrca = xmlBPtree.find_clade(v[i][1:-1])
       # print(mp, mrca,v[i],v[i+1]) 
        mrca.color = BPcolors[k%len(BPcolors)]
    k+=1     
    

with plt.rc_context({'lines.linewidth': 3}):
    matplotlib.rc('font', size=10)
    fig = plt.figure(figsize=(8,14), dpi=200)
    axes = fig.add_subplot(1, 1, 1)
    plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
    plt.rc('axes', titlesize=14)     # fontsize of the axes title
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    plt.rc('figure', titlesize=18)   # fontsize of the figure title
    Phylo.draw(BPtree,axes=axes,do_show=False,branch_labels=None)
    plt.savefig("../data/2021-01-06_EMClust_OMP_coloredtree.png")
    plt.show()#print(xmlBPtree)
