# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:05:17 2021

@author: aorti
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:05:52 2021

@author: aorti
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 16:13:27 2021

@author: aorti
"""
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

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=25)
        faces.add_face_to_node(N, node, 0, position="branch-right")



# =============================================================================
# TREEPATH = "../../../../ChromeDownloads/nextstrain_dengue_denv1_timetree.nwk" 
# OUTTREE = "../data/2021_flu.nwk"
# FAM = 'FLU'
# =============================================================================
# =============================================================================
# TREEPATH = "../../../Users/Adrian/Downloads/nextstrain_ncov_global_timetree.nwk" 
# OUTTREE = "../data/2021_rona.nwk"
# FAM = 'RONA'
# HEIGHT = 50
# =============================================================================

TREETYPE = 'ML'
TREEPATH = "../data/SPIKEtree/2021-01-24_SPIKE09_btsp50_"+TREETYPE+".nwk"
OUTTREE = "../data/SPIKEtree/clustered/2021-01-31_SPIKE09_"+TREETYPE+".nwk"
FAM = 'SPIKE'
HEIGHT = 22

# =============================================================================
# TREEPATH = "../data/SPIKEtree/2021-01-24_SPIKE09_btsp50_ML.nwk"
# OUTTREE = "../data/SPIKEtree/clustered/2021-02-02_SPIKE09_ML.nwk"
# FAM = 'SPIKE'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../../Taxon_Enzyme_Phylogeny/BSH_data/2020-03-05_MaxParsimony_BSH_id.nwk" 
# =============================================================================

TREEPATH = "../data/viraltrees/sarscov2sgly_muscle_x_cah.tree"
OUTTREE = "../data/sarscov2sgly_muscle_cah.new"
FAM = 'viral'
HEIGHT = 22

# =============================================================================
# TREEPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle_relaxed.tree"
# OUTTREE = "../data/2021-02-02_BSH_filterd_BEAST.nwk"
# FAM = 'BSH'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/OMPaln/2021-03-21_OMPuniprot_muscle.tree"
# OUTTREE = "../data/2021-03-21_OMPuniprot_muscle.tree"
# FAM = 'OMP'
# HEIGHT = 22
# =============================================================================


#boot for loop 


# =============================================================================
# #adjust 
# treeStr = open(TREEPATH,"r")
# treeStrLab = open(TREEPATH[:-4]+"LAB.nwk",'w')
# 
# print(treeStr)
# for line in treeStr:
#     treeStrLab.write(line.replace('_','.').replace(' ','.').replace('|org|','|').replace('|org-type|','|'))
#     print(line.replace('_','.').replace(' ','.'))
#     print()
# 
# treeStr.close()
# treeStrLab.close()
# =============================================================================


treeOR = dendropy.Tree.get(
    path=TREEPATH,
    schema="nexus",
    extract_comment_metadata=True)

ds = dendropy.DataSet.get_from_path(TREEPATH,
                                    "nexus",
                                    extract_comment_metadata=True)

# =============================================================================
# print(ds.as_ascii_plot())
# =============================================================================
# =============================================================================
# treeOR.ladderize()
# #treeOR.reroot_at_midpoint()
# #treeOR.reroot_at_midpoint(update_bipartitions=False)
# pdm = treeOR.phylogenetic_distance_matrix()
# print(pdm.as_data_table())
# print(treeOR.as_ascii_plot())
# =============================================================================

for a in ds.annotations:
    print("Data Set '%s': %s" % (treeOR.label, a))

for taxon_namespace in ds.taxon_namespaces:
    for a in taxon_namespace.annotations:
        print("Taxon Set '%s': %s" % (taxon_namespace.label, a))
    for taxon in taxon_namespace:
        for a in taxon.annotations:
            print("Taxon '%s': %s" % (taxon.label, a))
for tree_list in ds.tree_lists:
    for a in tree_list.annotations:
        print("Tree List '%s': %s" % (tree_list.label, a))
    for tree in tree_list:
        for a in tree.annotations:
            print("Tree '%s': %s" % (tree.label, a))
    
    
# =============================================================================
# for node in treeOR.postorder_node_iter():
#     print(dir(node))
#     print(node.label)
#     
# for edge in treeOR.postorder_edge_iter():
#     #print(dir(edge))
#     print(edge.annotations)
# =============================================================================
print(treeOR.as_string("newick"))


tree = dendropy.Tree.get(
    path=TREEPATH,
    schema="nexus",
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
tree.ladderize()
pdm = tree.phylogenetic_distance_matrix()
#print(pdm.as_data_table())


pdmA = np.zeros((len(tree.taxon_namespace),len(tree.taxon_namespace)+1))
counter = 0
print(pdmA.shape)
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon)[1:-1])
    print(taxon)
    for taxon2 in tree.taxon_namespace:
        edgecount = float(pdm.path_edge_count(taxon,taxon2))
        pdmA[i,j] =  float(pdm.distance(taxon,taxon2))
                
        j+=1
    i+=1
bootcount = 0
#i += 1
for node in tree.postorder_node_iter():#preorder_node_iter():
    if node.label != None:
        if float(node.label) > 0.8:
            #print(str(node.label))
            bootcount += 1#np.sqrt(float(node.label))
    if node.taxon != None:
        print(bootcount,str(node.taxon))
        index = taxons.index(str(node.taxon)[1:-1])
        #print(index,i)
        pdmA[index,i] = bootcount
i += 1
# =============================================================================
# for node in treeOR.preorder_node_iter():#preorder_node_iter():    if node.label != None:
#     if node.label != None:
#         if float(node.label) > 0.75:
#             #print(str(node.label))
#             bootcount += 1
#     if node.taxon != None:
#         print(bootcount,str(node.taxon).split(' ')[0])
#         index = taxons.index(str(node.taxon)[1:-1].split(' ')[0])
#         print(index,i)
#         pdmA[index,i] = bootcount
# =============================================================================

print(pdmA)

N=3
if FAM == 'PRACT':
    N=5
#pdist = ssd.pdist(pdmA)
#print(taxons)



metric = 'manhattan' #'euclidean'
metrics = ['manhattan','minkowski','canberra','braycurtis',#'haversine',
           'mahalanobis','wminkowski','cosine']
n_Hs = [int(len(taxons)*3/4)]#,100,125,150,200]# [20,75,100,222]#[20,100,150,180,220,1000,2000,3000]



#%matplotlib inline

#sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})
from Bio.Align import substitution_matrices
print(substitution_matrices.load())


print(substitution_matrices.load("JONES"))
matrix = substitution_matrices.load("JONES")

print(matrix.alphabet)
print(pdmA)
blosum62 = np.zeros((len(matrix.alphabet),len(matrix.alphabet)))
for i in range(0,len(matrix.alphabet)):
    print(matrix.alphabet[i])
    ia=matrix.alphabet[i]
    for j in range(0,len(matrix.alphabet)):
       # print(matrix.alphabet[j])
        ja=matrix.alphabet[j]
      #  print(matrix[ia][ja])
        blosum62[i][j] = matrix[i][j]
        
print(blosum62)
penguins = blosum62# pd.read_csv("https://github.com/allisonhorst/palmerpenguins/raw/5b5891f01b52ae26ad8cb9755ec93672f49328a8/data/penguins_size.csv")

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

linked = linkage(penguins, 'single')

labelList = matrix.alphabet# range(0, len())

plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()

import scipy.cluster.hierarchy as shc

plt.figure(figsize=(10, 7))
plt.title("Hierarchical Clustering Amino Acid Dendogram")
dend = shc.dendrogram(shc.linkage(penguins, method='weighted'),labels=labelList)
plt.show()


# =============================================================================
# PCoA(penguins)
# 
# =============================================================================

for n_H in n_Hs:
    print("metric:",n_H)
    reducer = umap.UMAP(n_components=N,    
                        n_neighbors=n_H,
                        min_dist=0.0,
                        random_state=42,
                        metric = metric)
    
    
    penguin_data = penguins
    scaled_penguin_data = StandardScaler().fit_transform(penguin_data)
    
    embedding = reducer.fit_transform(scaled_penguin_data)
    #print(embedding.shape)
    
    plt.scatter(
        embedding[:, 0],
        embedding[:, 1])#,
        #c=[sns.color_palette()[x] for x in penguins.species_short.map({"Adelie":0, "Chinstrap":1, "Gentoo":2})])
    #plt.gca().set_aspect('equal', 'box')
    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")
# =============================================================================
#     plt.grid()
# =============================================================================
    for i, txt in enumerate(list(matrix.alphabet)):
        print(txt, embedding[i,0],embedding[i,1])
        plt.annotate(txt, (float(embedding[i,0]),float(embedding[i,1])),size=15)
    plt.title('UMAP projection of the PDM', fontsize=24)
    plt.show()
    
    labels = ['PC' + str(x) for x in range(1,N+1) ]
    
    df = pd.DataFrame(list(matrix.alphabet),columns = ['target'])
    
    principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
    
    finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
       
    
    
    

    
    
    
    #######################CLUTSER
    
    X = embedding[:, :]
    
    print(pdmA,matrix.alphabet)
    lowest_bic = np.infty
    bic = []
    N_start = 4
    n_components_range = range(N_start,11)
    if FAM == 'PRACT':
        n_components_range = range(1,6)
    cv_types = ['diag', 'tied', 'full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          random_state = 111,
                                          covariance_type=cv_type,
                                          reg_covar=30e-06)
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
# =============================================================================
#     best_gmm = mixture.BayesianGaussianMixture(n_components=5,
#                                            random_state = 111,
#                                            covariance_type="full"
# # =============================================================================
# #                                            ,
# #                                            reg_covar=10e-06
# # =============================================================================
#                                            ).fit(X)
# =============================================================================
    bic = np.array(bic)
    BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7",
                 '#f46d43',"#abd9e9"]
                #'#a9a9a9',                '#ffffff']#, '#000000'
# =============================================================================
#                 ]
# =============================================================================
# =============================================================================
#                 'salmon','red','pink', 'orange',#'steelblue',#'cornflowerblue', 'turquoise', 'darkorange',
#                 'olive','teal','purple','gold',
#                 'cyan','tan','green',
#                 'silver','maroon',
#                 '#7B68EE','#87CEEB','#DDA0DD','brown','#D2691E','lime']#,'green','grey']
# =============================================================================
    if FAM == 'OMP':
        BPcolors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                '#ba3fdf', '#42d4f4', '#bfef45', '#fabed4', 
                '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', 
                '#aaffc3', '#808000', '#ffd8b1', '#0072B2']
    color_iter = itertools.cycle(BPcolors)
    clf = best_gmm
    bars = []
    
    # Plot the BIC scores
    plt.figure(figsize=(6,6))
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
    plt.text(xpos+N_start-1, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
    spl.set_xlabel('Number of Gaussian components')
    spl.legend([b[0] for b in bars], cv_types)
    plt.ylabel("BIC score")
    
    # Plot the winner
    color_iter = itertools.cycle(BPcolors)

    splot = plt.subplot(2, 1, 2)
    Y_ = clf.predict(X)
    print(clf.covariances_.shape,clf.covariances_,clf.means_)
    for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
                                               color_iter)):
# =============================================================================
#         v, w = linalg.eigh(cov)
# =============================================================================
        if not np.any(Y_ == i):
            continue
        
        print(i, mean, color)
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
    
# =============================================================================
#         # Plot an ellipse to show the Gaussian component
#         print(w[0][1], w[0][0])
#         angle = np.arctan2(-w[0][1]-0.01, -w[0][0]-0.01)
#         angle = 180. * angle / np.pi  # convert to degrees
#         v = 2. * np.sqrt(2.) * np.sqrt(v)
#         print("angle",angle)
#         ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. - angle, color=color)
#         ell.height=0.51
#         ell.width=0.1
#         if i == 2:
#             ell.height=0.3
#             ell.width=0.3
#         print(ell)
#         ell.set_clip_box(splot.bbox)
#         ell.set_alpha(.5)
#         splot.add_artist(ell)
# =============================================================================


    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")
# =============================================================================
#     plt.grid()
# =============================================================================
    plt.xlim(min(embedding[:,0])-0.2,max(embedding[:,0])+0.2)
    plt.ylim(min(embedding[:,1])-0.2,max(embedding[:,1])+0.2)
    for i, txt in enumerate(list(matrix.alphabet)):
        print(txt, embedding[i,0],embedding[i,1])
        plt.annotate(txt, (float(embedding[i,0]),float(embedding[i,1])),size=15)
    #plt.xticks(())
    #plt.yticks(())
    plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
    plt.subplots_adjust(hspace=.35, bottom=.02)
    plt.savefig('../data/EMMplots/2021-01-31_'+str(n_H)+'_GMMsweep_'+FAM+'_'+metric+'.png')
    plt.show()
    
    print(X)
    print(clf.predict(X))
    dictpred = {'0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[],
               '9':[],'10':[],'11':[],'12':[],'13':[],'14':[],'0':[],'0':[],'0':[],}
    for i in range(len(matrix.alphabet)):
        print(matrix.alphabet[i],clf.predict(X)[i])
        dictpred[str(clf.predict(X)[i])].append(matrix.alphabet[i]) 
    
    print(dictpred)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# =============================================================================
# import numpy as np
# import itertools
# 
# from scipy import linalg
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# 
# from sklearn import mixture
# 
# print(__doc__)
# 
# # Number of samples per component
# n_samples = 500
# 
# # Generate random sample, two components
# np.random.seed(0)
# C = np.array([[0., -0.1], [1.7, .4]])
# X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
#           .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]
# 
# lowest_bic = np.infty
# bic = []
# n_components_range = range(1, 7)
# cv_types = ['spherical', 'tied', 'diag', 'full']
# for cv_type in cv_types:
#     for n_components in n_components_range:
#         # Fit a Gaussian mixture with EM
#         gmm = mixture.GaussianMixture(n_components=n_components,
#                                       covariance_type=cv_type)
#         gmm.fit(X)
#         bic.append(gmm.bic(X))
#         if bic[-1] < lowest_bic:
#             lowest_bic = bic[-1]
#             best_gmm = gmm
# 
# bic = np.array(bic)
# color_iter = itertools.cycle(['navy', 'turquoise', 'cornflowerblue',
#                               'darkorange'])
# clf = best_gmm
# bars = []
# 
# # Plot the BIC scores
# plt.figure(figsize=(8, 6))
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
# plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
# spl.set_xlabel('Number of components')
# spl.legend([b[0] for b in bars], cv_types)
# 
# # Plot the winner
# splot = plt.subplot(2, 1, 2)
# Y_ = clf.predict(X)
# for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
#                                            color_iter)):
# # =============================================================================
# #     v, w = linalg.eigh(cov)
# # =============================================================================
#     if not np.any(Y_ == i):
#         continue
#     plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
# # =============================================================================
# # 
# #     # Plot an ellipse to show the Gaussian component
# #     print(w[0][1], w[0][0])
# #     angle = np.arctan2(w[0][1], w[0][0])
# #     angle = 180. * angle / np.pi  # convert to degrees
# #     v = 2. * np.sqrt(2.) * np.sqrt(v)
# #     print("angle",angle)
# #     ell = mpl.patches.Ellipse(mean, v[0], v[1], 180+angle, color=color)
# #     print(ell)
# #  
# #     ell.set_clip_box(splot.bbox)
# #     ell.set_alpha(.5)
# #     splot.add_artist(ell)
# # =============================================================================
# 
#     plt.xticks(())
#     plt.yticks(())
#     plt.title(f'Selected GMM: {best_gmm.covariance_type} model, '
#               f'{best_gmm.n_components} components')
#     plt.subplots_adjust(hspace=.35, bottom=.02)
#     plt.show()
# =============================================================================
# =============================================================================
#     print(plt.xlim(),plt.ylim())
#     print(clf.n_components)
#     
#     print(clf.predict(X))
#     print(clf.predict_proba(X))
#     for i in range(len(list(matrix.alphabet))):
#         #print(i)
#         print(matrix.alphabet[i],clf.predict(X)[i])
#         #print(np.around(clf.predict_proba(X)[i],decimals=10))
#         #print(clf.predict_proba(X)[i].dtype)
#     #print(taxons)
#     taxonsLabels = []
#     repDict={}
#     catTax = []
#     
#     for i in range(clf.n_components):
#         catTax.append([])#{0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[]}
#     for i in range(0,len(taxons)):
#         taxonsLabels.append(taxons[i]+"|"+str(clf.predict(X)[i]))
#         #print(clf.predict(X)[i],taxons[i],BPcolors[clf.predict(X)[i]%len(BPcolors)])
#         #for coloring tree
#         catTax[int(clf.predict(X)[i])].append(taxons[i]+"|"+str(clf.predict(X)[i]))
#         #for tree label replacement
#         repDict[taxons[i]] = taxons[i]+"|"+str(clf.predict(X)[i])
#     
#     strTree = tree.as_string(schema="newick")
#     # =============================================================================
#     print(strTree)
#     # =============================================================================
#     strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
#     strTree = re.sub('Inner[0-9][0-9]','', strTree)
#     strTree = re.sub('Inner[0-9]','', strTree)
#     
#     for key in repDict.items():
#         #print(key)
#         strTree = strTree.replace(key[0],key[1])
#         
#     # =============================================================================
#     print(strTree)
#     # =============================================================================
#     
#     import re
#     
#     #print(re.sub('\)\d:', '\).:',  strTree))
#     finalStrTree = strTree
#     for iter in re.finditer("\D(?<=:)",strTree):
#         a = iter.span()
#         b = a[1]
#         a = a[0]
#         if strTree[a-7] ==')':
#             print(strTree[a-7:b])
#             finalStrTree = finalStrTree.replace(strTree[a-7:b],'):')
#     print(finalStrTree)
#     EMtree_file = open(OUTTREE,'w+')
#     n = EMtree_file.write(finalStrTree)
#     EMtree_file.close()
#     
# =============================================================================
    
    
    
    
# =============================================================================
#     from Bio.Phylo.PhyloXML import Phylogeny
#     
#     BPtree = Phylo.read(OUTTREE, "newick")
#     BPtree.ladderize()
#     #print(BPtree)
#     #Phylo.draw_ascii(BPtree)
#     #Phylo.draw(BPtree)
#      
#     #make tree color-able
#     BPtree.rooted = True
#     BPtree = Phylogeny.from_tree(BPtree)
#     #BPtree = BPtree.as_phyloxml()
#     
#     #color tree
#     BPtree.root.color = (128, 128, 128)
#     Phylo.draw(BPtree)
#     
#     #BPcolors = list(color_iter)#["blue",'salmon','red','yellow','black','orange','purple','green','brown']
#     k=0
#     for v in catTax:
#         #print(k,v,BPcolors[k%len(BPcolors)])
#         if len(v)==1:
#             mrca = BPtree.common_ancestor(v)#[i], v[i+1])
#             mrca.color = BPcolors[k%len(BPcolors)]
#         for i in range(len(v)-1):
#             mp = BPtree.is_monophyletic(v)
#             try:
#                 mrca = BPtree.common_ancestor(v[i], v[i+1])
#                 mrca.color = BPcolors[k%len(BPcolors)]
#             except:
#                 pass
#         k+=1     
#         
#     
#     with plt.rc_context({'lines.linewidth': 4}):
#         matplotlib.rc('font', size=10)
#         fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
#         axes = fig.add_subplot(1, 1, 1)
#         plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
#         plt.rc('axes', titlesize=14)     # fontsize of the axes title
#         plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('figure', titlesize=18)   # fontsize of the figure title
#         Phylo.draw(BPtree,axes=axes,do_show=False,branch_labels=None)
#         plt.savefig("../data/"+FAM+"colored/2021-01-31_"+str(n_H)+"_EMClust_"+FAM+"_"+metric+"_coloredtree.png")
#         plt.show()
# 
# # =============================================================================
# #     #make target vector
# #     labs = []
# #     for ii in finalDf['target']:
# #         
# #         labs.append(BPcolors[int(repDict[ii].split('|')[-1])%len(BPcolors)])#print(ii,repDict[ii].split('|')[-1])
# #     finalDf['labels'] = labs
# #     print(finalDf)
# #     
# # 
# # 
# #     fig = plt.figure(figsize = (5,5))
# #     ax = fig.add_subplot(1,1,1) 
# #     
# #     ax.set_xlabel('Dimension 1 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
# #     ax.set_ylabel('Dimension 2 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
# #     ax.set_title('UMAP of selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
# # #    ax.set_title('2 Component UMAP', fontsize = 20)
# #     
# #     targets = [FAM]
# #     colors = [1,25,3]
# #     marks = ['o','v','*']
# # 
# #     ax.legend(targets)
# #     ax.grid()
# #     #plt.scatter(finalDf["PC1"],finalDf['PC2'],c=list(finalDf['labels']))
# #     
# #     
# #     for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
# #                                                color_iter)):
# # # =============================================================================
# # #         v, w = linalg.eigh(cov)
# # # =============================================================================
# #         if not np.any(Y_ == i):
# #             continue
# #         plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
# #     
# #     
# #     
# #     # Plot an ellipse to show the Gaussian component
# # # =============================================================================
# # #         angle = np.arctan2(w[0][1], w[0][0])
# # #         angle = 180. * angle / np.pi  # convert to degrees
# # #         v = 2. * np.sqrt(2.) * np.sqrt(v)
# # #         ell = mpl.patches.Ellipse(mean, v[0]*40, v[1]*40, 180. + angle, color=color)
# # #         ell.set_clip_box(splot.bbox)
# # #         ell.set_alpha(.45)
# # #         splot.add_artist(ell)
# # # =============================================================================
# #     
# #     
# #     
# #     
# #     
# #     
# #     
# #     
# #     
# #     
# #     plt.show()
# # 
# # #############################################################
# # 
# #     fig = plt.figure(figsize = (8,8))
# #     ax = plt.axes(projection = '3d')
# #     ax.set_xlabel('Dimension 1: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
# #     ax.set_ylabel('Dimension 2: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
# #     ax.set_title('3 Dimension UMAP', fontsize = 20)
# #     
# #     targets = [FAM]
# #     colors = ['r', 'g', 'b']
# #     marks = ['o','v','*']
# #     
# #     plt.scatter(finalDf['PC1']
# #                    , finalDf['PC2']
# #                    , finalDf['PC3']#,s=20
# #                    , color = list(finalDf['labels'])
# #                    , marker = 'D')
# #     ax.legend(targets)
# #     ax.grid()
# #     plt.show()
# #     
# # 
# # 
# # 
# #     eteColors = ["","","","","","",""]
# # 
# #     eteTree = Tree(strTree[4:].replace("'",""))
# #    # eteTree.convert_to_ultrametric()
# #     print(eteTree)
# #     
# #     ts = TreeStyle()
# #     ts.show_leaf_name = True
# #     #ts.mode = "c"
# #     ts.root_opening_factor = 1
# # 
# #     ts.min_leaf_separation = 1
# #     ts.optimal_scale_level = 10
# #     ts.layout_fn = layout
# #     ts.show_leaf_name = False
# #     print(dir(ts))
# # # =============================================================================
# # #     ts.arc_start = -180 # 0 degrees = 3 o'clock
# # #     ts.arc_span = 180
# # # =============================================================================
# #     cc=0
# #     for key in catTax:
# #         print(key)
# #         nst = NodeStyle()
# #         color = BPcolors[cc%len(BPcolors)]
# # # =============================================================================
# # #         nst["vt_line_color"] = color
# # #         nst["hz_line_color"] = color
# # # =============================================================================
# #         nst["vt_line_width"] = 12
# #         nst["hz_line_width"] = 12
# #         nst["bgcolor"] = color
# #         print(nst)
# #         
# #         for k in key:
# #             print("KEY:",k.replace(' ','_'))
# #             
# #             ke='|'.join(k.replace(' ','_').split('|')[:-1])
# #             leaf = eteTree.search_nodes(name=ke)
# #             leafb = eteTree.search_nodes(name="'"+ke+"'")
# #             leafc = eteTree.search_nodes(name=k)
# #             #print(leaf)
# #             if leaf != []:
# #                 #nst["fsize"] = 14
# #                 leaf[0].set_style(nst)
# # 
# #                 for v in key:
# #                     #print(v)
# #                     #v='|'.join(v.replace(' ','_').split('|')[:-1])
# #                     try:
# #                         val='|'.join(v.replace(' ','_').split('|')[:-1])
# #                         node = eteTree.get_common_ancestor(ke,val)
# #                     except:
# #                         v=v
# #                         node = eteTree.get_common_ancestor(ke,v)
# #                     node.set_style(nst)
# #             elif leafb != []:
# #                 leafb[0].set_style(nst)
# #                 for v in key:
# #                     #print(v)
# #                     node = eteTree.get_common_ancestor("'"+ke+"'",v)
# #                     node.set_style(nst)
# #                     
# #             else: 
# #                 leafc[0].set_style(nst)
# #                 for v in key:
# #                     #print(v)
# #                     try:
# #                         node = eteTree.get_common_ancestor(k,v)
# #                         
# #                     except:
# #                         val='|'.join(v.replace(' ','_').split('|')[:-1])
# #                         node = eteTree.get_common_ancestor(k,val)
# #                     node.set_style(nst)
# #         
# #         cc += 1
# #         
# #         
# #         
# #     
# #     eteTree.show(tree_style=ts)
# #     eteTree.render("../data/"+FAM+"colored/"+FAM+"_colored.png",w=183,units='mm', tree_style=ts,dpi=300)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # =============================================================================
# 
# =============================================================================
