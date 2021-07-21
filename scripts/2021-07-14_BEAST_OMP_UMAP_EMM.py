"""
Created on Thu Jan  7 16:13:27 2021

@author: aorti

Description

Usage

Arguments


"""

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

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=25)
        faces.add_face_to_node(N, node, 0, position="branch-right")



ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"
ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
TREEPATH = "../data/viraltrees/2021-06-23_sarscov2sgly_muscle_x.median.tree"
OUTTREE = "../data/viraltrees/clustered/2021-06-28_sarscov2sgly_umapemm.nwk"
TREEID='sarscov2sgly'
FAM = 'viral'
HEIGHT = 22

# =============================================================================
# =============================================================================
# # TREEPATH = "../data/SPIKEtree/2021-01-24_SPIKE09_btsp50_ML.nwk"
# # OUTTREE = "../data/SPIKEtree/clustered/2021-07-19_SPIKE09_ML.nwk"
# # FAM = 'SPIKE'
# # HEIGHT = 22
# # =============================================================================
# 
# 
# # =============================================================================
# # TREEPATH = "../../Taxon_Enzyme_Phylogeny/BSH_data/2020-03-05_MaxParsimony_BSH_id.nwk" 
# =============================================================================
# =============================================================================

TREEPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle_relaxed.tree"
OUTTREE = "../data/BSHtree/clustered/2021-07-19_BSH_filterd_BEAST_umapemm.nwk"
FAM = 'BSH'
TREEID='A'
HEIGHT = 22
# =============================================================================
# 
# TREEPATH = "../data/OMPtree/2021-06-16_OMPuniprot_muscle.tree"
# OUTTREE = "../data/OMPtree/2021-07-19_OMPuniprot_muscle_umapemm.tree"
# TREEID='A'
# FAM = 'OMP'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# 
# TREEPATH = "../data/pfamtrees/2021-06-14_PF00005_seed_strctclk.tree"
# OUTTREE = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# TREEID='PF00005'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

TREEPATH = "../data/pfamtrees/2021-06-18_PF00271_seed.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-06-20_PF00271_seed_umapemm.tree"
TREEID='PF00271'
FAM = 'pfam'
HEIGHT = 22

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-06-18_PF07690_seed.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-02_PF00067_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00067_seed_umapemm.tree"
# TREEID='PF00067'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00501_seed_umapemm.tree"
TREEID='PF00501'
FAM = 'pfam'
HEIGHT = 22

TREEPATH = "../data/pfamtrees/2021-07-15_PF00528_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00528_seed_umapemm.tree"
TREEID='PF00528'
FAM = 'pfam'
HEIGHT = 22

TREEPATH = "../data/pfamtrees/2021-07-15_PF00069_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00069_seed_umapemm.tree"
TREEID='PF00069'
FAM = 'pfam'
HEIGHT = 22

TREEPATH = "../data/pfamtrees/2021-07-15_PF00072_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00072_seed_umapemm.tree"
TREEID='PF00072'
FAM = 'pfam'
HEIGHT = 22

TREEPATH = "../data/pfamtrees/2021-07-15_PF00528_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00528_seed_umapemm.tree"
TREEID='PF00528'
FAM = 'pfam'
HEIGHT = 22

TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF00501_seed_umapemm.tree"
TREEID='PF00501'
FAM = 'pfam'
HEIGHT = 22
# =============================================================================
# 
# TREEPATH = "../data/pfamtrees/2021-07-15_PF02518_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-07-19_PF02518_seed_umapemm.tree"
# TREEID='PF02518'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================
# =============================================================================
# 
# TREEPATH = "../data/controltrees/2021-06-22_acyltransferase-filtered-homo-musclecorr.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-19_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
# TREEID='acyltrans'
# FAM = 'control'
# HEIGHT = 22
# 
# =============================================================================

# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_receptor-filtered-homo-musclecorr.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-19_receptor-filtered-homo-muscle_umapemm.mean.tree"
# TREEID='receptor'
# FAM = 'control'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-19_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.median.tree"
# TREEID='dehydrogenase' 
# FAM = 'control'
# HEIGHT = 22 
# =============================================================================


# =============================================================================
# 
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
# TREEPATH = "../data/nextstraintrees/nextstrain_ncov_global_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_ncov_global_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# TREEID='ncov'
# HEIGHT = 22
# =============================================================================


if (TREEPATH[-4:] == '.nhx') | (TREEPATH[-4:] == '.nwk'):
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

start_time = time.time()

pdmA = np.zeros((len(tree.taxon_namespace),len(tree.taxon_namespace)))
counter = 0
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon)[1:-1])
    for taxon2 in tree.taxon_namespace:
        if j>i:
            pdmA[i,j] =  float(pdm.distance(taxon,taxon2))
            pdmA[j,i] =  float(pdm.distance(taxon,taxon2))
# =============================================================================
#             if pdmA[i,j] < 0:
#                 pdmA[i,j] = abs(pdmA[i,j])
#                 pdmA[j,i] = abs(pdmA[j,i])
# =============================================================================
 
        j+=1
    i+=1


print(i,j)
print("--- %s seconds ---" % (time.time() - start_time))

N=7
if FAM == 'PRACT':
    N=5




metric = 'precomputed'

n_Hs = [int(len(taxons)*3/4)]


for n_H in n_Hs:
    print("metric:",metric)
    reducer = umap.UMAP(n_components=N,
                        n_neighbors=n_H,
                        random_state=222,
                        metric = 'precomputed')
    
    
    embedding = reducer.fit_transform(pdmA)
    
# =============================================================================
#     plt.scatter(
#         embedding[:, 0],
#         embedding[:, 1])
#     plt.ylabel("UMAP Dimension 2")
#     plt.xlabel("UMAP Dimension 1")
#     plt.title('UMAP projection of the PDM', fontsize=24)
#     plt.show()
# =============================================================================
    
    labels = ['PC' + str(x) for x in range(0,N) ]
    
    df = pd.DataFrame(taxons,columns = ['target'])
    
    principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
    
    finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
       
    #######################CLUTSER
    
    X = embedding[:, :]
    
    start = 1
    
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1,30)
    if (FAM == 'control'):
        n_components_range = range(1,50)
    if (FAM == 'OMP')|(FAM == 'nextstrain')|(FAM == 'viral'):
        start = 6
        n_components_range = range(start,30)
    if FAM == 'PRACT':
        n_components_range = range(1,6)
    cv_types = ['full']#'tied',
    
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          random_state = 111,
                                          covariance_type=cv_type,
                                          reg_covar=10e-06)
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
    print("--- %s seconds ---" % (time.time() - start_time))

 
    bic = np.array(bic)
    BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7",
                 '#f46d43',"#abd9e9"]
    
    if (FAM == 'OMP')|(FAM == 'nextstrain')|(FAM == 'pfam')|(FAM == 'control')|(FAM == 'viral'):
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
    plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
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
        if not np.any(Y_ == i):
            continue
        
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
    

    plt.ylabel("UMAP Dimension 2")
    plt.xlabel("UMAP Dimension 1")
    
    plt.xlim(min(embedding[:,0])-0.2,max(embedding[:,0])+0.2)
    plt.ylim(min(embedding[:,1])-0.2,max(embedding[:,1])+0.2)
    plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
    plt.subplots_adjust(hspace=.35, bottom=.02)
    plt.savefig('../data/EMMplots/2021-07-19_'+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.tiff',dpi=300)
    plt.show()
    
    
    
    

    taxonsLabels = []
    #for tree label replacement
    repDict={}
    #for coloring tree
    catTax = []
    
    for i in range(clf.n_components):
        catTax.append([])#{0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[]}
    for i in range(0,len(taxons)):
        ####
        taxons[i]=taxons[i].replace(' ','_')
        taxonsLabels.append(taxons[i]+"|"+str(clf.predict(X)[i]))
        
        #for coloring tree
        catTax[int(clf.predict(X)[i])].append(taxons[i]+"|"+str(clf.predict(X)[i]))
        
        #for tree label replacement
        repDict[taxons[i]] = taxons[i]+"|"+str(clf.predict(X)[i])
    
    strTree = tree.as_string(schema="newick")
    print(strTree)
    strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9]','', strTree)
    
    for key in repDict.items():
        strTree = strTree.replace(key[0],key[1])
        
    print(strTree)
    
    finalStrTree = strTree
    for iter in re.finditer("\D(?<=:)",strTree):
        a = iter.span()
        b = a[1]
        a = a[0]
        if strTree[a-7] ==')':
            print(strTree[a-7:b])
            finalStrTree = finalStrTree.replace(strTree[a-7:b],'):')
        
    print(finalStrTree)
    EMtree_file = open(OUTTREE,'w+')
    n = EMtree_file.write(finalStrTree)
    EMtree_file.close()
    
    #########################################################################
    from Bio.Phylo.PhyloXML import Phylogeny
    
    BPtree = Phylo.read(OUTTREE, "newick")
    BPtree.ladderize()
     
    #make tree color-able
    BPtree.rooted = True
    BPtree = Phylogeny.from_tree(BPtree)
    
    #color tree
    BPtree.root.color = (128, 128, 128)
    Phylo.draw(BPtree)
    
    k=0
    for v in catTax:
        if len(v)==1:
            mrca = BPtree.common_ancestor(v)
            mrca.color = BPcolors[k%len(BPcolors)]
        for i in range(len(v)-1):
            mp = BPtree.is_monophyletic(v)
            print(mp,v)
            try:
                mrca = BPtree.common_ancestor(v[i], v[i+1])
                mrca.color = BPcolors[k%len(BPcolors)]
            except:
                pass
        k+=1     
        
    print(mp)
    
# =============================================================================
#     with plt.rc_context({'lines.linewidth': 4}):
#         matplotlib.rc('font', size=0.0)
#         fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
#         axes = fig.add_subplot(1, 1, 1)
#         plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
#         plt.rc('axes', titlesize=14)     # fontsize of the axes title
#         plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('figure', titlesize=18)   # fontsize of the figure title
#         Phylo.draw(BPtree,axes=axes,do_show=False,branch_labels=None)
#         plt.savefig("../data/"+FAM+"colored/2021-07-19_"+str(n_H)+"_EMClust_"+FAM+TREEID+"_"+metric+"_coloredtree.tiff",dpi=300)
#         plt.show()
# 
#     #make target vector
#     labs = []
#     for ii in finalDf['target']:
#         ii=ii.replace(' ','_')
#         labs.append(BPcolors[int(repDict[ii].split('|')[-1])%len(BPcolors)])#print(ii,repDict[ii].split('|')[-1])
#     finalDf['labels'] = labs
# 
# 
#     fig = plt.figure(figsize = (5,5))
#     ax = fig.add_subplot(1,1,1) 
#     
#     ax.set_xlabel('Dimension 1 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
#     ax.set_ylabel('Dimension 2 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
#     ax.set_title('UMAP of selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
#     
#     targets = [FAM]
#     colors = [1,25,3]
#     marks = ['o','v','*']
# 
#     ax.legend(targets)
#     ax.grid()
#     
#     for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
#                                                color_iter)):
#         if not np.any(Y_ == i):
#             continue
#         plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
#     
#     plt.show()
# 
# #############################################################
#     eteColors = ["","","","","","",""]
# 
#     eteTree = Tree(strTree[4:].replace("'",""))
#     print(eteTree)
#     
#     ts = TreeStyle()
#     ts.show_leaf_name = True
#     ts.root_opening_factor = 1
# 
#     ts.min_leaf_separation = 1
#     ts.optimal_scale_level = 10
#     ts.layout_fn = layout
#     ts.show_leaf_name = False
#     print(dir(ts))
#     cc=0
#     for key in catTax:
#         print(key)
#         nst = NodeStyle()
#         color = BPcolors[cc%len(BPcolors)]
#         nst["vt_line_width"] = 12
#         nst["hz_line_width"] = 12
#         nst["bgcolor"] = color
#         print(nst)
#         
#         for k in key:
#             print("KEY:",k.replace(' ','_'))
#             
#             ke='|'.join(k.replace(' ','_').split('|')[:-1])
#             leaf = eteTree.search_nodes(name=ke)
#             leafb = eteTree.search_nodes(name="'"+ke+"'")
#             leafc = eteTree.search_nodes(name=k.replace(' ','_'))
#             if leaf != []:
#                 leaf[0].set_style(nst)
# 
#                 for v in key:
#                     try:
#                         val='|'.join(v.replace(' ','_').split('|')[:-1])
#                         node = eteTree.get_common_ancestor(ke,val)
#                     except:
#                         v=v
#                         node = eteTree.get_common_ancestor(ke,v)
#                     node.set_style(nst)
#             elif leafb != []:
#                 leafb[0].set_style(nst)
#                 for v in key:
#                     #print(v)
#                     node = eteTree.get_common_ancestor("'"+ke+"'",v)
#                     node.set_style(nst)
#                     
#             elif leafc != []: 
#                 leafc[0].set_style(nst)
#                 for v in key:
#                     try:
#                         node = eteTree.get_common_ancestor(k,v)
#                         
#                     except:
#                         break
#                         val='|'.join(v.replace(' ','_').split('|')[:-1])
#                         node = eteTree.get_common_ancestor(k,v)
#                     node.set_style(nst)
#             else:
#                 pass
#         cc += 1
#         
#         
#         
#     
#     eteTree.show(tree_style=ts)
#     eteTree.render("../data/"+FAM+"colored/2021-07-19"+FAM+TREEID+"_colored.tiff",w=183,units='mm', tree_style=ts,dpi=300)
# 
# 
# 
#     dictpred = {}
#     
#     for i in range(50):
#         dictpred[str(i)] = []
#     for i in range(len(taxons)):
#         #print(matrix.alphabet[i],clf.predict(X)[i])'
#         dictpred[str(clf.predict(X)[i])].append(taxons[i]) 
#     
#     pprint.pprint(dictpred)
# 
#     OUTF="../data/"+FAM+"colored/2021-07-19"+FAM+TREEID+"_colored.csv"
#     print(OUTF)
#     import csv
#     with open(OUTF,'w') as f :
#         writer=csv.writer(f)
#         for k,v in dictpred.items():
#             for val in v:
#                 writer.writerow([k,val])
#  
# 
# 
# 
# =============================================================================
