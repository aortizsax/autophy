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

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=25)
        faces.add_face_to_node(N, node, 0, position="branch-right")



ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"
ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
TREEPATH = "../data/viraltrees/2021-06-23_sarscov2sgly_muscle_x.median.tree"
OUTTREE = "../data/viraltrees/clustered/2021-08-04_sarscov2sgly_umapemm.nwk"
TREEID='sarscov2sgly'
FAM = 'viral'
HEIGHT = 22

# =============================================================================
# =============================================================================
# # TREEPATH = "../data/SPIKEtree/2021-01-24_SPIKE09_btsp50_ML.nwk"
# # OUTTREE = "../data/SPIKEtree/clustered/2021-08-04_SPIKE09_ML.nwk"
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
OUTTREE = "../data/BSHtree/clustered/2021-08-04_BSH_filterd_BEAST_umapemm.nwk"
FAM = 'BSH'
TREEID='A'
HEIGHT = 22

TREEPATH = "../data/OMPtree/2021-06-16_OMPuniprot_muscle.tree"
OUTTREE = "../data/OMPtree/2021-08-04_OMPuniprot_muscle_umapemm.tree"
TREEID='A'
FAM = 'OMP'
HEIGHT = 22
# =============================================================================
# 
# TREEPATH = "../data/panthertrees/2021-08-02_PTHR10000_aln.2222.mean.tree"
# OUTTREE = "../data/panthertrees/2021-08-04_PTHR10000_aln.umapemm.tree"
# TREEID='PTHR1000'
# FAM = 'panther'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-06-14_PF00005_seed_strctclk.tree"
# OUTTREE = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# TREEID='PF00005'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-06-18_PF00271_seed.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00271_seed_umapemm.tree"
# TREEID='PF00271'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================
# =============================================================================
# 
# TREEPATH = "../data/pfamtrees/2021-06-18_PF07690_seed.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# ###seccatTax 
# TREEPATH = "../data/pfamtrees/2021-07-02_PF00067_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00067_seed_umapemm.tree"
# TREEID='PF00067'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00501_seed_umapemm.tree"
# TREEID='PF00501'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00528_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00528_seed_umapemm.tree"
# TREEID='PF00528'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ##fix arrary seccattax
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00069_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00069_seed_umapemm.tree"
# TREEID='PF00069'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00072_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00528_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00528_seed_umapemm.tree"
# TREEID='PF00528'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF00501_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF00501_seed_umapemm.tree"
# TREEID='PF00501'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-15_PF02518_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-08-04_PF02518_seed_umapemm.tree"
# TREEID='PF02518'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_acyltransferase-filtered-homo-musclecorr.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-08-04_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
# TREEID='acyltrans'
# FAM = 'control'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_receptor-filtered-homo-musclecorr.mean.tree"
# OUTTREE = "../data/controltrees/clustered/2021-08-04_receptor-filtered-homo-muscle_umapemm.mean.tree"
# TREEID='receptor'
# FAM = 'control'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.1111.mean.tree"
# OUTTREE = "../data/controltrees/clustered/2021-08-04_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.mean.tree"
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
            # if negative branch legths
# =============================================================================
#             if pdmA[i,j] < 0:
#                 pdmA[i,j] = abs(pdmA[i,j])
#                 pdmA[j,i] = abs(pdmA[j,i])
# =============================================================================
 
        j+=1
    i+=1


print("--- %s seconds ---" % (time.time() - start_time))

N=6
N=6
if FAM == 'PRACT':
    N=5




metric = 'precomputed'

n_Hs = [int(len(taxons)/2)]#int(len(taxons)*3/4),


for n_H in n_Hs:
    print("metric:",metric)
    reducer = umap.UMAP(n_components=N,
                        n_neighbors=n_H,
                        random_state=222,
                        metric = 'precomputed')
    
    
    embedding = reducer.fit_transform(pdmA)
    
    
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
    plt.text(xpos-1+start, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
    spl.set_xlabel('Number of Gaussian components')
    spl.legend([b[0] for b in bars], cv_types)
    plt.ylabel("BIC score")
    
    # Plot the winner
    color_iter = itertools.cycle(BPcolors)

    splot = plt.subplot(2, 1, 2)
    Y_ = clf.predict(X)
    
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
    plt.savefig('../data/EMMplots/2021-08-04_'+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.tiff',dpi=300)
    plt.show()
    
    
    
    

    taxonsLabels = []
    #for tree label replacement
    repDict={}
    #for coloring tree
    catTax = []

    
    for i in range(clf.n_components):
        catTax.append([])
    for i in range(0,len(taxons)):
        taxons[i]=taxons[i].replace(' ','_')
        taxonsLabels.append(taxons[i]+"|"+str(clf.predict(X)[i]))
        
        #for coloring tree
        catTax[int(clf.predict(X)[i])].append(taxons[i]+"|"+str(clf.predict(X)[i]))
        
        #for tree label replacement
        repDict[taxons[i]] = taxons[i]+"|"+str(clf.predict(X)[i])
    
    strTree = tree.as_string(schema="newick")
    
    strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9]','', strTree)
    
    for key in repDict.items():
        strTree = strTree.replace(key[0],key[1])
        
    finalStrTree = strTree
    for iter in re.finditer("\D(?<=:)",strTree):
        a = iter.span()
        b = a[1]
        a = a[0]
        if strTree[a-7] ==')':
            
            finalStrTree = finalStrTree.replace(strTree[a-7:b],'):')
    print(finalStrTree)
    EMtree_file = open(OUTTREE,'w+')
    n = EMtree_file.write(finalStrTree)
    EMtree_file.close()

                

    #########################################################################
    from Bio.Phylo.PhyloXML import Phylogeny
    
    BPtree = Phylo.read(OUTTREE, "newick")
    

    
    terminals = BPtree.get_terminals()
    termdict = {}
    flag=[]
    
    for terminal in terminals:
        print(terminal)
        termdict[str(terminal)] = terminal

    for v in catTax:
        i=0
        taxa=[0]*len(v)
        for tax in v:
            print(len(tax),len('STRR6|EnsemblGenome.spr1170|UniProtKB'))
            if len(tax) > 39:                
                taxa[i] = termdict[tax[:37]+'...']
            else: 
                taxa[i] = termdict[tax]
            i+=1 
        mp = BPtree.is_monophyletic(taxa)
        print(mp)
        if (mp)==False:
            flag.append(v)
            


    #####recluster
    reclustCount = 0
    reclustlist = {}
    n_comp = clf.n_components
    for f in flag:    #flag:
        reclustCount+=1
        pdm2=np.zeros((len(f),len(f)))
        i=0
        for f1 in f:
            j=0
            for f2 in f:
                if i>j:
                    pdm2[i,j]=BPtree.distance(f1,f2)
                    pdm2[j,i]=BPtree.distance(f1,f2)
                j+=1
            i+=1
            

            
        
        penguins = pdm2
        clu = f2.split('|')[-1]

        ## Start HeirClustering
        linked = linkage(penguins, 'single')
        labelList = f
        plt.figure(figsize=(16, 4))
        plt.title('Hierarchical Clustering Dendrogram (subtree)', fontsize=20)
        plt.xlabel('Distance', fontsize=16)
        plt.yticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
        dendro = dendrogram(linked,leaf_rotation=30,
                            leaf_font_size = 10,
                    orientation='right',
                    labels=labelList,
                    distance_sort='descending',
                    show_leaf_counts=True)
        from datetime import date
        import datetime 
        print(datetime.datetime.now())
        print(date.today())
        daytime = str(datetime.datetime.now()).replace(':','.').replace(' ','_')
        plt.savefig('../data/EMMplots/'+str(daytime)+'_PASS2_'+f2.split('|')[-1]+'_dendrogramHC_'+FAM+TREEID+'_'+metric+'.tiff',dpi=300)
        plt.show()
# =============================================================================
#         
#         col_list = dendro['color_list']
#         
#                 
#         from sklearn.cluster import AgglomerativeClustering
#         maxcluster = 0
#         for i in range(len(col_list)): 
#             tempcol = int(col_list[i].replace('C','')) 
#             if tempcol > maxcluster: 
#                 maxcluster = tempcol
#         hc = AgglomerativeClustering(n_clusters=maxcluster+1, affinity = 'euclidean', linkage = 'single')
#         y_hc=hc.fit_predict(pdm2)
#         # print(dendro)        
#         # print(y_hc)
#         # print(f)
#         
#         #######################GMM UMAP
# # =============================================================================
# #         N_2 = 3
# #         reducer = umap.UMAP(n_components=N_2,
# #                             n_neighbors=n_H,
# #                             random_state=222,
# #                             metric = 'precomputed')
# #         
# #         
# #         embedding = reducer.fit_transform(pdm2)
# #         
# #         
# #         labels = ['PC' + str(x) for x in range(0,N_2) ]
# #         
# #         df = pd.DataFrame(taxons,columns = ['target'])
# #         
# #         principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
# #         
# #         finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
# #            
# #         #######################CLUTSER
# #         
# #         X = embedding[:, :]
# #         
# #         start = 1
# #         
# #         lowest_bic = np.infty
# #         bic = []
# #         n_components_range = range(1,int(0.75 * len(X)))
# # 
# #         cv_types = ['spherical','diag','tied','full']#'tied',
# #         
# #         for cv_type in cv_types:
# #             for n_components in n_components_range:
# #                 # Fit a Gaussian mixture with EM
# #                 gmm = mixture.GaussianMixture(n_components=n_components,
# #                                               random_state = 111,
# #                                               covariance_type=cv_type,
# #                                               reg_covar=10e-06)
# #                 gmm.fit(X)
# #                 bic.append(gmm.bic(X))
# #                 if bic[-1] < lowest_bic:
# #                     lowest_bic = bic[-1]
# #                     best_gmm = gmm
# #         print("--- %s seconds ---" % (time.time() - start_time))
# #     
# #      
# #         bic = np.array(bic)
# #         BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
# #                      "#D55E00", "#CC79A7",
# #                      '#f46d43',"#abd9e9"]
# #         
# #         if (FAM == 'OMP')|(FAM == 'nextstrain')|(FAM == 'pfam')|(FAM == 'control')|(FAM == 'viral'):
# #             BPcolors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
# #                     '#ba3fdf', '#42d4f4', '#bfef45', '#fabed4', 
# #                     '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', 
# #                     '#aaffc3', '#808000', '#ffd8b1', '#0072B2']
# #         color_iter = itertools.cycle(BPcolors)
# #         clf = best_gmm
# #         bars = []
# #         
# #         # Plot the BIC scores
# #         plt.figure(figsize=(6,6))
# #         spl = plt.subplot(2, 1, 1)
# #         for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
# #             xpos = np.array(n_components_range) + .2 * (i - 2)
# #             bars.append(plt.bar(xpos, bic[i * len(n_components_range):
# #                                           (i + 1) * len(n_components_range)],
# #                                 width=.2, color=color))
# #         plt.xticks(n_components_range)
# #         plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
# #         plt.title('BIC score per model')
# #         xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
# #             .2 * np.floor(bic.argmin() / len(n_components_range))
# #         plt.text(xpos-1+start, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
# #         spl.set_xlabel('Number of Gaussian components')
# #         spl.legend([b[0] for b in bars], cv_types)
# #         plt.ylabel("BIC score")
# #         
# #         # Plot the winner
# #         color_iter = itertools.cycle(BPcolors)
# #     
# #         splot = plt.subplot(2, 1, 2)
# #         Y_ = clf.predict(X)
# #         
# #         for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
# #                                                    color_iter)):
# #             if not np.any(Y_ == i):
# #                 continue
# #             
# #             plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
# #         
# #     
# #         plt.ylabel("UMAP Dimension 2")
# #         plt.xlabel("UMAP Dimension 1")
# #         
# #         plt.xlim(min(embedding[:,0])-0.2,max(embedding[:,0])+0.2)
# #         plt.ylim(min(embedding[:,1])-0.2,max(embedding[:,1])+0.2)
# #         plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
# #         plt.subplots_adjust(hspace=.35, bottom=.02)
# #         plt.savefig('../data/EMMplots/2021-08-04_PASS2'+f2.split('|')[-1]+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.tiff',dpi=300)
# #         plt.show()
# #         
# # =============================================================================
#             
#         secondpass = {}
#         maxcluster = 0
#         
#         for i in range(len(f)):
#             print(y_hc[i],clf.predict(X)[i], f[i])
#     
#         
#         for i in range(len(y_hc)): 
#             tempcol = y_hc[i]
#             print(tempcol,f[i]) 
#             if tempcol > maxcluster: 
#                 maxcluster = tempcol 
#         n_comp += maxcluster  
#             
#         for i in range(len(y_hc)):
#             col = str(y_hc[i])
#             tax = f[i]
#             lab = tax+'.'+col.replace('C','')
#             
#             secondpass[tax] = lab
# 
#                 
#         for k,v in secondpass.items():
#             finalStrTree = finalStrTree.replace(k,v)
#         reclustlist[str(clu)] = secondpass
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
#     print(clusthash)
#         
#     
# #########   WRITE TREE
#     EMtree_file = open(OUTTREE,'w+')
#     n = EMtree_file.write(finalStrTree)
#     EMtree_file.close()       
# 
#     
# 
#     ###OPEN reclustBPtree
#     from io import StringIO
# 
#     handle = StringIO(finalStrTree)
#     reclustBPtree = Phylo.read(OUTTREE,"newick")
#     reclustBPtree.ladderize()
#      
#     #make tree color-able
#     reclustBPtree.rooted = True
#     reclustBPtree = Phylogeny.from_tree(BPtree)
#     
#     #color tree
#     reclustBPtree.root.color = (128, 128, 128) 
# 
#             
#     mps=[]
#     k=0
#     for v in clusthash.values():
#         if len(v)==1:
#             try:
#                 mrca = reclustBPtree.common_ancestor(v)
#                 mrca.color = BPcolors[k%len(BPcolors)]
#             except:
#                 pass
#                 mrca = reclustBPtree.common_ancestor('.'.join(v[0].split('.')[:-1]))
#                 mrca.color = BPcolors[k%len(BPcolors)]
#         for i in range(len(v)-1):
# 
#             print(v[i])
#             try:
#                 mrca = reclustBPtree.common_ancestor(v[i], v[i+1])
#                 mrca.color = BPcolors[k%len(BPcolors)]
#             except ValueError:
#                 # Phylo.draw_ascii(reclustBPtree)  
#                 print((reclustBPtree))
#                 print(finalStrTree) 
#                 print(v[i], v[i+1])
#                 
#                 mrca = reclustBPtree.common_ancestor('.'.join(v[i].split('.')[:-1]), 
#                                                      '.'.join(v[i+1].split('.')[:-1]))
#                 mrca.color = BPcolors[k%len(BPcolors)]
#                 pass
#         k+=1     
# 
#     
#     with plt.rc_context({'lines.linewidth': 4}):
#         matplotlib.rc('font', size=0.0)
#         fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
#         axes = fig.add_subplot(1, 1, 1)
#         plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
#         plt.rc('axes', titlesize=14)     # fontsize of the axes title
#         plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
#         plt.rc('figure', titlesize=18)   # fontsize of the figure title
#         Phylo.draw(reclustBPtree,axes=axes,do_show=False,branch_labels=None)
#         plt.savefig("../data/"+FAM+"colored/2021-08-04_"+str(n_H)+"_EMClust_"+FAM+TREEID+"_"+metric+"_coloredtree.tiff",dpi=300)
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
# 
#     plt.show()
# 
# 
# 
# #############################################################
#     eteColors = ["","","","","","",""]
# 
#     eteTree = Tree(finalStrTree[4:].replace("'",""))
#     print(eteTree)
#     
#     ts = TreeStyle()
#     ts.show_leaf_name = True
#     
#     ts.mode = "c"
#     ts.arc_start = -180 # 0 degrees = 3 o'clock
#     ts.arc_span = 180
#     
#     ts.root_opening_factor = 1
# 
#     ts.min_leaf_separation = 1
#     ts.optimal_scale_level = 10
#     ts.layout_fn = layout
# # =============================================================================
# #     ts.show_leaf_name = False
# # =============================================================================
#     print(dir(ts))
#     cc=0
#     for key in clusthash.values():
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
# # =============================================================================
# #     eteTree.show(tree_style=ts)
# #     eteTree.render("../data/"+FAM+"colored/2021-08-04"+FAM+TREEID+"_colored.tiff",w=183,units='mm', tree_style=ts,dpi=300)
# # =============================================================================
# 
# 
# 
# 
# 
#     ### Change for secondpass
# 
#     
# # =============================================================================
# #     for i in range(50):
# #         dictpred[str(i)] = []
# #     for i in range(len(taxons)):
# #         print(clf,i)
# #         #print(matrix.alphabet[i],clf.predict(X)[i])'
# #         dictpred[str(clf.predict(X)[i])].append(taxons[i]) 
# # =============================================================================
#     dictpred = {}
#     for clustlist in seccatTax:
#         if len(clustlist)>0:
#             print(clustlist)
#             key = str(clustlist[0].split('|')[-1])
#             print(key)
#             dictpred[key] = []
#             for val in clustlist:
#                 print(val)
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
#    ################################################################################# 
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     ####################################################################################
#     pprint.pprint(dictpred)
# 
#     OUTF="../data/"+FAM+"colored/2021-08-04"+FAM+TREEID+"_colored.csv"
#     print(OUTF)
#     import csv
#     with open(OUTF,'w') as f :
#         writer=csv.writer(f)
#         writer.writerow(['Label','Cluster','Color'])
#         i = 0
#         for k,v in dictpred.items():
#             for val in v:
#                 print(k,val)
#                 writer.writerow([val,k,BPcolors[i%len(BPcolors)]])
#             i+=1
#  
# 
# 
# 
# print(clusthash)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# =============================================================================
