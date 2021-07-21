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
import re
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
import pprint

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
# TREEPATH = "../data/SPIKEtree/clustered/2021-02-02_SPIKE09_MP.nwk"
# =============================================================================

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
# # OUTTREE = "../data/SPIKEtree/clustered/2021-02-02_SPIKE09_ML.nwk"
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
# TREEPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle_relaxed.tree"
# OUTTREE = "../data/BSHtree/clustered/2021-07-13_BSH_filterd_BEAST_umapemm.nwk"
# FAM = 'BSH'
# TREEID='A'
# HEIGHT = 22
# =============================================================================
# =============================================================================
# 
# TREEPATH = "../data/OMPtree/2021-06-16_OMPuniprot_muscle.tree"
# OUTTREE = "../data/OMPtree/2021-07-13_OMPuniprot_muscle_umapemm.tree"
# TREEID='A'
# FAM = 'OMP'
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
# 
# TREEPATH = "../data/pfamtrees/2021-06-18_PF00271_seed.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-06-20_PF00271_seed_umapemm.tree"
# TREEID='PF00271'
# FAM = 'pfam'
# HEIGHT = 22
# 
# TREEPATH = "../data/pfamtrees/2021-06-18_PF07690_seed.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-06-18_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../data/pfamtrees/2021-07-02_PF00067_seed.1111.median.tree"
# OUTTREE = "../data/pfamtrees/clustered/2021-07-13_PF00067_seed_umapemm.tree"
# TREEID='PF00067'
# FAM = 'pfam'
# HEIGHT = 22
# 
# 
# TREEPATH = "../data/controltrees/2021-06-22_acyltransferase-filtered-homo-musclecorr.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-13_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
# TREEID='acyltrans'
# FAM = 'control'
# HEIGHT = 22
# =============================================================================


# =============================================================================
# TREEPATH = "../data/controltrees/2021-06-22_receptor-filtered-homo-musclecorr.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-13_receptor-filtered-homo-muscle_umapemm.mean.tree"
# TREEID='receptor'
# FAM = 'control'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# 
# TREEPATH = "../data/controltrees/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.1111.median.tree"
# OUTTREE = "../data/controltrees/clustered/2021-07-13_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.median.tree"
# TREEID='dehydrogenase' 
# FAM = 'control'
# HEIGHT = 22 
# =============================================================================

# =============================================================================
# TREEPATH = "../data/pfamtrees/PF02518.nhx"
# OUTTREE = "../data/pfamtrees/clustered/PF02518_umapgmmem.nhx"
# FAM = 'pfam'
# HEIGHT = 22
# 
# TREEPATH = "../data/pfamtrees/PF18885.nhx"
# OUTTREE = "../data/pfamtrees/clustered/PF18885_umapgmmem.nhx"
# FAM = 'pfam'
# HEIGHT = 22
# #boot for loop 
# TREEPATH = "../data/pfamtrees/PF18889.nhx"
# OUTTREE = "../data/pfamtrees/clustered/PF18889_umapgmmem.nhx"
# FAM = 'pfam'
# HEIGHT = 22
# #boot for loop 
# TREEPATH = "../data/pfamtrees/PF18886.nhx"
# OUTTREE = "../data/pfamtrees/clustered/PF18886_umapgmmem.nhx"
# FAM = 'pfam'
# HEIGHT = 22
# #boot for loop
# TREEPATH = "../data/pfamtrees/PF18886.nhx"
# OUTTREE = "../data/pfamtrees/clustered/PF18886_umapgmmem.nhx"
# FAM = 'pfam'
# HEIGHT = 22
# #boot for loop
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
# #boot for loop 
# TREEPATH = "../data/nextstraintrees/nextstrain_mumps_na_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_mumps_na_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# HEIGHT = 22
# TREEID="mumps"
# #boot for loop 
# =============================================================================

TREEPATH = "../data/nextstraintrees/nextstrain_ncov_global_timetree.nwk"
OUTTREE = "../data/nextstraintrees/clustered/nextstrain_ncov_global_timetree_umapgmmem.nwk"
FAM = 'nextstrain'
TREEID='ncov'
HEIGHT = 22
#boot for loop 

# =============================================================================
# TREEPATH = "../data/nextstraintrees/nextstrain_WNV_NA_timetree.nwk"
# OUTTREE = "../data/nextstraintrees/clustered/nextstrain_WNV_NA_timetree_umapgmmem.nwk"
# FAM = 'nextstrain'
# TREEID='WNV_NA'
# HEIGHT = 22
# #boot for loop 
# =============================================================================

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

print(TREEPATH[-4:])
if (TREEPATH[-4:] == '.nhx') | (TREEPATH[-4:] == '.nwk'):
    treeOR = dendropy.Tree.get(
        path=TREEPATH,
        schema="newick",
        extract_comment_metadata=True)
    
    ds = dendropy.DataSet.get_from_path(TREEPATH,
                                        "newick",
                                        extract_comment_metadata=True)
else:
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
# =============================================================================
#     if (TREEPATH[-4:] == '.nhx'):
#         tree.reroot_at_midpoint(update_bipartitions=False)
# =============================================================================

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
#print(pdm.as_data_table())


pdmA = np.zeros((len(tree.taxon_namespace),len(tree.taxon_namespace)))
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


print(pdmA)

N=5
if FAM == 'PRACT':
    N=5




metric = 'minkowski'
metrics = ['manhattan','minkowski','canberra','braycurtis',#'haversine',
           'mahalanobis','wminkowski','cosine']
n_Hs = [int(len(taxons)*3/4)]#,100,125,150,200]# [20,75,100,222]#[20,100,150,180,220,1000,2000,3000]



#%matplotlib inline
#sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

penguins = pdmA# pd.read_csv("https://github.com/allisonhorst/palmerpenguins/raw/5b5891f01b52ae26ad8cb9755ec93672f49328a8/data/penguins_size.csv")


for n_H in n_Hs:
    print("metric:",metric)
    reducer = umap.UMAP(n_components=N,
                        n_neighbors=n_H,
                        random_state=222,
                        metric = 'precomputed')
    
    
    penguin_data = pdmA
    
    embedding = reducer.fit_transform(penguin_data)
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
# =============================================================================
#     for i, txt in enumerate(taxons):
#         print(txt, embedding[i,0],embedding[i,1])
#         plt.annotate(txt, (float(embedding[i,0]),float(embedding[i,1])),size=15)
# =============================================================================
    plt.title('UMAP projection of the PDM', fontsize=24)
    plt.show()
    
    labels = ['PC' + str(x) for x in range(0,N) ]
    
    df = pd.DataFrame(taxons,columns = ['target'])
    
    principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
    
    finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
       
    
    
    

    
    
    
    #######################CLUTSER
    
    X = embedding[:, :]
    
    print(pdmA,taxons)
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1,19)
    if (FAM == 'control'):
        n_components_range = range(1,26)
    if (FAM == 'OMP')|(FAM == 'nextstrain'):
        n_components_range = range(1,21)
    if FAM == 'PRACT':
        n_components_range = range(1,6)
    cv_types = ['tied','full']#'spherical', 'diag', 'tied',
    
    
# =============================================================================
#     delta =100
#     bic = 100
#     while delta > 15:
#         
# =============================================================================
    
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
                
    
    
    bic = np.array(bic)
    BPcolors = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7",
                 '#f46d43',"#abd9e9"]
                #'#a9a9a9',                '#ffffff']#, '#000000'
# =============================================================================
#                 'salmon','red','pink', 'orange',#'steelblue',#'cornflowerblue', 'turquoise', 'darkorange',
#                 'olive','teal','purple','gold',
#                 'cyan','tan','green',
#                 'silver','maroon',
#                 '#7B68EE','#87CEEB','#DDA0DD','brown','#D2691E','lime']#,'green','grey']
# =============================================================================
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
# =============================================================================
#     for i, txt in enumerate(taxons):
#         print(txt, embedding[i,0],embedding[i,1])
#         plt.annotate(txt, (float(embedding[i,0]),float(embedding[i,1])),size=15)
# =============================================================================
    #plt.xticks(())
    #plt.yticks(())
    plt.title('Selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
    plt.subplots_adjust(hspace=.35, bottom=.02)
    plt.savefig('../data/EMMplots/2021-07-13_'+str(n_H)+'_GMMsweep_'+FAM+TREEID+'_'+metric+'.png')
    plt.show()
    
    
    
    
  
    print(clf.predict(X))
    print(clf.predict_proba(X))
    for i in range(len(taxons)):
        print(i)
        print(taxons[i])
        print(clf.predict(X)[i])
        print(np.around(clf.predict_proba(X)[i],decimals=10))
        print(clf.predict_proba(X)[i].dtype)
    
    
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
        #print(key)
        strTree = strTree.replace(key[0],key[1])
        
    print(strTree)
    
    #print(re.sub('\)\d:', '\).:',  strTree))
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
        #print(k,v,BPcolors[k%len(BPcolors)])
        if len(v)==1:
            mrca = BPtree.common_ancestor(v)#[i], v[i+1])
            mrca.color = BPcolors[k%len(BPcolors)]
        for i in range(len(v)-1):
            mp = BPtree.is_monophyletic(v)
            try:
                mrca = BPtree.common_ancestor(v[i], v[i+1])
                mrca.color = BPcolors[k%len(BPcolors)]
            except:
                pass
        k+=1     
        
    
    with plt.rc_context({'lines.linewidth': 4}):
        matplotlib.rc('font', size=0.5)
        fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
        axes = fig.add_subplot(1, 1, 1)
        plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
        plt.rc('axes', titlesize=14)     # fontsize of the axes title
        plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
        plt.rc('figure', titlesize=18)   # fontsize of the figure title
        Phylo.draw(BPtree,axes=axes,do_show=False,branch_labels=None)
        plt.savefig("../data/"+FAM+"colored/2021-07-13_"+str(n_H)+"_EMClust_"+FAM+TREEID+"_"+metric+"_coloredtree.png")
        plt.show()

    #make target vector
    labs = []
    for ii in finalDf['target']:
        ii=ii.replace(' ','_')
        labs.append(BPcolors[int(repDict[ii].split('|')[-1])%len(BPcolors)])#print(ii,repDict[ii].split('|')[-1])
    finalDf['labels'] = labs
    print(finalDf)
    


    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(1,1,1) 
    
    ax.set_xlabel('Dimension 1 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
    ax.set_ylabel('Dimension 2 ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
    ax.set_title('UMAP of selected GMM: full model, '+ str(clf.n_components) +' Gaussian components')
#    ax.set_title('2 Component UMAP', fontsize = 20)
    
    targets = [FAM]
    colors = [1,25,3]
    marks = ['o','v','*']

    ax.legend(targets)
    ax.grid()
    #plt.scatter(finalDf["PC1"],finalDf['PC2'],c=list(finalDf['labels']))
    
    
    for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
                                               color_iter)):
# =============================================================================
#         v, w = linalg.eigh(cov)
# =============================================================================
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 22, color=color,marker='o')
    
    
    
    # Plot an ellipse to show the Gaussian component
# =============================================================================
#         angle = np.arctan2(w[0][1], w[0][0])
#         angle = 180. * angle / np.pi  # convert to degrees
#         v = 2. * np.sqrt(2.) * np.sqrt(v)
#         ell = mpl.patches.Ellipse(mean, v[0]*40, v[1]*40, 180. + angle, color=color)
#         ell.set_clip_box(splot.bbox)
#         ell.set_alpha(.45)
#         splot.add_artist(ell)
# =============================================================================
    
    
    
    
    
    
    
    
    
    
    plt.show()

#############################################################

# =============================================================================
#     fig = plt.figure(figsize = (8,8))
#     ax = plt.axes(projection = '3d')
#     ax.set_xlabel('Dimension 1: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
#     ax.set_ylabel('Dimension 2: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
#     ax.set_title('3 Dimension UMAP', fontsize = 20)
#     
#     targets = [FAM]
#     colors = ['r', 'g', 'b']
#     marks = ['o','v','*']
#     
#     plt.scatter(finalDf['PC0']
#                    , finalDf['PC1']
#                    , finalDf['PC2']#,s=20
#                    , color = list(finalDf['labels'])
#                    , marker = 'D')
#     ax.legend(targets)
#     ax.grid()
#     plt.show()
# =============================================================================
    



    eteColors = ["","","","","","",""]

    eteTree = Tree(strTree[4:].replace("'",""))
   # eteTree.convert_to_ultrametric()
    print(eteTree)
    
    ts = TreeStyle()
    ts.show_leaf_name = True
    #ts.mode = "c"
    ts.root_opening_factor = 1

    ts.min_leaf_separation = 1
    ts.optimal_scale_level = 10
    ts.layout_fn = layout
    ts.show_leaf_name = False
    print(dir(ts))
# =============================================================================
#     ts.arc_start = -180 # 0 degrees = 3 o'clock
#     ts.arc_span = 180
# =============================================================================
    cc=0
    for key in catTax:
        print(key)
        nst = NodeStyle()
        color = BPcolors[cc%len(BPcolors)]
# =============================================================================
#         nst["vt_line_color"] = color
#         nst["hz_line_color"] = color
# =============================================================================
        nst["vt_line_width"] = 12
        nst["hz_line_width"] = 12
        nst["bgcolor"] = color
        print(nst)
        
        for k in key:
            print("KEY:",k.replace(' ','_'))
            
            ke='|'.join(k.replace(' ','_').split('|')[:-1])
            leaf = eteTree.search_nodes(name=ke)
            leafb = eteTree.search_nodes(name="'"+ke+"'")
            leafc = eteTree.search_nodes(name=k.replace(' ','_'))
            #print(leaf)
            if leaf != []:
                #nst["fsize"] = 14
                leaf[0].set_style(nst)

                for v in key:
                    #print(v)
                    #v='|'.join(v.replace(' ','_').split('|')[:-1])
                    try:
                        val='|'.join(v.replace(' ','_').split('|')[:-1])
                        node = eteTree.get_common_ancestor(ke,val)
                    except:
                        v=v
                        node = eteTree.get_common_ancestor(ke,v)
                    node.set_style(nst)
            elif leafb != []:
                leafb[0].set_style(nst)
                for v in key:
                    #print(v)
                    node = eteTree.get_common_ancestor("'"+ke+"'",v)
                    node.set_style(nst)
                    
            elif leafc != []: 
                leafc[0].set_style(nst)
                for v in key:
                    #print(v)
                    try:
                        node = eteTree.get_common_ancestor(k,v)
                        
                    except:
                        break
                        val='|'.join(v.replace(' ','_').split('|')[:-1])
                        node = eteTree.get_common_ancestor(k,v)
                    node.set_style(nst)
            else:
                pass
        cc += 1
        
        
        
    
    eteTree.show(tree_style=ts)
    eteTree.render("../data/"+FAM+"colored/2021-07-13"+FAM+TREEID+"_colored.png",w=183,units='mm', tree_style=ts,dpi=300)






    print(clf.predict(X))
    print(clf.predict_proba(X))
    for i in range(len(taxons)):
        print(taxons[i])
        print(clf.predict(X)[i])


    dictpred = {'0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[],
               '9':[],'10':[],'11':[],'12':[],'13':[],'14':[],'15':[],'16':[],
               '17':[],'18':[],'19':[],'20':[],'21':[],'22':[],'23':[],'24':[],
               '25':[]}
    for i in range(len(taxons)):
        #print(matrix.alphabet[i],clf.predict(X)[i])'
        dictpred[str(clf.predict(X)[i])].append(taxons[i]) 
    
    pprint.pprint(dictpred)

    OUTF="../data/"+FAM+"colored/2021-07-13"+FAM+TREEID+"_colored.csv"
    print(OUTF)
    import csv
    with open(OUTF,'w') as f :
        writer=csv.writer(f)
        for k,v in dictpred.items():
            for val in v:
                writer.writerow([k,val])


