# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 16:13:27 2021

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
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap


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
TREEPATH = "../../Taxon_Enzyme_Phylogeny/BSH_data/2020-03-05_MaxParsimony_BSH_id.nwk" 
TREEPATH = "../data/OMPtree/2021-01-15_BSH_filterd_ML.nwk"
OUTTREE = "../data/2021-01-05_MP_BSH_idShort.nwk"
FAM = 'BSH'
HEIGHT = 22
# =============================================================================
# TREEPATH = "../../PhyloClust/data/2020-04-06_MaxParsimony_BSH_idFull.nwk" 
# TREEPATH = "../data/OMPtree/uniprot-gene_omp+reviewed_yes-8924-5104.nwk"#"2021_ML_OMP.nwk" 
# OUTTREE = "../data/2021-01-14_MP_OMP.nwk"
# FAM = 'OMP'
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
pdm = treeOR.phylogenetic_distance_matrix()
#print(pdm.as_data_table())


pdmA = np.zeros((len(treeOR.taxon_namespace),len(treeOR.taxon_namespace)+2))
counter = 0
print(pdmA.shape)


i = 0
taxons = []
for taxon in treeOR.taxon_namespace:
    j = 0
    taxons.append(str(taxon)[1:-1])
   
    for taxon2 in treeOR.taxon_namespace:
        edgecount = float(pdm.path_edge_count(taxon,taxon2))
        pdmA[i,j] = edgecount#float(pdm.distance(taxon,taxon2)) 
                
        j+=1
    i+=1
    
    




taxonsLabels = []
repDict={}

for i in range(0,len(taxons)):
    taxonsLabels.append(taxons[i].split(' ')[0])
    #print(clf.predict(X)[i],taxons[i],BPcolors[clf.predict(X)[i]%len(BPcolors)]
    #for tree label replacement
    repDict[taxons[i].replace(' ', '_')] = taxons[i].split(' ')[0].split('_')[0]
    #print(i,taxons[i],repDict[taxons[i]])
strTree = treeOR.as_string(schema="newick")
strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
strTree = re.sub('Inner[0-9][0-9]','', strTree)
strTree = re.sub('Inner[0-9]','', strTree)

for key in repDict.items():
    print(key)
    strTree = strTree.replace(key[0],key[1])
    strTree = strTree.replace(key[0].replace('_',' '),key[1])
    

EMtree_file = open(OUTTREE.replace('.nwk', 'shortid.nwk'),'w+')
n = EMtree_file.write(strTree)
EMtree_file.close()


tree = dendropy.Tree.get(
    path=OUTTREE.replace('.nwk', 'shortid.nwk'),
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
tree.ladderize()
pdm = tree.phylogenetic_distance_matrix()
#print(pdm.as_data_table())


pdmA = np.zeros((len(tree.taxon_namespace),len(tree.taxon_namespace)+3))
counter = 0
print(pdmA.shape)
i = 0
taxons = []
for taxon in tree.taxon_namespace:
    j = 0
    taxons.append(str(taxon)[1:-1].split(' ')[0])
    for taxon2 in tree.taxon_namespace:
        edgecount = float(pdm.path_edge_count(taxon,taxon2))
        pdmA[i,j] =  float(pdm.distance(taxon,taxon2))
                
        j+=1
    i+=1
bootcount = 0
i += 1
for node in treeOR.postorder_node_iter():#preorder_node_iter():
    if node.label != None:
        if float(node.label) > 0.75:
            #print(str(node.label))
            bootcount += 1
    if node.taxon != None:
        print(bootcount,str(node.taxon).split(' ')[0])
        index = taxons.index(str(node.taxon)[1:-1].split(' ')[0])
        print(index,i)
        pdmA[index,i] = bootcount
i += 1
for node in treeOR.preorder_node_iter():#preorder_node_iter():    if node.label != None:
    if node.label != None:
        if float(node.label) > 0.75:
            #print(str(node.label))
            bootcount += 1
    if node.taxon != None:
        print(bootcount,str(node.taxon).split(' ')[0])
        index = taxons.index(str(node.taxon)[1:-1].split(' ')[0])
        print(index,i)
        pdmA[index,i] = bootcount

print(pdmA)

N=10
#pdist = ssd.pdist(pdmA)
print(taxons)



metric = 'minkowski'
metrics = ['manhattan','minkowski','canberra','braycurtis',#'haversine',
           'mahalanobis','wminkowski','cosine']
n_Hs = [20,75,100,222]#[20,100,150,180,220,1000,2000,3000]



#%matplotlib inline

#sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

penguins = pdmA# pd.read_csv("https://github.com/allisonhorst/palmerpenguins/raw/5b5891f01b52ae26ad8cb9755ec93672f49328a8/data/penguins_size.csv")


#sns.pairplot(penguins, hue='species_short')

for n_H in n_Hs:
    print("metric:",n_H)
    reducer = umap.UMAP(n_components=N,    
                        n_neighbors=n_H,
                        min_dist=0.0,
                        random_state=42,
                        metric = metric)
    
    
    penguin_data = pdmA
    scaled_penguin_data = StandardScaler().fit_transform(penguin_data)
    
    embedding = reducer.fit_transform(scaled_penguin_data)
    #print(embedding.shape)
    
    plt.scatter(
        embedding[:, 0],
        embedding[:, 1])#,
        #c=[sns.color_palette()[x] for x in penguins.species_short.map({"Adelie":0, "Chinstrap":1, "Gentoo":2})])
    #plt.gca().set_aspect('equal', 'box')
    plt.title('UMAP projection of the Penguin dataset', fontsize=24)
    #plt.show()
    
    labels = ['PC' + str(x) for x in range(1,N+1) ]
    
    df = pd.DataFrame(taxons,columns = ['target'])
    
    principalDf = pd.DataFrame(data = embedding[:, :], columns = labels)
    
    finalDf = pd.concat([principalDf,df[['target']]], axis = 1)
       
    
    
    

    
    
    
    #######################CLUTSER
    
    X = embedding[:, :]
    
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
                                          covariance_type=cv_type,
                                          reg_covar=10e-06)
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
    #print(clf.covariances_.shape,clf.covariances_)
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
    plt.savefig('../data/EMMplots/2021-01-15a_'+str(n_H)+'_GMMsweep_'+FAM+'_'+metric+'.png')
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
        #print(clf.predict(X)[i],taxons[i],BPcolors[clf.predict(X)[i]%len(BPcolors)])
        #for coloring tree
        catTax[int(clf.predict(X)[i])].append(taxons[i]+"|"+str(clf.predict(X)[i]))
        #for tree label replacement
        repDict[taxons[i]] = taxons[i]+"|"+str(clf.predict(X)[i])
    
    strTree = tree.as_string(schema="newick")
    # =============================================================================
    print(strTree)
    # =============================================================================
    strTree = re.sub('Inner[0-9][0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9][0-9]','', strTree)
    strTree = re.sub('Inner[0-9]','', strTree)
    
    for key in repDict.items():
        #print(key)
        strTree = strTree.replace(key[0],key[1])
        
    # =============================================================================
    print(strTree)
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
        #print(k,v,BPcolors[k%len(BPcolors)])
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
        fig = plt.figure(figsize=(8,HEIGHT), dpi=200)
        axes = fig.add_subplot(1, 1, 1)
        plt.rc('font', size=6)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
        plt.rc('axes', titlesize=14)     # fontsize of the axes title
        plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
        plt.rc('figure', titlesize=18)   # fontsize of the figure title
        Phylo.draw(BPtree,axes=axes,do_show=False,branch_labels=None)
        plt.savefig("../data/"+FAM+"colored/2021-01-15a_"+str(n_H)+"_EMClust_"+FAM+"_"+metric+"_coloredtree.png")
        plt.show()

    #make target vector
    labs = []
    for ii in finalDf['target']:
        
        labs.append(BPcolors[int(repDict[ii].split('|')[-1])%len(BPcolors)])#print(ii,repDict[ii].split('|')[-1])
    finalDf['labels'] = labs
    print(finalDf)
    


    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    
    ax.set_xlabel('Component 1: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
    ax.set_ylabel('Component 2: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
    ax.set_title('2 Component UMAP', fontsize = 20)
    
    targets = ['OMP', 'PMP', 'PAL']
    colors = [1,25,3]
    marks = ['o','v','*']

    ax.legend(targets)
    ax.grid()
    plt.scatter(finalDf["PC1"],finalDf['PC2'],c=list(finalDf['labels']))
    plt.show()

#############################################################

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('Component 1: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[0], 2))
    ax.set_ylabel('Component 2: ', fontsize = 15)#+str(round(pca.explained_variance_ratio_[1], 2))
    ax.set_title('3 Components UMAP', fontsize = 20)
    
    targets = ['OMP', 'PMP', 'PAL']
    colors = ['r', 'g', 'b']
    marks = ['o','v','*']
    
    plt.scatter(finalDf['PC1']
                   , finalDf['PC2']
                   , finalDf['PC3']
                   , c = list(finalDf['labels']))
    ax.legend(targets)
    ax.grid()
    plt.show()
    


print(pdmA)
