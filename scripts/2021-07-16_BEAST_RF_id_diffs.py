# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:37:51 2021

@author: aorti
"""

import dendropy

import statistics 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from collections import OrderedDict

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc

from prettytable import PrettyTable


#open aln do random forest ################get pdist for location
# investigate gaps 

#make matrix 
#row taxons
#columns are locations of gaps


ALNPATH = "../data/viralaln/sarscov2sgly_muscle_x.nexus"
ALNPATH = "../data/viralaln/sarscov2sgly_muscle.fasta"

#TREEPATH = "../data/viraltrees/sarscov2sgly_muscle_x_cah.tree"
TREEPATH = "../data/viraltrees/clustered/2021-07-15_sarscov2sgly_umapemm.nwk"
N = 1300
FAM = 'viral'
TREEID='sarscov2sgly'
REP='YP '

ALNPATH = "../data/BSHaln/2021-03-21_BSH_filtered_muscle.fasta"
TREEPATH ="../data/BSHtree/clustered/2021-07-15_BSH_filterd_BEAST_umapemm.nwk"
N=417
FAM = 'BSH'


# =============================================================================
# ALNPATH = "../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
# TREEPATH ="../data/OMPtree/2021-07-15_OMPuniprot_muscle_umapemm.tree"
# N=1230
# FAM = 'OMP'
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00005_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/PF00005_seed_umapemm.tree"
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00271_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-15_PF00271_seed_umapemm.tree"
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# =============================================================================
# =============================================================================
# # # ALNPATH = "../data/pfamaln/PF00067_seed.txt"
# # # TREEPATH ="../data/pfamtrees/clustered/2021-07-15_PF00067_seed_umapemm.tree"
# # # TREEID='PF00067'
# # # FAM = 'pfam'
# # # HEIGHT = 22
# =============================================================================
# =============================================================================
# =============================================================================

ALNPATH = "../data/pfamaln/PF00501_seed.txt"
TREEPATH = "../data/pfamtrees/clustered/2021-07-16_PF00501_seed_umapemm.tree"
TREEID='PF00501'
FAM = 'pfam'
N=835
HEIGHT = 22

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00528_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-16_PF00528_seed_umapemm.tree"
# TREEID='PF00528'
# REP=' '
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00069_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-16_PF00069_seed_umapemm.tree"
# TREEID='PF00069'
# REP=' '
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00072_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-16_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/pfamaln/PF00072_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-16_PF00072_seed_umapemm.tree"
# TREEID='PF00072'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================
# # =============================================================================
# 
# ALNPATH = "../data/pfamaln/PF07690_seed.txt"
# TREEPATH = "../data/pfamtrees/clustered/2021-07-15_PF07690_seed_umapemm.tree"
# TREEID='PF07690'
# FAM = 'pfam'
# HEIGHT = 22
# =============================================================================

# =============================================================================
# ALNPATH = "../data/controlaln/2021-06-22_acyltransferase-filtered-homo-muscle.fasta"
# TREEPATH ="../data/controltrees/clustered/2021-07-15_acyltransferase-filtered-homo-muscle_umapemm.1111.mean.tree"
# TREEID='acyltrans'
# N=1000
# FAM = 'control'
# HEIGHT = 22
# REP='.'
# =============================================================================

# =============================================================================
# ALNPATH = "../data/controlaln/2021-06-22_receptor-filtered-homo-muscle.fasta"
# TREEPATH = "../data/controltrees/clustered/2021-07-15_receptor-filtered-homo-muscle_umapemm.mean.tree"
# TREEID='receptor'
# N=1289
# FAM = 'control'
# HEIGHT = 22
# REP='.'
# =============================================================================

# =============================================================================
# ALNPATH= '../data/controlaln/2021-06-22_uniprot-dehydrogenase-filtered-homo-muscle.fasta'
# TREEPATH = "../data/controltrees/clustered/2021-07-15_uniprot-dehydrogenase-filtered-homo-muscleumapemm.1111.median.tree"
# TREEID='dehydrogenase' 
# REP=' '
# FAM = 'control'
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
    finish_node_fn=None,
    case_sensitive_taxon_labels=False,
    preserve_underscores=False,
    suppress_internal_node_taxa=True,
    suppress_leaf_node_taxa=False,
    terminating_semicolon_required=True,
    ignore_unrecognized_keyword_arguments=False,
    )
treeOR.ladderize()

fileType = 'a'
i = 0
taxons = []
catDict = {}
for taxon in treeOR.taxon_namespace:
    j = 0
    print(taxon)
    print('|'.join(str(taxon)[1:-1].split('|')[:]))
    print()
    taxons.append(str(taxon)[1:-1])
    
    if FAM == 'control':
        catDict['|'.join(str(taxon)[1:-1].split('|')[:-1])] = str(taxon)[1:-1]
    elif TREEID=='PF00528':
        catDict['|'.join(str(taxon)[1:-1].split('|')[:])] = str(taxon)[1:-1]        
    else:
        catDict['|'.join(str(taxon)[1:-1].split('|')[:-1])] = str(taxon)[1:-1]
    
    


#'0': ['A', 'P', 'T'], '1': ['C', 'F', 'W', 'Y'], '2': ['D', 'E'], \
#'3': ['I', 'L', 'M', 'V'], '4': ['H', 'K'], '5': ['R', 'Q'], 
#'6': ['N', 'G', 'S']

#'0': ['A', 'P', 'T'],
dict0 = {}
#'1': ['C', 'F', 'W', 'Y'], 
dict1 = {}
#'2': ['D', 'E'], 
dict2 = {}
#'3': ['I', 'L', 'M', 'V'], 
dict3 = {}
#'4': ['H', 'K'], 
dict4 = {}
#'5': ['R', 'Q'], 
dict5 = {}
#'6': ['N', 'G', 'S']
dict6 = {}

gapsDict = {}
seqCount = -1
taxa = []
alnDict = {}

gapsDict['clust'] = []
for i in range(N+1):
    gapsDict[i]=[0]*(len(taxons))
    dict0[str(i)+'APT']=[0]*(len(taxons))
    dict1[str(i)+'CFWY']=[0]*(len(taxons))
    dict2[str(i)+'DE']   =[0]*(len(taxons))
    dict3[str(i)+'ILMV']=[0]*(len(taxons))
    dict4[str(i)+'HK']=[0]*(len(taxons))
    dict5[str(i)+'RQ']=[0]*(len(taxons))
    dict6[str(i)+'NGS']=[0]*(len(taxons))



if (ALNPATH[-5:] == "fasta")|(ALNPATH[-3:] == "txt"):
        
    with open(ALNPATH,'r') as reader:
        for l in reader:
            l=l.strip()
            print(l)

            if l[0] == '>':
                seqCount += 1
                print(catDict)
                print(l)
                if (FAM =='control') | (TREEID == 'sarscov2sgly') | (TREEID == 'PF00528'):
                    try:
                        taxon = catDict[l.strip()[1:].split(" ")[0].replace('_',REP)]
                    except:
                        taxon = catDict[l.strip()[1:].split(" ")[0].replace('YP',REP)]
                else:
                    print(l.strip()[1:].split(" ")[0])
                    taxon = catDict[l.strip()[1:].split(" ")[0]]    
                gapsDict['clust'].append(taxon.split("|")[-1])
                print(taxon,seqCount)
                taxa.append(taxon)
                
                counter = 0
                line = ''
            else:
                print(l)
                line = l
            
            gapCounter = 0
            
            for char in line:
                if char == '-':
                #if (char == 'H') | (char == 'K'):
                    gapCounter+=1
                    gapsDict[counter][seqCount] = gapCounter
                elif (char == 'A') | (char == 'P') | (char == 'T'):
                    gapCounter = 0
                    dict0[str(counter)+'APT'][seqCount] = 1

                elif (char == 'C') | (char == 'F') | (char == 'W') | (char == 'Y'):
                    gapCounter = 0
                    dict1[str(counter)+'CFWY'][seqCount] = 1

                elif (char == 'D') | (char == 'E'):
                    gapCounter = 0
                    dict2[str(counter)+'DE'][seqCount] = 1

                elif (char == 'I') | (char == 'L') | (char == 'V') | (char == 'M'):
                    gapCounter = 0
                    dict3[str(counter)+'ILMV'][seqCount] = 1

                elif (char == 'H') | (char == 'K'):
                    gapCounter = 0
                    dict4[str(counter)+'HK'][seqCount] = 1

                elif (char == 'R') | (char == 'Q'):
                    gapCounter = 0
                    dict5[str(counter)+'RQ'][seqCount] = 1
                
                elif (char == 'N') | (char == 'G') | (char == 'S'):
                    gapCounter = 0
                    dict6[str(counter)+'NGS'][seqCount] = 1

                else:
                    gapCounter = 0
                    
                counter += 1

            
                
d = {}
d['clust'] = gapsDict['clust']
for k in gapsDict.keys():
    if k != 'clust':
        if sum(gapsDict[k]) > 3 :
            print(k,len(gapsDict[k]))
            medianGap = statistics.median(gapsDict[k])
            d[str(k)] = gapsDict[k]

            for i in range(len(gapsDict[k])):
               d[str(k)][i] = int(d[str(k)][i] > 0)#medianGap)
            #d[str(k)] = gapsDict[k]
            if len(gapsDict[k]) != len(taxa):
                print('woah')
                print(k,len(gapsDict[k]))

d0={}
for k in dict0.keys():
    if k != 'clust':
        if sum(dict0[k]) > 3:
            d0[str(k)] = dict0[k]

d1={}
for k in dict1.keys():
    if k != 'clust':
        if sum(dict1[k]) > 3:
            d1[str(k)] = dict1[k]
d2={}
for k in dict2.keys():
    if k != 'clust':
        if sum(dict2[k]) > 3:
            d2[str(k)] = dict2[k]

d3={}
for k in dict3.keys():
    if k != 'clust':
        if sum(dict3[k]) > 3:
            d3[str(k)] = dict3[k]
   
d4={}
for k in dict4.keys():
    if k != 'clust':
        if sum(dict4[k]) > 3:
            d4[str(k)] = dict4[k]
   
d5={}
for k in dict5.keys():
    if k != 'clust':
        if sum(dict5[k]) > 3:
            d5[str(k)] = dict5[k]

d6={}
for k in dict6.keys():
    if k != 'clust':
        if sum(dict6[k]) > 3:
            d6[str(k)] = dict6[k]

print(len(taxa),taxa)
df = pd.DataFrame(data=d,index=taxa)
df0 = pd.DataFrame(data=d0,index=taxa)
df1 = pd.DataFrame(data=d1,index=taxa)
df2 = pd.DataFrame(data=d2,index=taxa)
df3 = pd.DataFrame(data=d3,index=taxa)
df4 = pd.DataFrame(data=d4,index=taxa)
df5 = pd.DataFrame(data=d5,index=taxa)
df6 = pd.DataFrame(data=d6,index=taxa)
df = df.merge(df0,left_index=True,right_index=True)
df = df.merge(df1,left_index=True,right_index=True)
df = df.merge(df2,left_index=True,right_index=True)
df = df.merge(df3,left_index=True,right_index=True)
df = df.merge(df4,left_index=True,right_index=True)
df = df.merge(df5,left_index=True,right_index=True)
df = df.merge(df6,left_index=True,right_index=True)

print(df.index)
print(df.columns)



# Split the Groups from the dataset where y is category and x is data with species
y = df.iloc[:,0]
x = df.iloc[:,1:]

print(y)
print(x)


# Split the data into training and test data for the categories(y) and dataset(x)
# Here we are spliting it 65% training and 35% test
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=.35, random_state=42)


ensemble_clfs = [
    ("RandomForestClassifier, max_features='sqrt'",
        RandomForestClassifier(warm_start=True, 
                               oob_score=True,
                               max_features="sqrt", 
                               random_state=42)),
    ("RandomForestClassifier, max_features='log2'",
        RandomForestClassifier(warm_start=True, 
                               max_features='log2',
                               oob_score=True, 
                               random_state=42)),
    ("RandomForestClassifier, max_features=None",
        RandomForestClassifier(warm_start=True, 
                               max_features=None,
                               oob_score=True, 
                               random_state=42))
]

# Map a classifier name to a list of (<n_estimators>, <error rate>) pairs.
error_rate = OrderedDict((label, []) for label, _ in ensemble_clfs)

# Range of `n_estimators` values to explore.
min_estimators = 50
max_estimators = 100

for label, clf in ensemble_clfs:
    for i in range(min_estimators, max_estimators + 1):
        clf.set_params(n_estimators=i)
        clf.fit(X_train, y_train)

        # Record the OOB error for each `n_estimators=i` setting.
        oob_error = 1 - clf.oob_score_
        error_rate[label].append((i, oob_error))

# Generate the "OOB error rate" vs. "n_estimators" plot.
for label, clf_err in error_rate.items():
    xs, ys = zip(*clf_err)
    plt.plot(xs, ys, label=label)

plt.xlim(min_estimators, max_estimators)
plt.xlabel("n_estimators")
plt.ylabel("OOB error rate")
plt.legend(loc="upper right")
plt.show()

clf = RandomForestClassifier(n_estimators=70, max_features=None, random_state=42)
all_accuracies = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=2)
print('Mean Validation Scores: ' ,end='')
print(np.mean(all_accuracies))

clf_final = RandomForestClassifier(n_estimators=70, bootstrap=True,max_features=None,oob_score= True,
                                   random_state= 42)
clf_final.fit(X_train,y_train)
y_pred = clf_final.predict(X_test)
print("Test Set Accuracy:",metrics.accuracy_score(y_test, y_pred))


print(clf_final.oob_score_)

# Visualise classical Confusion M0atrix
from sklearn.metrics import confusion_matrix
CM = confusion_matrix(y_test, y_pred)
print(CM)

# Visualize it as a heatmap
import seaborn
seaborn.heatmap(CM, cmap="YlGnBu")
plt.show()

feats = {} # a dict to hold feature_name: feature_importance
for feature, importance in zip(df.columns, clf_final.feature_importances_):
    feats[feature] = importance #add the name/value pair 

importances = pd.DataFrame.from_dict(feats, orient='index').rename(columns={0: 'Gini-importance'})
imp=importances.sort_values(by='Gini-importance',ascending=False)
print(imp.head(20))
print(list(imp.index)[:20])
alnImp = list(imp.index)[:20]
imp.head(20).to_csv("../data/RF/2021-07-15_RF_"+TREEID+FAM+".csv")




alnDict = {}

#gapsDict['clust'] = []
for i in range(len(taxons)):
    alnDict[i]=[0]*N


with open(ALNPATH,'r') as reader:
    for l in reader:
       #print(l.strip())
        line = l.strip()
        if len(line)==0:
            pass
        elif line[0] == '#' or line[0] == '>':
            taxa = line[1:].replace('_','.')
            alnDict[taxa] = ''
        else:
            alnDict[taxa] += line

c=0
for impo in alnImp:
    alnImp[c] = str(impo)
    c+=1 

alnImp.sort()
l=[]
lDict={}
for taxon, aln in alnDict.items():
    if type(taxon) == int or taxon == 'mega':
        pass
    else:
        taxon = taxon.strip().split(' ')[0]      
        toAppend=[]        
        
        if TREEID == 'dehydrogenase':
            toAppend.append('>'+catDict[taxon[:].replace('.',' ')])
        elif TREEID == 'sarscov2sgly':
            toAppend.append('>'+catDict[taxon[:].replace('YP','YP ')])
        else:
            toAppend.append('>'+catDict[taxon[:].replace('..','')])
        for feat in alnImp:
            number = ''
            for word in feat:
                if word.isdigit():
                    number = number + word
            number = str(number)
            
            
            toAppend.append(aln[int(number)-5:int(number)+4])
        l.append(toAppend)

        #fig gini imp 
        #string order by cluster number


table = PrettyTable(['Name', alnImp[0],alnImp[1],alnImp[2], alnImp[3], alnImp[4],
                     alnImp[5], alnImp[6], alnImp[7], alnImp[8], alnImp[9],
                     alnImp[10], alnImp[11], alnImp[12], alnImp[13], alnImp[14],
                     alnImp[15], alnImp[16], alnImp[17], alnImp[18], alnImp[19]
                     ])

for rec in l:
    table.add_row(rec)
    
OUTF="../data/"+FAM+"colored/2021-07-15"+FAM+TREEID+"_RF.csv"
print(OUTF)

import csv
with open(OUTF,'w') as f :
    writer=csv.writer(f)
    for k in l:
        writer.writerow(k)

counter = 0
print(alnImp)
for i in alnImp:
    print(i)