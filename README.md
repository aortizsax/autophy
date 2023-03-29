# Autophy

AutoPhy is a P


The primary home page for 

## Requirements and Installation

AutoPhy runs under Python 3 ( > 3.7, <=3.9)

You can install AutoPhy by running the following commands. To clone the repo:
```bash
git clone https://github.com/aortizsax/autophy.git
```

Change directory in to cloned repo:

```bash
cd autophy
```

In the future, you can update AutoPhy by running the following commands after changing directory into repo folder. To update the repo run this:
```bash
git fetch
```

Create the conda enviroment:

```bash
conda create -n autophy python==3.8
```

Activate conda enviroment:

```bash
conda activate autophy
```

Install repo as python package:

```bash
pip install .
```

Check build

```bash
autophy -h
```

## Documentation

### Introduction
In the following section we will walk the user through the analysis using the OMP data outlined in the paper (DOI: ), as well as detailing where input your specific data  

### Your external data
For this analysis the user will need their alignment file in fasta format, tree in newick format, and two string identifiers (Family name or acronym; orinating database)


### Usage and Analysis

Usage:
```bash
autophy -h
```
returns:
```{bash}
usage: autophy [-h] [-a ALIGNMENT] [-t TREE] [-id TREEID] [-d DATABASE]
               [-H HEIGHT] [-o OUTPUT_PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -a ALIGNMENT, --alignment ALIGNMENT
                        Suffix for output files [default=../data/OMPaln/2021-0
                        6-11_OMPuniprot_muscle.fasta].
  -t TREE, --tree TREE  Suffix for output files [default=../data/uniprotbeast/
                        uniprotOMP_muscle_1647964589147_cah.tree].
  -id TREEID, --treeid TREEID
                        Suffix for output files [default=OMP].
  -d DATABASE, --database DATABASE
                        Suffix for output files [default=uniprot].
  -H HEIGHT, --height HEIGHT
                        Suffix for output files [default=22].
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Suffix for output files [default=umapgmmemm].
```

To run:
```{bash}
cd data
autophy -a ../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta -t ../data/uniprotbeast/uniprotOMP_muscle_1647964589147_cah.tree -id OMP -d uniprot
```
Or for your data
```{bash}
autophy -a ALIGNMENTFILE.fasta -t TREE.tree -id ID -d DATASET 
```

First UMAP projection. Parameters used 

![](./data/output/2023-03-28_10_ogUMAPproguniprotOMP_precomputed.svg)

Second, we fit Gaussian Mixture Model(GMM) using Expectation Maximization(EM) to the UMAP projection. GMM fits Gaussian peaks to any dimension. We score GMM models of increasing peak number using BIC score and pick the model with the lowest BIC score.

![](./data/output/2023-03-28_10_GMMsweep_uniprotOMP_precomputed.svg)


Third, in the case that when GMM peaks in the model pertain to paraphyletic clusters, we slpit those GMM peaks into simple monoplytic clusters and add a decimal at the end of the label so the user knows which GMM peak they originated from. Finally, the program prints the final clustered and colored tree 

![](./data/output/2023-03-28_10_EMClust_uniprotOMP_precomputed_coloredtree.svg)


## License and Warranty
Please see the file "LICENSE" for details.
