#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
## Copyright (c) 2023 Adrian Ortiz.
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
##############################################################################

import os
import pathlib
import sys
import argparse
from autophy.analyze.analysis import cluster

def main():
    parser = argparse.ArgumentParser(description=None)
#    parser.add_argument(
 #           "src_paths",
  #          action="store",
   #         nargs="+",
    #        metavar="FILE",
     #       help="Path to source file(s).")
    parser.add_argument(
            "-a", "--alignment",
            action="store",
            default="../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta",
            help="Suffix for output files [default=%(default)s].")        
    parser.add_argument(
            "-t", "--tree",
            action="store",
            default="../data/uniprotbeast/uniprotOMP_muscle_1647964589147_cah.tree",
            help="Suffix for output files [default=%(default)s].")
    parser.add_argument(
            "-id", "--treeid",
            action="store",
            default="OMP",
            help="Suffix for output files [default=%(default)s].")
    parser.add_argument(
            "-d", "--database",
            action="store",
            default="OMP",
            help="Suffix for output files [default=%(default)s].")
    parser.add_argument(
            "-H", "--height",
            action="store",
            default="22",
            help="Suffix for output files [default=%(default)s].")
    parser.add_argument(
            "-o", "--output-prefix",
            action="store",
            default="umapgmmemm",
            help="Suffix for output files [default=%(default)s].")
    args = parser.parse_args()
    print("Hello, world.")
    print(dir(args))
    
#    'alignment', 'database', 'height', 'output_prefix', 'tree', 'treeid'
    
    alignment = args.alignment
    database = args.database
    height = args.height
    output_prefix = args.output_prefix
    treefile = args.tree
    treeid = args.treeid
     
    cluster(alignment, database, height, output_prefix, treefile, treeid)
    #ALNPATH = "../data/OMPaln/2021-06-11_OMPuniprot_muscle.fasta"
    #TREEPATH = "../data/uniprotbeast/uniprotOMP_muscle_1647964589147_cah.tree"
    #OUTTREE = "../data/uniprottrees/uniprotOMP_muscle_cah_umapemm.tree"
    #Naln=1230
    #TREEID='OMP'#'uniprot'
    #FAM = 'uniprot'#'OMP'
    #HEIGHT = 22


if __name__ == '__main__':
    main()
