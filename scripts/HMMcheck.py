# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 15:48:07 2021

@author: aorti
"""

import re
import string
import sys
import os
from optparse import OptionParser
from SBL import HMMER_Wrapper
from HMMER_Wrapper import *

name_line_re = re.compile(">([A-Za-z0-9\.-]+)")
fasta_line_re = re.compile("^[A-Za-z]+")
parser = OptionParser()
parser.add_option(
    "-f",
    "--file",
    dest="file_name",
    type="string",
    help="Input a file containing a list of fasta sequences.",
)
parser.add_option(
    "-d",
    "--database",
    dest="database",
    type="string",
    help="Input a FASTA database to perform the query. ",
)
(options, args) = parser.parse_args()
if not options.file_name:
    sys.exit("You must provide a fasta file. ")
elif not options.database:
    sys.exit("You must provide a database. ")
else:
    fasta_sequences = dict()
    fasta_file = open(options.file_name)
    name = ""
    for line in fasta_file:
        name_line = name_line_re.search(line)
        fasta_line = fasta_line_re.search(line)
        if name_line:
            name = name_line.group(1)
        if fasta_line:
            if not name:
                sys.exit("Badly formatted fasta file, each sequence requires a name...")
            fasta_sequences[name] = line
            name = ""

    hmm = HMM("ex_hmm", fasta_sequences, options.database)

    hmm.run()

    for accession in hmm.get_results(0.1):
        print(accession)
    fasta_file.close()
