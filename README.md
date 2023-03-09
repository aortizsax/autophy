AutoPhy is a Python library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and
manipulation of phylogenetic trees and character matrices, and supports the
reading and writing of phylogenetic data in a range of formats, such as NEXUS,
NEWICK, NeXML, Phylip, FASTA, etc.  Application scripts for performing some
useful phylogenetic operations, such as data conversion and tree posterior
distribution summarization, are also distributed and installed as part of the
libary.  DendroPy can thus function as a stand-alone library for phylogenetics,
a component of more complex multi-library phyloinformatic pipelines, or as a
scripting "glue" that assembles and drives such pipelines.

The primary home page for 

## Requirements and Installation
=============================

AutoPhy runs under Python 3 ( > 3.7, <=3.10)

You can install AutoPhy by running::

```bash
$ gh repo clone aortizsax/autophy
```

```bash
$ cd autophy
```

```bash
$ conda create -n autophy
```

```bash
$ conda activate autophy
```

```bash
$ pip install .
```

```bash
$ autophy -h
```

## Documentation
=============
OLD
Full documentation is available here:

    http://dendropy.org/

This includes:

    -   `A comprehensive "getting started" primer <http://dendropy.org/primer/index.html>`_ .
    -   `API documentation <http://dendropy.org/library/index.html>`_ .
    -   `Descriptions of data formats supported for reading/writing <http://dendropy.org/schemas/index.html>`_ .

and more.

License and Warranty
====================

Please see the file "LICENSE" for details.
