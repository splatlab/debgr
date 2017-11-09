# debgr
deBGR: An Efficient and Near-Exact Representation of the Weighted de Bruijn Graph

This work appeared at ISMB 2017. If you use this software please cite us.
```
@article{doi:10.1093/bioinformatics/btx261,
author = {Pandey, Prashant and Bender, Michael A. and Johnson, Rob and Patro, Rob},
title = {deBGR: an efficient and near-exact representation of the weighted de Bruijn graph},
journal = {Bioinformatics},
volume = {33},
number = {14},
pages = {i133-i141},
year = {2017},
doi = {10.1093/bioinformatics/btx261},
URL = { + http://dx.doi.org/10.1093/bioinformatics/btx261},
eprint = {/oup/backfile/content_public/journal/bioinformatics/33/14/10.1093_bioinformatics_btx261/2/btx261.pdf}
}
```
Overview
--------

deBGR a memory-efficient and near-exact representation of the weighted de
Bruijn Graph. Our representation is based upon a recently-introduced counting
filter data structure which, itself, provides an approximate representation of
the weighted de Bruijn Graph.

deBGR representation yields a four-orders-of-magnitude reduction in the number
of errors made by a state-of-the-art (Squeakr) approximate weighted de Bruijn
Graph representation, while increasing the space requirements by less than
10%. Our technique is based on a simple invariant that all weighted de Bruijn
Graphs must satisfy, and hence is likely to be of general interest and
applicable in most weighted de Bruijn Graph-based systems.

API
--------
* 'main': count k-mers in a read dataset and build a k-mer representation and auxillary data structures
* 'debruijn_graph': correct k-mer abundances in the representation

Build
-------
This library depends on libssl and boost. To build a compact and near-exact dBG representation you first need to run "main" and then "debruijn_graph" on the read dataset. The CQF code uses two new instructions to implement select on machine words introduced in intel's Haswell line of CPUs. However, there is also an alternate implementation of select on machine words to work on CPUs older than Haswell. To build on an older hardare (older than Haswell) use "NH=1" as a make argument.

```bash
 $ make main
 $ ./main 0 20 12 1 test.fastq
```

 Following are the argumenrs to main:
 - file format: 0 - plain fastq, 1 - gzip compressed fastq, 2 - bzip2 compressed fastq
 - CQF size: the log of the number of slots in the CQF
 - exact CQF size: the log of the number of slots in exact CQF to count start and end occurrences of k-mers
 - num of threads: number of threads to count
 - file(s): "filename" or "dirname/*" for all the files in a directory

```bash
 $ make debruijn_graph
 $ ./main test.fastq
```

 Following are the argumenrs to debruijn_graph:
 - file: "filename" if only one file is used with "main" or name of the first file in the directory if a directory is used with "main"

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
