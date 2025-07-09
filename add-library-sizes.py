#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import gzip

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <raw-input.gz> <lib-sizes.txt> <output.txt.gz>\n")
(rawFile,libSizesFile,outFile)=sys.argv[1:]

# Read library sizes
sizes=[]; labels=[]
with open(libSizesFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        labels.append(fields[0]+"_libsize")
        sizes.append(fields[1])
sizes="\t".join(sizes)
labels="\t".join(labels)

# Process the raw data file and append library sizes
OUT=gzip.open(outFile,"wt")
with gzip.open(rawFile,"rt") as IN:
    header=IN.readline()
    header=header.rstrip()+"\t"+labels
    print(header,file=OUT)
    for line in IN:
        line.rstrip()
        print(line,sizes,sep="\t",file=OUT)
OUT.close()


