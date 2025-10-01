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
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <in-counts.txt.gz> <in-libs.txt> <out-counts.txt.gz\n")
(inCountsFile,inLibFile,outCountsFile)=sys.argv[1:]

newLibs=[]
with open(inLibFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        newLibs.append(fields[1])
numLibSizes=len(newLibs)

OUT=gzip.open(outCountsFile,"wt")
with gzip.open(inCountsFile,"rt") as IN:
    header=IN.readline()
    print(header,file=OUT,end="")
    for line in IN:
        fields=line.rstrip().split()
        numFields=len(fields)
        numKeep=numFields-numLibSizes
        newFields=fields[:numKeep]
        newFields.extend(newLibs)
        print("\t".join(newFields),file=OUT)
OUT.close()



