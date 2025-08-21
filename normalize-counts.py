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
import os
import ProgramName
import gzip
from Rex import Rex
rex=Rex()
BIG_NUMBER=1e7

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <infile> <outdir>\n")
(infile,outdir)=sys.argv[1:]

outfile=outdir+"/"+infile
OUT=gzip.open(outfile,"wt")
IN=gzip.open(infile,"rt")
header=IN.readline()
print(header,end="",file=OUT)
rex.findOrDie("DNA=(\\d+)\\s+RNA=(\\d+)",header)
numDNA=int(rex[1]); numRNA=int(rex[2])
for line in IN:
    fields=line.rstrip().split()
    fields=[int(x) for x in fields]
    DNAcounts=fields[:numDNA]
    RNAcounts=fields[numDNA:(numDNA+numRNA)]
    DNAlibs=fields[(numDNA+numRNA):(2*numDNA+numRNA)]
    RNAlibs=fields[(2*numDNA+numRNA):]
    for i in range(numDNA):
        DNAcounts[i]=DNAcounts[i]/DNAlibs[i]*BIG_NUMBER
    for i in range(numRNA):
        RNAcounts[i]=RNAcounts[i]/RNAlibs[i]*BIG_NUMBER
    fields=DNAcounts; fields.extend(RNAcounts)
    print("\t".join([str(x) for x in fields]),file=OUT)
IN.close(); OUT.close()



