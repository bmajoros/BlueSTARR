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
import numpy as np
import gzip
from Rex import Rex
import TempFilename
rex=Rex()
rng = np.random.default_rng()
alpha=1e-5
beta=1e-5
#scalingFactor=1.0/7.0
tempfile=TempFilename.generate()

def downsample(X):
    shape=alpha+X
    rate=beta+1
    scale=1/rate
    g=rng.gamma(shape, scale * scalingFactor)
    newX=rng.poisson(g)
    return newX

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <IN:counts.txt.gz> <outdir> <scaling-factor>\n")
(infile,outdir,scalingFactor)=sys.argv[1:]
scalingFactor=float(scalingFactor)

# Downsample into a temporary file without library sizes
if(outdir=="."): raise Exception("operation will overwite input file")
OUT=gzip.open(tempfile,"wt")
IN=gzip.open(infile,"rt")
header=IN.readline()
print(header,end="",file=OUT)
rex.findOrDie("DNA=(\\d+)\\s+RNA=(\\d+)",header)
numDNA=int(rex[1]); numRNA=int(rex[2])
newLibSizes=[0]*(numDNA+numRNA)
for line in IN:
    fields=line.rstrip().split()
    fields=[int(x) for x in fields]
    DNAcounts=fields[:numDNA]
    RNAcounts=fields[numDNA:(numDNA+numRNA)]
    DNAlibs=fields[(numDNA+numRNA):(2*numDNA+numRNA)]
    RNAlibs=fields[(2*numDNA+numRNA):]
    newX=[downsample(x) for x in DNAcounts]
    newY=[downsample(y) for y in RNAcounts]
    fields=newX; fields.extend(newY)
    for i in range(numDNA+numRNA):
        newLibSizes[i]+=fields[i]
    print("\t".join([str(x) for x in fields]),file=OUT)
IN.close(); OUT.close()

# Add library sizes and move into target output file
outfile=outdir+"/"+infile
OUT=gzip.open(outfile,"wt")
IN=gzip.open(tempfile,"rt")
header=IN.readline()
print(header,end="",file=OUT)
for line in IN:
    fields=line.rstrip().split()
    fields=[int(x) for x in fields]
    DNAcounts=fields[:numDNA]
    RNAcounts=fields[numDNA:(numDNA+numRNA)]
    fields.extend(newLibSizes)
    print("\t".join([str(x) for x in fields]),file=OUT)
IN.close(); OUT.close()
os.remove(tempfile)
