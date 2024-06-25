#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2023 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import gzip
import math
import random
from FastaReader import FastaReader
from FastaWriter import FastaWriter
import ProgramName
from Rex import Rex
rex=Rex()

def computeAccept(N,M):
    alpha=float(N)/float(M)
    return alpha

def countCounts(countsFile):
    IN_COUNTS=gzip.open(countsFile,"rt")
    header=IN_COUNTS.readline()
    M=0
    for line in IN_COUNTS: M+=1
    return M
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <in:fasta.gz> <in:counts.txt.gz> <desired-sample-size> <out-dir> <train|test|validation>")
(fastaFile,countsFile,N,outDir,fileLabel)=sys.argv[1:]
N=int(N)
outCountFile=outDir+"/"+fileLabel+"-counts.txt.gz"
outFastaFile=outDir+"/"+fileLabel+".fasta.gz"

# Compute acceptance probability for each bin
M=countCounts(countsFile)
alpha=computeAccept(N,M)

# Create output files
OUT_COUNTS=gzip.open(outCountFile,"wt")
OUT_FASTA=gzip.open(outFastaFile,"wt")

# Downsample the data
IN_COUNTS=gzip.open(countsFile,"rt")
header=IN_COUNTS.readline()
if(not rex.find("DNA=(\d+)\s+RNA=(\d+)",header)):
    raise Exception("Can't parse header in counts file: "+header)
print(header,end="",file=OUT_COUNTS)
DNA_REPS=int(rex[1]); RNA_REPS=int(rex[2])
fastaReader=FastaReader(fastaFile)
fastaWriter=FastaWriter()
for line in IN_COUNTS:
    pair=fastaReader.nextSequence() # returns None at eof
    if(pair is None): raise Exception("Reached eof of FASTA file prematurely")
    (defline,sequence)=pair
    if(random.random()<=alpha):
        print(line,end="",file=OUT_COUNTS)
        fastaWriter.addToFasta(defline,sequence,OUT_FASTA)
IN_COUNTS.close(); OUT_COUNTS.close(); OUT_FASTA.close()

