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
from FastaReader import FastaReader
from Rex import Rex
rex=Rex()

def getKeeperList(inFastaFile):
    keep=set()
    with gzip.open(inFastaFile,"rt") as IN:
        for line in IN:
            if(rex.find("> /coord=(\S+:\d+-\d+)",line)):
                keep.insert(rex[1])
    return keep

def process(filestem,inFile,inDir,outDir,dnaReps,rnaReps):
    inFastaFile=inDir+"/"+filestem+".fasta.gz"
    #inCountsFile=inDir+"/"+filestem+"-counts.txt.gz"
    outFastaFile=outDir+"/"+filestem+".fasta.gz"
    outCountsFile=outDir+"/"+filestem+"-counts.txt.g"
    keep=getKeeperList(inFastaFile)
    COUNTS=gzip.open(outCountsFile,"wt")
    FASTA=gzip.open(outFastaFile,"wt")
    print("DNA=",dnaReps," RNA=",rnaReps,sep="",file=COUNTS)
    writeRecords(inFile,keep,COUNTS,FASTA)
    COUNTS.close(); FASTA.close()

def writeRecords(inFile,keep,COUNTS,FASTA):
    headerFields=None
    with gzip.open(inFile,"rt") as IN:
        header=IN.readline()
        headerFields=header.rstrip().split()
        for line in IN:
            fields=line.rstrip().split()
            key=fields[0]+":"+fields[1]+"-"+fields[2]
            if key not in keep: continue
            endIdx=headerFields.index("end")
            fcIdx=headerFields.index("log2FC")
            seqIdx=headerFields.index("sequence")
            counts=fields[(endIdx+1):fcIdx]
            libsizes=fields[(seqIdx+1):]
            counts.extend(libsizes)
            print("\t".join(counts)) #,file=COUNTS)


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <#DNA-reps> <#RNA-reps> <all-data.txt.gz> <in-dir> <out-dir>\n")
(dnaReps,rnaReps,inFile,inDir,outDir)=sys.argv[1:]

if(inDir==outDir): exit("Input and output directories cannot be the same")

process("train",inFile,inDir,outDir,dnaReps,rnaReps)
#process("validation",inFile,inDir,outDir,dnaReps,rnaReps)
#process("test",inFile,inDir,outDir,dnaReps,rnaReps)



