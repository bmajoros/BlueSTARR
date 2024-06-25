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

RANGE=0.99

def getRange(pct,sortedThetas):
    M=len(sortedThetas)
    bothTails=(1-pct)*M
    eachTail=int(bothTails/2)
    lower=sortedThetas[eachTail]
    upper=sortedThetas[M-1-eachTail]
    return (lower,upper)

def computeDistribution(thetas,numBins):
    sortedThetas=sorted(thetas)
    M=len(sortedThetas)
    (lower,upper)=getRange(RANGE,sortedThetas)
    #print("range=",lower,upper)
    thetaRange=upper-lower
    binSize=thetaRange/numBins
    distr=[]
    pos=0
    for i in range(1,numBins+1):
        binEnd=lower+i*binSize
        if(i==numBins): binEnd=sortedThetas[M-1]+1
        (count,pos)=countInBin(pos,sortedThetas,binEnd)
        p=float(count)/float(M)
        distr.append([p,binEnd])
    #print("distr=",distr)
    return distr

def countInBin(pos,sortedThetas,binEnd):
    M=len(sortedThetas)
    firstPos=pos
    while(pos<M and sortedThetas[pos]<binEnd): pos+=1
    count=pos-firstPos
    return (count,pos)

def loadThetas(thetasFile):
    IN=open(thetasFile,"rt")
    thetas=[]
    for line in IN:
        if(not rex.find("\d",line)): continue
        theta=float(line)
        thetas.append(math.log(theta))
    IN.close()
    return thetas

def computeAccept(distr,N,M,B):
    alpha=[]
    for i in range(B):
        (pi,binEnd)=distr[i]
        a=min(1,float(N)/(float(M)*float(B)*pi))
        alpha.append(a)
    return alpha

def nextCountsLine(line,DNA_REPS,RNA_REPS):
    fields=line.rstrip().split()
    fields=[int(x) for x in fields]
    DNA=fields[:DNA_REPS]; RNA=fields[DNA_REPS:]
    sumDNA=sum(DNA);  sumRNA=sum(RNA)
    naiveTheta=(sumRNA+RNA_REPS)/(sumDNA+DNA_REPS)
    return math.log(naiveTheta)

def getBin(distr,theta):
    B=len(distr)
    for i in range(B):
        (prop,end)=distr[i]
        if(theta<end): return i
    return B-1
    #raise Exception("No bin found for theta=",theta)

def countCounts(countsFile):
    IN_COUNTS=gzip.open(countsFile,"rt")
    header=IN_COUNTS.readline()
    M=0
    for line in IN_COUNTS: M+=1
    return M
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <in:naive-thetas.txt> <in:fasta.gz> <in:counts.txt.gz> <#bins> <desired-sample-size> <out-dir> <train|test|validation>")
(thetasFile,fastaFile,countsFile,numBins,N,outDir,fileLabel)=sys.argv[1:]
numBins=int(numBins)
N=int(N)
outCountFile=outDir+"/"+fileLabel+"-counts.txt.gz"
outFastaFile=outDir+"/"+fileLabel+".fasta.gz"

# Compute empirical distribution of theta
thetas=loadThetas(thetasFile)
distr=computeDistribution(thetas,numBins)

# Compute acceptance probability for each bin
M=countCounts(countsFile)
alpha=computeAccept(distr,N,M,numBins)

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
    theta=nextCountsLine(line,DNA_REPS,RNA_REPS)
    pair=fastaReader.nextSequence() # returns None at eof
    if(pair is None): raise Exception("Reached eof of FASTA file prematurely")
    (defline,sequence)=pair
    b=getBin(distr,theta)
    accept=alpha[b]
    if(random.random()<=accept):
        print(line,end="",file=OUT_COUNTS)
        fastaWriter.addToFasta(defline,sequence,OUT_FASTA)
IN_COUNTS.close(); OUT_COUNTS.close(); OUT_FASTA.close()

