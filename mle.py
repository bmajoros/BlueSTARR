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
import math
import numpy as np
from scipy.optimize import minimize
import ProgramName
import gzip
from Rex import Rex
rex=Rex()
global numDna
global numRna

def log(x):
    return math.log(x)

def lgamma(x):
    return math.lgamma(x)

def logLik(sumX,numX,rna,theta,alpha,beta,sumDnaLibs,RnaLibs):
    total=0
    for i in range(numRna):
        Yj=rna[i]
        libRatio=RnaLibs[i]/sumDnaLibs
        thetaL=theta*libRatio
        LL=(sumX+alpha)*log(beta+numX)+lgamma(Yj+sumX+alpha)+Yj*log(thetaL)\
            -lgamma(sumX+alpha)-lgamma(Yj+1)-\
            (Yj+sumX+alpha)*log(thetaL+beta+numX)
        total+=LL
        #print(LL)
    return total

def getClosure(data):
    def f(theta):
        theta=theta[0]
        totalLogLik=0
        for case in data:
            totalLogLik+=logLik(case.dna,numDna,case.rna,theta,1e-10,
                                1e-10,case.dnaLibSum,case.rnaLibs)
        return -totalLogLik
    return f

class Case:
    def __init__(self,dna,rna,dnaLibSum,rnaLibs):
        self.dna=dna
        self.rna=rna
        self.dnaLibSum=dnaLibSum
        self.rnaLibs=rnaLibs

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <in-counts.txt.gz> <max-data>\n")
(inCountsFile,maxData)=sys.argv[1:]
maxData=int(maxData)

# Read data
data=[]
lineNum=1
with gzip.open(inCountsFile,"rt") as IN:
    header=IN.readline()
    if(not rex.find("DNA=(\\d+)\\s+RNA=(\\d+)",header)):
        raise Exception("expecting RNA=# DNA=# header")
    numDna=int(rex[1]); numRna=int(rex[2])
    firstLib=numDna+numRna
    for line in IN:
        fields=line.rstrip().split()
        fields=[int(x) for x in fields]
        dnaSum=float(sum(fields[:numDna]))
        dnaLibSum=sum(fields[firstLib:(firstLib+numDna)])
        rnaLibs=fields[(firstLib+numDna):]
        rna=fields[numDna:(numDna+numRna)]
        #print("rna=",rna,"dnaSum=",dnaSum,"rnaLibs=",rnaLibs,"dnaLibSum=",dnaLibSum)
        #naive=(rna/rnaLibs[i])/(dnaSum/dnaLibSum)
        case=Case(dnaSum,rna,dnaLibSum,rnaLibs)
        data.append(case)
        lineNum+=1
        if(lineNum>=maxData): break

# Minimize the negative log likelihood
theta0 = np.array([0.5])
opt = minimize(getClosure(data), theta0, method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True})
print("theta =",opt.x[0])


