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
import random
import ProgramName
import TempFilename
from FastaReader import FastaReader
from Pipe import Pipe
from Rex import Rex
rex=Rex()

# THESE SHOULD BE MOVED INTO A CONFIGURATION FILE:
TWO_BIT_TO_FA="/hpc/home/bmajoros/twobit/twoBitToFa"
TWO_BIT="/datacommons/allenlab/hg38/hg38.2bit"
CHILD="/hpc/group/majoroslab/deepstarr/git/mutator-child.py"
MAX_N=-1 # -1 = no maximum

# GLOBALS
ALPHA="ACGT"

# CLASSES
class CRE:
    def __init__(self,ID):
        self.ID=ID
        if(not rex.find("(\S+):(\d+)-(\d+)",ID)):
            raise Exception("Can't parse CRE: "+ID)
        self.chrom=rex[1]; self.begin=int(rex[2]); self.end=int(rex[3])
        self.positions=[]
        self.window=None
    def addPos(self,pos):
        self.positions.append(pos)

class CrePosition:
    def __init__(self,text):
        if(not rex.find("(\d+):ref=(.):(.),(.)",text)):
            raise Exception("Can't parse pos: "+text)
        self.pos=int(rex[1])
        self.ref=rex[2];
        self.alleles=[rex[3],rex[4]]

# FUNCTIONS
def loadCREs(filename,MAX_N):
    n=0
    recs=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)<2): continue
            rec=CRE(fields[0])
            for site in fields[1:]:
                rec.addPos(CrePosition(site))
            recs.append(rec)
            n+=1
            if(MAX_N>0 and n>=MAX_N): break
    return recs

def writeCoords(recs,seqLen,filename):
    halfLen=int(seqLen/2)
    COORD=open(filename,"wt")
    for i in range(len(recs)):
        rec=recs[i]
        chrom=rec.chrom; oldBegin=rec.begin; oldEnd=rec.end
        center=int((oldBegin+oldEnd)/2)
        begin=center-halfLen
        end=begin+seqLen
        print(chrom+":"+str(begin)+"-"+str(end),file=COORD)
        localBegin=None; localEnd=None
        if(oldBegin<=begin):
            localBegin=0; localEnd=seqLen
        else:
            localBegin=oldBegin-begin; localEnd=oldEnd-begin
        ID=chrom+":"+str(oldBegin)+"-"+str(oldEnd)
        #recs[i]=[chrom,ID,begin,end,localBegin,localEnd]
        rec.window=[begin,end]
    COORD.close()

def makeFasta(coordFile,fastaFile):
    # twoBitToFa needs 0-based half-open coordinates
    cmd=TWO_BIT_TO_FA+" -noMask -seqList="+coordFile+" "+TWO_BIT+" "+fastaFile
    Pipe.run(cmd)

# NOTE: the "ID" is the cCRE interval from ENCODE; whereas the
#       "actualInterval" is the actual interval fed to the model
#       (we grow or shrink to make all sequences the same length)
def runModel(recs,seqLen,fastaFile,inputFile,outputFile,JOBSIZE,model):
    halfLen=int(seqLen/2)
    OUT=open(outputFile,"wt")
    INPUTS=open(inputFile,"wt")
    reader=FastaReader(fastaFile)
    reader.doUppercase()
    creIndex=0
    n=0
    while(True):
        pair=reader.nextSequence()
        if(pair is None): break
        (defline,seq)=pair
        (actualInterval,attr)=FastaReader.parseDefline(defline)
        cre=recs[creIndex]
        chrom=cre.chrom
        (globalBegin,globalEnd)=cre.window
        positions=cre.positions
        for Pos in positions:
            pos=Pos.pos; ref=Pos.ref; alleles=Pos.alleles
            for allele in alleles:
                local=pos-globalBegin
                if(seq[local]!=ref):
                    raise Exception("Ref mismatch:"+seq[local]+" vs "+ref)
                altSeq=seq[:local]+allele+seq[(local+1):]
                print(cre.ID+"\t"+actualInterval+"\tpos="+str(pos)+"\tref="+\
                      ref+"\t"+allele+"\t"+altSeq,file=INPUTS)
                n+=1
                if(n>=JOBSIZE):
                    INPUTS.close()
                    cmd=CHILD+" "+model+" "+inputFile
                    print(cmd)
                    pipe=Pipe(cmd)
                    while(True):
                        line=pipe.readline()
                        if(line==None): break
                        if(rex.find("ref=",line)): print(line,file=OUT)
                    pipe.close()
                    n=0
                    INPUTS=open(inputFile,"wt")
        creIndex+=1
    if(n>0):
        INPUTS.close()
        cmd=CHILD+" "+model+" "+inputFile
        pipe=Pipe(cmd)
        while(True):
            line=pipe.readline()
            if(line==None): break
            if(rex.find("ref=",line)): print(line,file=OUT)
        pipe.close()
    reader.close()
    OUT.close(); INPUTS.close

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <model> <CREs.txt> <seq-len> <seqs-per-job> <out-file>\n")
(model,creFile,seqLen,JOBSIZE,outputFile)=sys.argv[1:]
seqLen=int(seqLen)
JOBSIZE=int(JOBSIZE)

# Make some temp files
coordFile=TempFilename.generate(".coords")
fastaFile=TempFilename.generate(".fasta")
inputFile=TempFilename.generate(".inputs")

print("Loading CREs",flush=True)
recs=loadCREs(creFile,MAX_N)
print("Writing coords file",flush=True)
localPos=writeCoords(recs,seqLen,coordFile)
print("Extracting genomic sequences",flush=True)
makeFasta(coordFile,fastaFile)
print("Running the model",flush=True)
runModel(recs,seqLen,fastaFile,inputFile,outputFile,JOBSIZE,model)

# Clean up
os.remove(coordFile); os.remove(fastaFile); os.remove(inputFile)



