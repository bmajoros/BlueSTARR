#!/usr/bin/env python
#========================================================================
# BlueSTARR Version 0.1
#
# Adapted from DeepSTARR by Bill Majoros (bmajoros@alumni.duke.edu)
#
#========================================================================
import gc
import gzip
import time
import math
import tensorflow as tf
import keras
import keras.layers as kl
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras.layers import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers import BatchNormalization, InputLayer, Input, LSTM, GRU, Bidirectional, Add, Concatenate, LayerNormalization, MultiHeadAttention
import keras_nlp
from keras_nlp.layers import SinePositionEncoding
from keras import models
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint
import keras.backend as backend
from keras.backend import int_shape
import pandas as pd
import numpy as np
import ProgramName
import sys
import IOHelper
import SequenceHelper
import random
from scipy import stats
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
from NeuralConfig import NeuralConfig
from Rex import Rex
rex=Rex()


#========================================================================
#                                GLOBALS
#========================================================================
config=None
NUM_DNA=None # number of DNA replicates
NUM_RNA=None # number of RNA replicates
#RANDOM_SEED=1234
ALPHA={"A":0,"C":1,"G":2,"T":3}
BATCH_SIZE=1

#=========================================================================
#                                main()
#=========================================================================
def main(infile,modelFilestem):
    #startTime=time.time()

    # Load model
    model=None
    with open(modelFilestem+'.json', "r") as json_file:
        model_json=json_file.read()
        model = tf.keras.models.model_from_json(model_json)
        model.load_weights(modelFilestem+'.h5')

    # Load data
    IN=open(infile,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)<6): continue
        ID=fields[0]; ref=fields[1]
        if(not rex.find("ref=(.)",ref)):
            raise Exception("Can't parse ref: "+ref)
        ref=rex[1]
        #alleles=[fields[2],fields[4],fields[6],fields[8]]
        #seqs=[fields[3],fields[5],fields[7],fields[9]]
        alleles=[]; seqs=[]
        i=2
        while(i<len(fields)):
            alleles.append(fields[i])
            seqs.append(fields[i+1])
            i+=2
        Y=[]
        for seq in seqs:
            X=oneHot(seq)
            X=X.reshape((1,X.shape[0],X.shape[1]))
            pred=model.predict(X,batch_size=1,verbose=0)
            Y.append(pred[0][0][0])
            del X
        recs=getScores(ref,alleles,Y)
        line=[ID]
        for rec in recs: line.extend([str(x) for x in rec])
        print("\t".join(line))
        del recs; del fields; del line; del Y; del seqs; del alleles
        del ref; del ID
        gc.collect()

    # Report elapsed time
    #endTime=time.time()
    #seconds=endTime-startTime
    #minutes=seconds/60
    #print("Elapsed time:",round(minutes,2),"minutes")

#========================================================================
#                               FUNCTIONS
#========================================================================
def oneHot(seq):
    L=len(seq)
    X=np.zeros((L,4))
    for i in range(L):
        c=seq[i]
        cat=ALPHA.get(c,-1)
        if(cat>=0): X[i,cat]=1
    return X

def findRef(ref,alleles):
    n=len(alleles)
    for i in range(n):
        if(alleles[i]==ref): return i
    raise Exception("Can't find ref allele")

def getScores(ref,alleles,scores):
    r=findRef(ref,alleles)
    refScore=scores[r]
    n=len(alleles)
    recs=[]
    for i in range(n):
        if(i==r): continue
        log2FC=scores[i]-refScore
        log2FC=round(log2FC,2)
        rec=[alleles[i],log2FC]
        recs.append(rec)
    return recs

#=========================================================================
#                         Command Line Interface
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <model-filestem> <data>\n")
(modelFilestem,infile)=sys.argv[1:]
main(infile,modelFilestem)


