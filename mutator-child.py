#!/usr/bin/env python
#========================================================================
# BlueSTARR Version 0.1
# Adapted from DeepSTARR by Bill Majoros (bmajoros@alumni.duke.edu)
#========================================================================
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import gc
import gzip
import time
import math
import tensorflow as tf
tf.autograph.set_verbosity(1)
import logging
logging.getLogger("tensorflow").setLevel(logging.ERROR)
logging.getLogger('tensorflow').disabled = True
from silence_tensorflow import silence_tensorflow
silence_tensorflow()
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
import sys
import random
from scipy import stats
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
import ProgramName
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
tf.compat.v1.disable_eager_execution()

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
    recs=[]
    for line in IN:
        rec=line.rstrip().split()
        if(len(rec)<6): continue
        recs.append(rec)
    X=oneHot(recs)
    batchSize=len(recs)
    pred=model.predict(X,batch_size=batchSize,verbose=0)
    numRecs=len(recs)
    for i in range(numRecs):
        y=pred[i][0][0] #.numpy()
        rec=recs[i]
        (ID,actualInterval,pos,ref,allele,seq)=rec
        print(ID,actualInterval,pos,ref,allele,y,sep="\t")
        
#========================================================================
#                               FUNCTIONS
#========================================================================
def oneHot(recs):
    firstRec=recs[0]
    (ID,actualInterval,pos,ref,allele,seq)=firstRec
    L=len(seq)
    numRecs=len(recs)
    X=np.zeros((numRecs,L,4))
    for j in range(numRecs):
        rec=recs[j]
        (ID,actualInterval,pos,ref,allele,seq)=rec
        for i in range(L):
            c=seq[i]
            cat=ALPHA.get(c,-1)
            if(cat>=0): X[j,i,cat]=1
    return X

#=========================================================================
#                         Command Line Interface
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <model-filestem> <data>\n")
(modelFilestem,infile)=sys.argv[1:]
main(infile,modelFilestem)


