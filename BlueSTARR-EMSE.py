#!/usr/bin/env python
#========================================================================
# BlueSTARR-multitask Version 0.1
#
# Adapted from DeepSTARR by Bill Majoros (bmajoros@alumni.duke.edu)
# and Alexander Thomson.
#========================================================================
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
from keras_nlp.layers import SinePositionEncoding, TransformerEncoder, RotaryEmbedding
from keras import models
from keras.models import Sequential, Model
#from keras.layers import RotaryEmbedding
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint
import keras.backend as backend
from keras.backend import int_shape
import pandas as pd
import numpy as np
import ProgramName
import sys
#import IOHelper
import SequenceHelper
import random
from scipy import stats
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
from NeuralConfig import NeuralConfig
from Rex import Rex
rex=Rex()
#from rotary_embedding_tensorflow import apply_rotary_emb, RotaryEmbedding



#========================================================================
#                                GLOBALS
#========================================================================
config=None
NUM_DNA=None # array: numbers of DNA replicates in each cell type
NUM_RNA=None # array: numbers of RNA replicates in each cell type
#RANDOM_SEED=1234
EPSILON=tf.cast(1e-10,tf.float32)


#=========================================================================
#                                main()
#=========================================================================
def main(configFile,subdir,modelFilestem):
    startTime=time.time()
    #random.seed(RANDOM_SEED)
    
    # Load hyperparameters from configuration file
    global config
    config=NeuralConfig(configFile)

    # Load data
    print("loading data",flush=True)
    shouldRevComp=config.RevComp==1
    (X_train_sequence, X_train_seq_matrix, X_train, Y_train, idx_train) = \
        prepare_input("train",subdir,shouldRevComp,config.MaxTrain,config)
    (X_valid_sequence, X_valid_seq_matrix, X_valid, Y_valid, idx_val) = \
        prepare_input("validation",subdir,shouldRevComp,config.MaxTrain,config)
    (X_test_sequence, X_test_seq_matrix, X_test, Y_test, idx_test) = \
        prepare_input("test",subdir,shouldRevComp,config.MaxTest,config) \
        if(config.ShouldTest!=0) else (None, None, None, None)
    seqlen=X_train.shape[1]

    # Build model
    model=BuildModel(seqlen)
    model.summary()

    # Train
    if(config.Epochs>0):
        print("Training...",flush=True)
        print("Training set:",X_train.shape) #,Y_train.shape)
        (model,history)=train(model,X_train,Y_train,X_valid,Y_valid)
        print(history.history)
        print("Done training",flush=True)
        print("loss",history.history['loss'])
        print("val_loss",history.history['val_loss'])
    
    # Save model to a file
    model_json=model.to_json()
    with open(modelFilestem+".json","w") as json_file:
        json_file.write(model_json)
    model.save_weights(modelFilestem+".h5")

    # Test and report accuracy
    if(config.ShouldTest!=0):
        numTasks=len(config.Tasks)
        for i in range(numTasks):
            summary_statistics(X_test,Y_test,"Test",i,numTasks,
                               config.Tasks[i],model,idx_test,modelFilestem)
    print('Min validation loss:', round(min(history.history['val_loss']), 4))

    # Report elapsed time
    endTime=time.time()
    seconds=endTime-startTime
    minutes=seconds/60
    print("Elapsed time:",round(minutes,2),"minutes")


    
def summary_statistics(X, Y, set, taskNum, numTasks, taskName, model, idx, modelFilestem):
    pred = model.predict(X, batch_size=config.BatchSize)
    if (config.useCustomLoss) :
        naiveTheta, cor=naiveCorrelation(Y,pred,taskNum,numTasks) # naiveTheta: normal scale, pred: log scale
        df = pd.DataFrame({'idx':idx, 'true':tf.math.log(naiveTheta),'predicted': pred.squeeze()}) # log scale
        mse = np.mean((df['true'] - df['predicted'])**2)
    else:
        cor=stats.spearmanr(tf.math.exp(pred.squeeze()),tf.math.exp(Y))
        df = pd.DataFrame({'idx':idx, 'true':Y.numpy().ravel(),'predicted': pred.squeeze()}) # log scale
        mse = np.mean((df['true'] - df['predicted'])**2)
    
    print(taskName+" rho=",cor.statistic,"p=",cor.pvalue)
    print(taskName+' mse=', mse)


    
def naiveCorrelation(y_true, y_pred, taskNum, numTasks):
    a=0
    for i in range(taskNum): a+=NUM_DNA[i]+NUM_RNA[i]
    b=a+NUM_DNA[taskNum]
    c=b+NUM_RNA[taskNum]
    DNA=y_true[:,a:b] #+1
    RNA=y_true[:,b:c] #+1
    avgX = tf.reduce_mean(DNA, axis=1)
    avgY = tf.reduce_mean(RNA, axis=1)
    naiveTheta = avgY / avgX
    cor=None
    if(numTasks==1):
        cor=stats.spearmanr(tf.math.exp(y_pred.squeeze()),naiveTheta)
    else:
        cor=stats.spearmanr(tf.math.exp(y_pred[taskNum].squeeze()),naiveTheta)
    return naiveTheta, cor




#========================================================================
#                               FUNCTIONS
#========================================================================
def log(x):
    return tf.math.log(x)

def exp(x):
    return tf.math.exp(x)

def logGam(x):
    return tf.math.lgamma(x)

def logLik(sumX,numX,Yj,logTheta,alpha,beta,numRNA,sumDnaLibs,RnaLibs):
    Y=tf.reduce_sum(Yj,axis=1)
    print("Y=",Y)
    n=numRNA
    print("n=",n)
    #L=tf.reduce_sum(RnaLibs,axis=1)/sumDnaLibs
    #print("L=",L)
    sumRnaLibs=tf.reduce_sum(RnaLibs,axis=1)
    print("sumRnaLibs=",sumRnaLibs)
    X=sumX
    print("X=",X)
    theta=exp(logTheta)
    print("theta=",theta)
    firstRatio=log(sumDnaLibs)+logGam(alpha+Y+1)+logGam(alpha+X-1) -\
        ( log(sumRnaLibs)+logGam(alpha+Y)+logGam(alpha+X) )
    print("firstRatio=",firstRatio)
    secondRatio=2*log(sumDnaLibs)+logGam(alpha+Y+2)+logGam(alpha+X-2) -\
        ( 2*log(sumRnaLibs)+logGam(alpha+Y)+logGam(alpha+X) )
    print("secondRatio=",secondRatio)
    squaredError=theta*theta - 2*theta*exp(firstRatio) + exp(secondRatio)
    print("sqErr=",squaredError)
    return squaredError
#    n=tf.shape(sumX)[0]
#    sumX=tf.tile(tf.reshape(sumX,[n,1]),[1,numRNA])
#    theta=tf.math.exp(logTheta) # assume model is predicting log(theta)
#    n=tf.shape(sumDnaLibs)[0]
#    sumDnaLibs=tf.tile(tf.reshape(sumDnaLibs,[n,1]),[1,numRNA])
#    libRatio=RnaLibs/sumDnaLibs
#    theta=theta*libRatio
#    LL=(sumX+alpha)*log(beta+numX)+logGam(Yj+sumX+alpha)+Yj*log(theta)\
#        -logGam(sumX+alpha)-logGam(Yj+1)-(Yj+sumX+alpha)*log(theta+beta+numX)
#    reduced=tf.reduce_sum(LL,axis=1) # sum logLik across iid replicates
#    return reduced

@tf.autograph.experimental.do_not_convert
def makeClosure(taskNum):
    a=0 # first column of DNA counts
    for i in range(taskNum): a+=NUM_DNA[i]+NUM_RNA[i] # skip previous tasks
    b=a+NUM_DNA[taskNum] # first column of RNA counts
    c=b+NUM_RNA[taskNum] # first column of DNA lib sizes
    d=c+NUM_DNA[taskNum] # first column of RNA lib sizes
    @tf.autograph.experimental.do_not_convert
    def loss(y_true, y_pred):
        global EPSILON
        DNA=y_true[:,a:b]
        RNA=y_true[:,b:c]
        DnaLibs=y_true[:,c:d]
        RnaLibs=y_true[:,d:]
        sumX=tf.reduce_sum(DNA,axis=1)
        sumDnaLibs=tf.reduce_sum(DnaLibs,axis=1)
        LL=-logLik(sumX,b-a,RNA,y_pred,EPSILON,EPSILON,NUM_RNA[taskNum],
                   sumDnaLibs,RnaLibs)
        return LL
    return loss

@tf.autograph.experimental.do_not_convert
def mseClosure(taskNum):
    a=0
    for i in range(taskNum): a+=NUM_DNA[i]+NUM_RNA[i]
    b=a+NUM_DNA[taskNum]
    c=b+NUM_RNA[taskNum]
    @tf.autograph.experimental.do_not_convert
    def loss(y_true, y_pred):
        global EPSILON
        DNA=y_true[:,a:b]+1
        RNA=y_true[:,b:c]+1
        sumX=tf.reduce_sum(DNA,axis=1) # axis=1: sum across replicates
        sumY=tf.reduce_sum(RNA,axis=1) # axis=1: sum across replicates
        naiveTheta=sumY/sumX
        mse=tf.math.reduce_mean(tf.math.square(y_pred-tf.math.log(naiveTheta)),
                                axis=1)
        return mse
    return loss

def generate_complementary_sequence(sequence):
    comp_seq = []
    for b in sequence:
        if b == "A":
            comp_seq.append("T")
        elif b == "T":
            comp_seq.append("A")
        elif b == "C":
            comp_seq.append("G")
        elif b == "G":
            comp_seq.append("C")
        elif b == "N":
            comp_seq.append("N")
        else:
            raise ValueError("Cannot convert base {0} to complement base!".format(b))
    return ''.join(comp_seq)



def loadFasta(fasta_path, as_dict=False,uppercase=False, stop_at=None,
              revcomp=False):
    fastas = []
    seq = None
    header = None
    for r in (gzip.open(fasta_path) if fasta_path.endswith(".gz") else open(fasta_path)):
        if type(r) is bytes: r = r.decode("utf-8")
        r = r.strip()
        if r.startswith(">"):
            if seq != None and header != None:
                fastas.append([header, seq])
                if stop_at != None and len(fastas) >= stop_at:
                    break
            seq = ""
            header = r[1:]
        else:
            if seq != None:
                seq += r.upper() if uppercase else r
            else:
                seq = r.upper() if uppercase else r
    if stop_at != None and len(fastas) < stop_at:
        fastas.append([header, seq])
    elif stop_at == None:
        fastas.append([header, seq])
    if as_dict:
        return {h: s for h, s in fastas}
    if(revcomp):
        for rec in fastas:
            rc=generate_complementary_sequence(rec[1])
            rec[1]=rec[1]+"NNNNNNNNNNNNNNNNNNNN"+rc
    return pd.DataFrame({'location': [e[0] for e in fastas],
                         'idx': [e[0].split(' ')[0] for e in fastas],
                         'sequence': [e[1] for e in fastas]})


def loadCounts(filename,maxCases,config):
    IN=gzip.open(filename) if filename.endswith(".gz") else open(filename)
    header=IN.readline()
    if type(header) is bytes: header = header.decode("utf-8")
    if(not rex.find("DNA=([,\d]+)\s+RNA=([,\d]+)",header)):
        raise Exception("Can't parse counts file header: "+header)
    DNAreps=[int(x) for x in rex[1].split(",")]
    RNAreps=[int(x) for x in rex[2].split(",")]
    numTasks=len(DNAreps)
    linesRead=0
    lines=[]
    for line in IN:
        if type(line) is bytes: line = line.decode("utf-8")
        fields=line.rstrip().split()
        fields=[int(x) for x in fields]
        if(config.useCustomLoss): lines.append(fields)
        else: lines.append(computeNaiveTheta(fields,DNAreps,RNAreps))
        linesRead+=1
        if(linesRead>=maxCases): break
    lines=np.array(lines)
    return (DNAreps,RNAreps,lines)


def computeNaiveTheta(line,DNAreps,RNAreps):
    numTasks=len(DNAreps)
    a=0; rec=[]
    for i in range(numTasks):
        b=a+DNAreps[i]
        c=b+RNAreps[i]
        DNA=line[a:b] 
        RNA=line[b:c]
        avgX=sum(DNA)/DNAreps[i]
        avgY=sum(RNA)/RNAreps[i]  # normalized data
        naiveTheta=float(avgY)/float(avgX)
        rec.append(tf.math.log(naiveTheta)) # log-scale
        a=c
    return rec


def prepare_input(set,subdir,shouldRevComp,maxCases,config):
    # Convert sequences to one-hot encoding matrix
    file_seq = str(subdir+"/" + set + ".fasta.gz")
    input_fasta_data_A=loadFasta(file_seq,uppercase=True,revcomp=shouldRevComp,
                                 stop_at=maxCases)
    sequence_length = len(input_fasta_data_A.sequence.iloc[0])
    seq_matrix_A = SequenceHelper.do_one_hot_encoding(input_fasta_data_A.sequence, sequence_length, SequenceHelper.parse_alpha_to_seq)
    X = np.nan_to_num(seq_matrix_A) # Replace NaN with 0 and inf w/big number
    X_reshaped = X.reshape((X.shape[0], X.shape[1], X.shape[2]))
    (DNAreps,RNAreps,Y)=loadCounts(subdir+"/"+set+"-counts.txt.gz",
                                   maxCases,config)
    global NUM_DNA
    global NUM_RNA
    NUM_DNA=DNAreps
    NUM_RNA=RNAreps
    matrix=pd.DataFrame(Y)
    matrix=tf.cast(matrix,tf.float32)
    return (input_fasta_data_A.sequence, seq_matrix_A, X_reshaped, matrix, input_fasta_data_A.idx)

 
def BuildModel(seqlen):
    # Build model

    # Input layer
    inputLayer=kl.Input(shape=(seqlen,4))
    x=inputLayer

    # Optional convolutional layers
    skip=None
    for i in range(config.NumConv):
        skip=x
        if(config.KernelSizes[i]>=seqlen): continue
        dilation=1 if i==0 else config.DilationFactor
        if(i>0 and config.ConvDropout!=0): x=Dropout(config.DropoutRate)(x)
        x=kl.Conv1D(config.NumKernels[i],
                    kernel_size=config.KernelSizes[i],
                    padding=config.ConvPad,
                    dilation_rate=dilation)(x)
        x=BatchNormalization()(x)
        x=Activation('relu')(x)
        if(config.ConvResidualSkip!=0 and
           i-1>=0 and
           config.NumKernels[i-1]==config.NumKernels[i]):
            #skip=tf.tile(skip,config.NumKernels[i])
            x=Add()([x,skip])
        if(config.ConvPoolSize>1 and seqlen>config.ConvPoolSize):
            x=MaxPooling1D(config.ConvPoolSize)(x)
            seqlen/=config.ConvPoolSize
            
    # Optional Transformer encoder layers
    if(config.NumAttn>0):
        #x=x+keras_nlp.layers.SinePositionEncoding()(x)
        x=x+keras_nlp.layers.RotaryEmbedding()(x)
    for i in range(config.NumAttn):
        skip=x
        x=LayerNormalization()(x)
        #x=MultiHeadAttention(num_heads=config.AttnHeads[i],
        #                     key_dim=config.AttnKeyDim[i])(x,x)
        x = TransformerEncoder(intermediate_dim=config.AttnKeyDim[i],
                               num_heads=config.AttnHeads[i],
                               dropout=config.DropoutRate)(x)
        x=Dropout(config.DropoutRate)(x)
        if(config.AttnResidualSkip!=0):
            x=Add()([x,skip])

    # Global pooling
    if(config.GlobalMaxPool!=0):
        x=MaxPooling1D(int_shape(x)[1])(x)
    if(config.GlobalAvePool!=0):
        x=AveragePooling1D(int_shape(x)[1])(x)
    
    # Flatten
    if(config.Flatten!=0):
        x = Flatten()(x) # Commented out on 3/22/2023

    # dense layers
    if(config.NumDense>0):
        x=Dropout(config.DropoutRate)(x)
    for i in range(config.NumDense):
        x=kl.Dense(config.DenseSizes[i])(x)
        x=BatchNormalization()(x)
        x=Activation('relu')(x)
        x=Dropout(config.DropoutRate)(x)
    
    # Heads per cell type
    tasks=config.Tasks
    outputs=[]; losses=[]
    weights=[float(x) for x in config.TaskWeights]
    numTasks=len(tasks)
    for i in range(numTasks):
        task=tasks[i]
        outputs.append(kl.Dense(1,activation='linear',name=task)(x))
        loss=makeClosure(i) if config.useCustomLoss else "mse" #mseClosure(i)
        losses.append(loss)
    model = keras.models.Model([inputLayer], outputs)
    model.compile(keras.optimizers.Adam(learning_rate=config.LearningRate),
                  run_eagerly=True,
                  #metrics=losses,
                  loss=losses,
                  loss_weights=weights)
    return model



def train(model,X_train,Y_train,X_valid,Y_valid):
    earlyStop=EarlyStopping(patience=config.EarlyStop,monitor="val_loss",
                            restore_best_weights=True)
    history=model.fit(X_train,Y_train,verbose=config.Verbose,
              validation_data=(X_valid,Y_valid),batch_size=config.BatchSize,
              epochs=config.Epochs,callbacks=[earlyStop,History()])
    return (model,history)


#=========================================================================
#                         Command Line Interface
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <parms.config> <data-subdir> <out:model-filestem>\n")
(configFile,subdir,modelFilestem)=sys.argv[1:]
main(configFile,subdir,modelFilestem)

