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
from ConfigFile import ConfigFile

class NeuralConfig:
    def __init__(self,filename):
        config=ConfigFile(filename)

        # Slurm parameters
        #self.NumRestarts=int(config.lookupOrDie("NumRestarts"))
        
        # Parameters pertaining to the input
        self.RevComp=int(config.lookupOrDie("RevComp"))

        # Data
        self.MaxTrain=config.lookup("MaxTrain")
        if(self.MaxTrain is None or int(self.MaxTrain)<=0):
            self.MaxTrain=sys.maxsize
        else: self.MaxTrain=int(self.MaxTrain)
        self.MaxTest=config.lookup("MaxTest")
        if(self.MaxTest is None or int(self.MaxTest)<=0):
            self.MaxTtest=sys.maxsize
        else: self.MaxTest=int(self.MaxTest)
        self.ShouldTest=config.lookup("ShouldTest")
        if(self.ShouldTest is None): self.ShouldTest=0
        else: self.ShouldTest=int(self.ShouldTest)
        
        # Convolutional parameters
        self.NumConv=int(config.lookupOrDie("NumConvLayers"))
        self.NumKernels=\
            [int(x) for x in config.lookupListOrDie("NumKernels")]
        self.checkSize("NumKernels",self.NumKernels,self.NumConv)
        self.KernelSizes=\
            [int(x) for x in config.lookupListOrDie("KernelSizes")]
        self.checkSize("KernelSizes",self.KernelSizes,self.NumConv)
        self.ConvResidualSkip=int(config.lookupOrDie("ConvResidualSkip"))
        self.ConvPoolSize=int(config.lookupOrDie("ConvPoolSize"))
        self.DilationFactor=int(config.lookupOrDie("DilationFactor"))
        self.ConvPad=config.lookupOrDie("ConvPad")
        self.ConvDropout=int(config.lookupOrDie("ConvDropout"))

        # Attention parameters
        self.NumAttn=int(config.lookupOrDie("NumAttentionLayers"))
        self.AttnHeads=\
            [int(x) for x in config.lookupListOrDie("AttentionHeads")]
        self.checkSize("AttentionHeads",self.AttnHeads,self.NumAttn)
        self.AttnKeyDim=\
            [int(x) for x in config.lookupListOrDie("AttentionKeyDim")]
        self.checkSize("AttentionKeyDim",self.AttnKeyDim,self.NumAttn)
        self.AttnResidualSkip=int(config.lookupOrDie("AttentionResidualSkip"))

        # Global pool & flatten
        self.GlobalMaxPool=int(config.lookupOrDie("GlobalMaxPool"))
        self.GlobalAvePool=int(config.lookupOrDie("GlobalAvePool"))
        self.Flatten=config.lookup("Flatten")
        if(self.Flatten is None or self.Flatten==""): self.Flatten=0
        else: self.Flatten=int(self.Flatten)
        
        # Dense parameters
        self.NumDense=int(config.lookupOrDie("NumDense"))
        self.DenseSizes=[int(x) for x in config.lookupListOrDie("DenseSizes")]
        self.checkSize("DenseSizes",self.DenseSizes,self.NumDense)

        # Parameters pertaining to output and loss
        self.Tasks=config.lookupList("Tasks")
        self.TaskWeights=config.lookupList("TaskWeights")
        self.useCustomLoss=int(config.lookupOrDie("UseCustomLoss"))
        self.useCustomLoss=self.useCustomLoss!=0
        if(self.useCustomLoss):
            print("Using negative log likelihood loss function")
        else: print("Using MSE loss function with naive estimator")
        
        # Training hyperparameters
        self.Verbose=int(config.lookupOrDie("Verbose"))
        self.BatchSize=int(config.lookupOrDie("BatchSize"))
        self.Epochs=int(config.lookupOrDie("Epochs"))
        self.EarlyStop=int(config.lookupOrDie("EarlyStop"))
        self.DropoutRate=float(config.lookupOrDie("DropoutRate"))
        self.LearningRate=float(config.lookupOrDie("LearningRate"))

    def checkSize(self,name,array,size):
        if(len(array)<size):
            raise Exception(name+" wrong size: expecting "+str(size)+
                            ", found "+str(len(array)))
        
    def dump(self):
        print("GlobalMaxPool=",self.GlobalMaxPool,sep="")
        print("GlobalAvePool=",self.GlobalAvePool,sep="")
        #print("NumRestarts=",self.NumRestarts,sep="")
        print("ConvResidualSkip=",self.ConvResidualSkip,sep="")
        print("ConvPad=",self.ConvPad,sep="")
        print("DilationFactor=",self.DilationFactor,sep="")
        print("ConvPoolSize=",self.ConvPoolSize,sep="")
        print("NumDense=",self.NumDense,sep="")
        print("DenseSizes=",self.DenseSizes,sep="")
        print("NumConv=",self.NumConv,sep="")
        print("KernelSizes=",self.KernelSizes,sep="")
        print("NumKernels=",self.NumKernels,sep="")
        print("NumAttn=",self.NumAttn,sep="")
        print("AttnHeads=",self.AttnHeads,sep="")
        print("AttnKeyDim=",self.AttnKeyDim,sep="")
        print("AttnResidualSkip=",self.AttnResidualSkip,sep="")
        print("Verbose=",self.Verbose,sep="")
        print("RevComp=",self.RevComp,sep="")
        print("BatchSize=",self.BatchSize,sep="")
        print("Epochs=",self.Epochs,sep="")
        print("EarlyStop=",self.EarlyStop,sep="")
        print("DropoutRate=",self.DropoutRate,sep="")
        print("LearningRate=",self.LearningRate,sep="")
        print("ConvDropout=",self.ConvDropout,sep="")
        
#=========================================================================
# main()
#=========================================================================
#if(len(sys.argv)!=2):
#    exit(ProgramName.get()+" <in.config>\n")
#(infile,)=sys.argv[1:]
#config=NeuralConfig(infile)
#config.dump()
