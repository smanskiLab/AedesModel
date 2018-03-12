# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 08:03:46 2018

@author: msmanski

This code should only be used for SSIMS/SGI/FL experiments (i.e. anything with a 6-letter genotype code).
It will take as input a .pkl file generated from the simulation model and return as outputs (i) graphs of time steps versus
population number for each genotype (a line graph for each subpopulation in the experiment) and (ii) a collection of .csv files 
that can be used to plot individual time steps, etc.

Note: for a 49-grid, 300 time step simulation model, this script takes ~12 hours to complete. It requires an amount of RAM that
is approximately 10x larger than the .pkl file size.

A similar file 'GridProcessorGD.py' should be used for SGD experiments.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = "Figure_4aSSIMS_repeatReleaseGridREDOexperiment0.pkl"                                               # Replace with .pkl filename
StateGenPosData = pd.read_pickle(filename)
experiment = 11                                                         # Used in naming .csv output files
typeExp = '49-grid'                                                     # Used in naming .csv output files
cells = 49                                                              # Number of populations being simulated
height = 7                                                              # Number of rows in population matrix


def count_non_Eggs_genotypes(genotypeList,StateGenPosData, x, y):
    """Counts how many of each genotype are present at each step"""
    a = 0
    for i in range(len(genotypeList)):
        gt = genotypeList[i]
        b = sum(1 for item in StateGenPosData if not 'egg' in item[0] and not 'gestating' in item[0] and gt in item[1] and item[2]==[x,y])
        a = a + b
    return a

def processModelData(StateGenPosData, typeExp, experiment, cells, height):
    Xaxis = np.zeros((1,1))
    YaxisTotal = np.zeros((1,cells))
    YaxisWT = np.zeros((1,cells))
    YaxisSSIMS = np.zeros((1,cells))
    YaxisSS = np.zeros((1,cells))
    YaxisFL = np.zeros((1,cells))
    YaxisWT_FLhet = np.zeros((1,cells))
    YaxisEmbryonicLethal = np.zeros((1,cells))
    YaxisSSIMSRes = np.zeros((1,cells))
    YaxisSS_FLhet = np.zeros((1,cells))
    YaxisOtherGenotype = np.zeros((1,cells))
    for i in range(len(StateGenPosData.axes[0].levels[0])):
        Xaxis = np.vstack((Xaxis,i))
        totalRow = np.zeros((1,cells))
        WTRow = np.zeros((1,cells))
        SSIMSRow = np.zeros((1,cells))
        SSRow = np.zeros((1,cells))
        FLRow = np.zeros((1,cells))
        WT_FLhetRow = np.zeros((1,cells))
        EmbryonicLethalRow = np.zeros((1,cells))
        SSIMSResRow = np.zeros((1,cells))
        SS_FLhetRow = np.zeros((1,cells))
        OtherGenotypeRow = np.zeros((1,cells))
        for a in range(7):
            for b in range(7):
                col = (height*a)+b
                totalRow[0,col] = count_non_Eggs_genotypes(['ppttll','PPTTLL','PPTTll','ppttLL','ppttLl','ppttlL','pPtTll','PptTll','pPTtll','PpTtll','pPtTLl','PptTLl','pPTtLl','PpTtLl','pPtTlL','PptTlL','pPTtlL','PpTtlL','pPtTLL','PptTLL','pPTtLL','PpTtLL','PPTTLl','PPTTlL','PPtTll','pptTll','PPTtll','ppTtll','PPtTLl','pptTLl','PPTtLl','ppTtLl','PPtTlL','pptTlL','PPTtlL','ppTtlL','PPtTLL','pptTLL','PPTtLL','ppTtLL','pPTTll','PpTTll','pPttll','Ppttll','pPTTLl','PpTTLl','pPttLl','PpttLl','pPTTlL','PpTTlL','pPttlL','PpttlL','pPTTLL','PpTTLL','pPttLL','PpttLL','PPttll','ppTTll','PPttLl','ppTTLl','PPttlL','ppTTlL','PPttLL','ppTTLL'],StateGenPosData.loc[i].values, a, b)
                WTRow[0,col] = count_non_Eggs_genotypes(['ppttll'],StateGenPosData.loc[i].values, a, b)
                SSIMSRow[0,col] = count_non_Eggs_genotypes(['PPTTLL'],StateGenPosData.loc[i].values, a, b)
                SSRow[0,col] = count_non_Eggs_genotypes(['PPTTll'],StateGenPosData.loc[i].values, a, b)
                FLRow[0,col] = count_non_Eggs_genotypes(['ppttLL'],StateGenPosData.loc[i].values, a, b)
                WT_FLhetRow[0,col] = count_non_Eggs_genotypes(['ppttLl','ppttlL'],StateGenPosData.loc[i].values, a, b)
                EmbryonicLethalRow[0,col] = count_non_Eggs_genotypes(['pPtTll','PptTll','pPTtll','PpTtll','pPtTLl','PptTLl','pPTtLl','PpTtLl','pPtTlL','PptTlL','pPTtlL','PpTtlL','pPtTLL','PptTLL','pPTtLL','PpTtLL'],StateGenPosData.loc[i].values, a, b)
                SSIMSResRow[0,col] = count_non_Eggs_genotypes(['q'],StateGenPosData.loc[i].values, a, b)
                SS_FLhetRow[0,col] = count_non_Eggs_genotypes(['PPTTLl','PPTTlL'],StateGenPosData.loc[i].values, a, b)
                OtherGenotypeRow[0,col] = count_non_Eggs_genotypes(['PPtTll','pptTll','PPTtll','ppTtll','PPtTLl','pptTLl','PPTtLl','ppTtLl','PPtTlL','pptTlL','PPTtlL','ppTtlL','PPtTLL','pptTLL','PPTtLL','ppTtLL','pPTTll','PpTTll','pPttll','Ppttll','pPTTLl','PpTTLl','pPttLl','PpttLl','pPTTlL','PpTTlL','pPttlL','PpttlL','pPTTLL','PpTTLL','pPttLL','PpttLL','PPttll','ppTTll','PPttLl','ppTTLl','PPttlL','ppTTlL','PPttLL','ppTTLL'],StateGenPosData.loc[i].values, a, b)


        YaxisTotal = np.vstack((YaxisTotal,totalRow))
        YaxisWT = np.vstack((YaxisWT,WTRow))
        YaxisSSIMS = np.vstack((YaxisSSIMS,SSIMSRow))
        YaxisSS = np.vstack((YaxisSS,SSRow))
        YaxisFL = np.vstack((YaxisFL,FLRow))
        YaxisWT_FLhet = np.vstack((YaxisWT_FLhet,WT_FLhetRow))
        YaxisEmbryonicLethal = np.vstack((YaxisEmbryonicLethal,EmbryonicLethalRow))
        YaxisSSIMSRes = np.vstack((YaxisSSIMSRes,SSIMSResRow))
        YaxisSS_FLhet = np.vstack((YaxisSS_FLhet,SS_FLhetRow))
        YaxisOtherGenotype = np.vstack((YaxisOtherGenotype,OtherGenotypeRow))

    i=0
    while i < cells:    
        plt.figure()
        plt.plot(Xaxis,YaxisTotal[:,i])
        plt.plot(Xaxis,YaxisWT[:,i])
        plt.plot(Xaxis,YaxisSSIMS[:,i])
        plt.plot(Xaxis,YaxisSS[:,i])
        plt.plot(Xaxis,YaxisFL[:,i])
        plt.plot(Xaxis,YaxisWT_FLhet[:,i])
        plt.plot(Xaxis,YaxisEmbryonicLethal[:,i])
        plt.plot(Xaxis,YaxisSSIMSRes[:,i])
        plt.plot(Xaxis,YaxisSS_FLhet[:,i])
        plt.plot(Xaxis,YaxisOtherGenotype[:,i])
        plt.legend(['Total', 'WT', 'SSIMS', 'SS','FL','WT_FLhet','Embryonic Lethal','SSIMSResistant','SS_FLhet','Other'], loc='upper left')
        plt.title('Genotype non-eggs at each time step, cell = '+ str(i))
        plt.show()
        i=i+1

    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsTotal.csv'
    np.savetxt(filename,YaxisTotal, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsWT.csv'
    np.savetxt(filename,YaxisWT, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsSSIMS.csv'
    np.savetxt(filename,YaxisSSIMS, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsSS.csv'
    np.savetxt(filename,YaxisSS, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsFL.csv'
    np.savetxt(filename,YaxisFL, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsWT_FLhet.csv'
    np.savetxt(filename,YaxisWT_FLhet, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsEmbryonicLethal.csv'
    np.savetxt(filename,YaxisEmbryonicLethal, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsSSIMSRes.csv'
    np.savetxt(filename,YaxisSSIMSRes, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsSS_FLhet.csv'
    np.savetxt(filename,YaxisSS_FLhet, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsOtherGenotype.csv'
    np.savetxt(filename,YaxisOtherGenotype, delimiter=',')




processModelData(StateGenPosData, typeExp, experiment, cells, height)