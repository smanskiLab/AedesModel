# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 08:03:46 2018

@author: msmanski

This code should only be used for SGD experiments (i.e. anything with a 2-letter genotype code).
It will take as input a .pkl file generated from the simulation model and return as outputs (i) graphs of time steps versus
population number for each genotype (a line graph for each subpopulation in the experiment) and (ii) a collection of .csv files 
that can be used to plot individual time steps, etc.

Note: for a 49-grid, 300 time step simulation model, this script takes ~12 hours to complete. It requires an amount of RAM that
is approximately 10x larger than the .pkl file size.

A similar file 'GridProcessorSSIMS.py' should be used for SSIMS/SGI/FL experiments.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = "Filename.pkl"                                               # Replace with .pkl filename
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
    YaxisGD = np.zeros((1,cells))
    YaxisRes = np.zeros((1,cells))
    YaxisEmbryonicLethal = np.zeros((1,cells))
    for i in range(len(StateGenPosData.axes[0].levels[0])):
        Xaxis = np.vstack((Xaxis,i))
        totalRow = np.zeros((1,cells))
        WTRow = np.zeros((1,cells))
        GDRow = np.zeros((1,cells))
        ResRow = np.zeros((1,cells))
        EmbryonicLethalRow = np.zeros((1,cells))
        for a in range(7):
            for b in range(7):
                col = (height*a)+b
                totalRow[0,col] = count_non_Eggs_genotypes(['WW','GW','WG','GG','RR','WR','RW','GR','RG','RL','LR','WL','LW','GL','LG','LL'],StateGenPosData.loc[i].values, a, b)
                WTRow[0,col] = count_non_Eggs_genotypes(['WW'],StateGenPosData.loc[i].values, a, b)
                GDRow[0,col] = count_non_Eggs_genotypes(['GW','WG'],StateGenPosData.loc[i].values, a, b)
                ResRow[0,col] = count_non_Eggs_genotypes(['WR','RW','RR','GR','RG'],StateGenPosData.loc[i].values, a, b)
                EmbryonicLethalRow[0,col] = count_non_Eggs_genotypes(['WL','LW','GG','LL','GL','LG','RL','LR'],StateGenPosData.loc[i].values, a, b)

        YaxisTotal = np.vstack((YaxisTotal,totalRow))
        YaxisWT = np.vstack((YaxisWT,WTRow))
        YaxisGD = np.vstack((YaxisGD,GDRow))
        YaxisRes = np.vstack((YaxisRes,ResRow))
        YaxisEmbryonicLethal = np.vstack((YaxisEmbryonicLethal,EmbryonicLethalRow))

    i=0
    while i < cells:    
        plt.figure()
        plt.plot(Xaxis,YaxisTotal[:,i])
        plt.plot(Xaxis,YaxisWT[:,i])
        plt.plot(Xaxis,YaxisGD[:,i])
        plt.plot(Xaxis,YaxisRes[:,i])
        plt.plot(Xaxis,YaxisEmbryonicLethal[:,i])
        plt.legend(['Total', 'WT', 'GD', 'Res','Embryonic Lethal'], loc='upper left')
        plt.title('Genotype of hatched Aedes at each time step, cell = '+ str(i))
        plt.show()
        i=i+1

    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsTotal.csv'
    np.savetxt(filename,YaxisTotal, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsWT.csv'
    np.savetxt(filename,YaxisWT, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsGD.csv'
    np.savetxt(filename,YaxisGD, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsRes.csv'
    np.savetxt(filename,YaxisRes, delimiter=',')
    filename = typeExp + 'experiment' + str(experiment) +'_non_eggsEmbryonicLethal.csv'
    np.savetxt(filename,YaxisEmbryonicLethal, delimiter=',')
    

processModelData(StateGenPosData, typeExp, experiment, cells, height)