# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 21:19:29 2018

This script will generate CSV files from a compressed .pkl dataframe for All Mosquitoes,
all hatched mosquitoes (no eggs or gestating agents), or all adult mosquitoes. Line28 can be modified 
to report only adult females, provided the agent.sex information was recorded in the .pkl file

Genotypes reported include total, wild-type, SGD, Resistant, and embryonic lethal

@author: msmanski
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def count_genotypes(genotypeList,StateGenPosData, x, y):
    """Counts how many of each genotype are present at each step"""
    allMos = 0
    nonEggs = 0
    Adults = 0
    for i in range(len(genotypeList)):
        gt = genotypeList[i]
        b = sum(1 for item in StateGenPosData if not 'new' in item[0] and not 'gestating' in item[0] and gt in item[1] and item[2]==[x,y])
        c = sum(1 for item in StateGenPosData if not 'new' in item[0]  and not 'egg' in item[0] and not 'gestating' in item[0] and gt in item[1] and item[2]==[x,y])
        d = sum(1 for item in StateGenPosData if 'A' in item[0] and gt in item[1] and item[2]==[x,y])
        allMos = allMos + b
        nonEggs = nonEggs + c
        Adults = Adults + d
    return allMos, nonEggs, Adults

def processModelData(StateGenPosData, filename, cells, height):
    Xaxis = np.zeros((1,1))
    YaxisTotal_allMos = np.zeros((1,cells))
    YaxisWT_allMos = np.zeros((1,cells))
    YaxisGD_allMos = np.zeros((1,cells))
    YaxisRes_allMos = np.zeros((1,cells))
    YaxisEmbryonicLethal_allMos = np.zeros((1,cells))
    YaxisTotal_nonEggs = np.zeros((1,cells))
    YaxisWT_nonEggs = np.zeros((1,cells))
    YaxisGD_nonEggs = np.zeros((1,cells))
    YaxisRes_nonEggs = np.zeros((1,cells))
    YaxisEmbryonicLethal_nonEggs = np.zeros((1,cells))
    YaxisTotal_Adults = np.zeros((1,cells))
    YaxisWT_Adults = np.zeros((1,cells))
    YaxisGD_Adults = np.zeros((1,cells))
    YaxisRes_Adults = np.zeros((1,cells))
    YaxisEmbryonicLethal_Adults = np.zeros((1,cells))
    for i in range(len(StateGenPosData.axes[0].levels[0])):
        Xaxis = np.vstack((Xaxis,i))
        totalRow_allMos = np.zeros((1,cells))
        WTRow_allMos = np.zeros((1,cells))
        GDRow_allMos = np.zeros((1,cells))
        ResRow_allMos = np.zeros((1,cells))
        EmbryonicLethalRow_allMos = np.zeros((1,cells))
        totalRow_nonEggs = np.zeros((1,cells))
        WTRow_nonEggs = np.zeros((1,cells))
        GDRow_nonEggs = np.zeros((1,cells))
        ResRow_nonEggs = np.zeros((1,cells))
        EmbryonicLethalRow_nonEggs = np.zeros((1,cells))
        totalRow_Adults = np.zeros((1,cells))
        WTRow_Adults = np.zeros((1,cells))
        GDRow_Adults = np.zeros((1,cells))
        ResRow_Adults = np.zeros((1,cells))
        EmbryonicLethalRow_Adults = np.zeros((1,cells))
        for a in range(7):
            for b in range(7):
                col = (height*a)+b
                totalRow_allMos[0,col], totalRow_nonEggs[0,col], totalRow_Adults[0,col] = count_genotypes(['WW','GW','WG','GG','RR','WR','RW','GR','RG','RL','LR','WL','LW','GL','LG','LL'],StateGenPosData.loc[i].values, a, b)
                WTRow_allMos[0,col], WTRow_nonEggs[0,col], WTRow_Adults[0,col] = count_genotypes(['WW'],StateGenPosData.loc[i].values, a, b)
                GDRow_allMos[0,col], GDRow_nonEggs[0,col], GDRow_Adults[0,col] = count_genotypes(['GW','WG'],StateGenPosData.loc[i].values, a, b)
                ResRow_allMos[0,col], ResRow_nonEggs[0,col], ResRow_Adults[0,col] = count_genotypes(['WR','RW','RR','GR','RG'],StateGenPosData.loc[i].values, a, b)
                EmbryonicLethalRow_allMos[0,col], EmbryonicLethalRow_nonEggs[0,col], EmbryonicLethalRow_Adults[0,col] = count_genotypes(['WL','LW','GG','LL','GL','LG','RL','LR'],StateGenPosData.loc[i].values, a, b)

        YaxisTotal_allMos = np.vstack((YaxisTotal_allMos,totalRow_allMos))
        YaxisWT_allMos = np.vstack((YaxisWT_allMos,WTRow_allMos))
        YaxisGD_allMos = np.vstack((YaxisGD_allMos,GDRow_allMos))
        YaxisRes_allMos = np.vstack((YaxisRes_allMos,ResRow_allMos))
        YaxisEmbryonicLethal_allMos = np.vstack((YaxisEmbryonicLethal_allMos,EmbryonicLethalRow_allMos))
        
        YaxisTotal_nonEggs = np.vstack((YaxisTotal_nonEggs,totalRow_nonEggs))
        YaxisWT_nonEggs = np.vstack((YaxisWT_nonEggs,WTRow_nonEggs))
        YaxisGD_nonEggs = np.vstack((YaxisGD_nonEggs,GDRow_nonEggs))
        YaxisRes_nonEggs = np.vstack((YaxisRes_nonEggs,ResRow_nonEggs))
        YaxisEmbryonicLethal_nonEggs = np.vstack((YaxisEmbryonicLethal_nonEggs,EmbryonicLethalRow_nonEggs))

        YaxisTotal_Adults = np.vstack((YaxisTotal_Adults,totalRow_Adults))
        YaxisWT_Adults = np.vstack((YaxisWT_Adults,WTRow_Adults))
        YaxisGD_Adults = np.vstack((YaxisGD_Adults,GDRow_Adults))
        YaxisRes_Adults = np.vstack((YaxisRes_Adults,ResRow_Adults))
        YaxisEmbryonicLethal_Adults = np.vstack((YaxisEmbryonicLethal_Adults,EmbryonicLethalRow_Adults))


    filenameb = filename.split('.')[0] +'_allMos_Total.csv'
    np.savetxt(filenameb,YaxisTotal_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_WT.csv'
    np.savetxt(filenameb,YaxisWT_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_GD.csv'
    np.savetxt(filenameb,YaxisGD_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_Res.csv'
    np.savetxt(filenameb,YaxisRes_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_allMos, delimiter=',')
    
    filenameb = filename.split('.')[0] +'_nonEggs_Total.csv'
    np.savetxt(filenameb,YaxisTotal_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_nonEggs_WT.csv'
    np.savetxt(filenameb,YaxisWT_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_nonEggs_GD.csv'
    np.savetxt(filenameb,YaxisGD_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_nonEggs_Res.csv'
    np.savetxt(filenameb,YaxisRes_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_nonEggs_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_nonEggs, delimiter=',')
    
    filenameb = filename.split('.')[0] +'_Adults_Total.csv'
    np.savetxt(filenameb,YaxisTotal_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_WT.csv'
    np.savetxt(filenameb,YaxisWT_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_GD.csv'
    np.savetxt(filenameb,YaxisGD_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_Res.csv'
    np.savetxt(filenameb,YaxisRes_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_Adults, delimiter=',')
    

filename = "SGDexperiment__.pkl"                     # Replace with .pkl filename for SGD experiment
StateGenPosData = pd.read_pickle(filename)
cells = 3                                                              # Number of populations being simulated
height = 1                                                              # Number of rows in population matrix


processModelData(StateGenPosData, filename, cells, height)

    
