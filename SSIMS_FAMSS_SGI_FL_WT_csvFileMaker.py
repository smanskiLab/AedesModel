# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 08:03:46 2018

@author: msmanski

This script will generate CSV files from a compressed .pkl dataframe for All Mosquitoes,
all hatched mosquitoes (no eggs or gestating agents), or all adult mosquitoes. Line30 can be modified 
to report only adult females, provided the agent.sex information was recorded in the .pkl file

Genotypes reported include total, wild-type, SSIMS/FAMSS(same genotype), SGI (SS), FL, WT-FL hybrid, embyonic lethal,
promoter conversion resistant mutants (SSIMSRes), SGI-FL hybrids, and 'other genotypes'

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



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
    YaxisSSIMS_allMos = np.zeros((1,cells))
    YaxisSS_allMos = np.zeros((1,cells))
    YaxisFL_allMos = np.zeros((1,cells))
    YaxisWT_FLhet_allMos = np.zeros((1,cells))
    YaxisEmbryonicLethal_allMos = np.zeros((1,cells))
    YaxisSSIMSRes_allMos = np.zeros((1,cells))
    YaxisSS_FLhet_allMos = np.zeros((1,cells))
    YaxisOtherGenotype_allMos = np.zeros((1,cells))
    YaxisTotal_nonEggs = np.zeros((1,cells))
    YaxisWT_nonEggs = np.zeros((1,cells))
    YaxisSSIMS_nonEggs = np.zeros((1,cells))
    YaxisSS_nonEggs = np.zeros((1,cells))
    YaxisFL_nonEggs = np.zeros((1,cells))
    YaxisWT_FLhet_nonEggs = np.zeros((1,cells))
    YaxisEmbryonicLethal_nonEggs = np.zeros((1,cells))
    YaxisSSIMSRes_nonEggs = np.zeros((1,cells))
    YaxisSS_FLhet_nonEggs = np.zeros((1,cells))
    YaxisOtherGenotype_nonEggs = np.zeros((1,cells))    
    YaxisTotal_Adults = np.zeros((1,cells))
    YaxisWT_Adults = np.zeros((1,cells))
    YaxisSSIMS_Adults = np.zeros((1,cells))
    YaxisSS_Adults = np.zeros((1,cells))
    YaxisFL_Adults = np.zeros((1,cells))
    YaxisWT_FLhet_Adults = np.zeros((1,cells))
    YaxisEmbryonicLethal_Adults = np.zeros((1,cells))
    YaxisSSIMSRes_Adults = np.zeros((1,cells))
    YaxisSS_FLhet_Adults = np.zeros((1,cells))
    YaxisOtherGenotype_Adults = np.zeros((1,cells))
    for i in range(len(StateGenPosData.axes[0].levels[0])):
        Xaxis = np.vstack((Xaxis,i))
        totalRow_allMos = np.zeros((1,cells))
        WTRow_allMos = np.zeros((1,cells))
        SSIMSRow_allMos = np.zeros((1,cells))
        SSRow_allMos = np.zeros((1,cells))
        FLRow_allMos = np.zeros((1,cells))
        WT_FLhetRow_allMos = np.zeros((1,cells))
        EmbryonicLethalRow_allMos = np.zeros((1,cells))
        SSIMSResRow_allMos = np.zeros((1,cells))
        SS_FLhetRow_allMos = np.zeros((1,cells))
        OtherGenotypeRow_allMos = np.zeros((1,cells))
        totalRow_nonEggs = np.zeros((1,cells))
        WTRow_nonEggs = np.zeros((1,cells))
        SSIMSRow_nonEggs = np.zeros((1,cells))
        SSRow_nonEggs = np.zeros((1,cells))
        FLRow_nonEggs = np.zeros((1,cells))
        WT_FLhetRow_nonEggs = np.zeros((1,cells))
        EmbryonicLethalRow_nonEggs = np.zeros((1,cells))
        SSIMSResRow_nonEggs = np.zeros((1,cells))
        SS_FLhetRow_nonEggs = np.zeros((1,cells))
        OtherGenotypeRow_nonEggs = np.zeros((1,cells))
        totalRow_Adults = np.zeros((1,cells))
        WTRow_Adults = np.zeros((1,cells))
        SSIMSRow_Adults = np.zeros((1,cells))
        SSRow_Adults = np.zeros((1,cells))
        FLRow_Adults = np.zeros((1,cells))
        WT_FLhetRow_Adults = np.zeros((1,cells))
        EmbryonicLethalRow_Adults = np.zeros((1,cells))
        SSIMSResRow_Adults = np.zeros((1,cells))
        SS_FLhetRow_Adults = np.zeros((1,cells))
        OtherGenotypeRow_Adults = np.zeros((1,cells))
        for a in range(height):
            for b in range(np.floor_divide(cells,height)):
                col = (height*a)+b
                totalRow_allMos[0,col], totalRow_nonEggs[0,col], totalRow_Adults[0,col] = count_genotypes(['ppttll','PPTTLL','PPTTll','ppttLL','ppttLl','ppttlL','pPtTll','PptTll','pPTtll','PpTtll','pPtTLl','PptTLl','pPTtLl','PpTtLl','pPtTlL','PptTlL','pPTtlL','PpTtlL','pPtTLL','PptTLL','pPTtLL','PpTtLL','PPTTLl','PPTTlL','PPtTll','pptTll','PPTtll','ppTtll','PPtTLl','pptTLl','PPTtLl','ppTtLl','PPtTlL','pptTlL','PPTtlL','ppTtlL','PPtTLL','pptTLL','PPTtLL','ppTtLL','pPTTll','PpTTll','pPttll','Ppttll','pPTTLl','PpTTLl','pPttLl','PpttLl','pPTTlL','PpTTlL','pPttlL','PpttlL','pPTTLL','PpTTLL','pPttLL','PpttLL','PPttll','ppTTll','PPttLl','ppTTLl','PPttlL','ppTTlL','PPttLL','ppTTLL'],StateGenPosData.loc[i].values, a, b)
                WTRow_allMos[0,col], WTRow_nonEggs[0,col], WTRow_Adults[0,col] = count_genotypes(['ppttll'],StateGenPosData.loc[i].values, a, b)
                SSIMSRow_allMos[0,col], SSIMSRow_nonEggs[0,col], SSIMSRow_Adults[0,col] = count_genotypes(['PPTTLL'],StateGenPosData.loc[i].values, a, b)
                SSRow_allMos[0,col], SSRow_nonEggs[0,col], SSRow_Adults[0,col] = count_genotypes(['PPTTll'],StateGenPosData.loc[i].values, a, b)
                FLRow_allMos[0,col], FLRow_nonEggs[0,col], FLRow_Adults[0,col] = count_genotypes(['ppttLL'],StateGenPosData.loc[i].values, a, b)
                WT_FLhetRow_allMos[0,col], WT_FLhetRow_nonEggs[0,col], WT_FLhetRow_Adults[0,col] = count_genotypes(['ppttLl','ppttlL'],StateGenPosData.loc[i].values, a, b)
                EmbryonicLethalRow_allMos[0,col], EmbryonicLethalRow_nonEggs[0,col], EmbryonicLethalRow_Adults[0,col] = count_genotypes(['pPtTll','PptTll','pPTtll','PpTtll','pPtTLl','PptTLl','pPTtLl','PpTtLl','pPtTlL','PptTlL','pPTtlL','PpTtlL','pPtTLL','PptTLL','pPTtLL','PpTtLL'],StateGenPosData.loc[i].values, a, b)
                SSIMSResRow_allMos[0,col], SSIMSResRow_nonEggs[0,col], SSIMSResRow_Adults[0,col] = count_genotypes(['q'],StateGenPosData.loc[i].values, a, b)
                SS_FLhetRow_allMos[0,col], SS_FLhetRow_nonEggs[0,col], SS_FLhetRow_Adults[0,col] = count_genotypes(['PPTTLl','PPTTlL'],StateGenPosData.loc[i].values, a, b)
                OtherGenotypeRow_allMos[0,col], OtherGenotypeRow_nonEggs[0,col], OtherGenotypeRow_Adults[0,col] = count_genotypes(['PPtTll','pptTll','PPTtll','ppTtll','PPtTLl','pptTLl','PPTtLl','ppTtLl','PPtTlL','pptTlL','PPTtlL','ppTtlL','PPtTLL','pptTLL','PPTtLL','ppTtLL','pPTTll','PpTTll','pPttll','Ppttll','pPTTLl','PpTTLl','pPttLl','PpttLl','pPTTlL','PpTTlL','pPttlL','PpttlL','pPTTLL','PpTTLL','pPttLL','PpttLL','PPttll','ppTTll','PPttLl','ppTTLl','PPttlL','ppTTlL','PPttLL','ppTTLL'],StateGenPosData.loc[i].values, a, b)


        YaxisTotal_allMos = np.vstack((YaxisTotal_allMos,totalRow_allMos))
        YaxisWT_allMos = np.vstack((YaxisWT_allMos,WTRow_allMos))
        YaxisSSIMS_allMos = np.vstack((YaxisSSIMS_allMos,SSIMSRow_allMos))
        YaxisSS_allMos = np.vstack((YaxisSS_allMos,SSRow_allMos))
        YaxisFL_allMos = np.vstack((YaxisFL_allMos,FLRow_allMos))
        YaxisWT_FLhet_allMos = np.vstack((YaxisWT_FLhet_allMos,WT_FLhetRow_allMos))
        YaxisEmbryonicLethal_allMos = np.vstack((YaxisEmbryonicLethal_allMos,EmbryonicLethalRow_allMos))
        YaxisSSIMSRes_allMos = np.vstack((YaxisSSIMSRes_allMos,SSIMSResRow_allMos))
        YaxisSS_FLhet_allMos = np.vstack((YaxisSS_FLhet_allMos,SS_FLhetRow_allMos))
        YaxisOtherGenotype_allMos = np.vstack((YaxisOtherGenotype_allMos,OtherGenotypeRow_allMos))
        
        YaxisTotal_nonEggs = np.vstack((YaxisTotal_nonEggs,totalRow_nonEggs))
        YaxisWT_nonEggs = np.vstack((YaxisWT_nonEggs,WTRow_nonEggs))
        YaxisSSIMS_nonEggs = np.vstack((YaxisSSIMS_nonEggs,SSIMSRow_nonEggs))
        YaxisSS_nonEggs = np.vstack((YaxisSS_nonEggs,SSRow_nonEggs))
        YaxisFL_nonEggs = np.vstack((YaxisFL_nonEggs,FLRow_nonEggs))
        YaxisWT_FLhet_nonEggs = np.vstack((YaxisWT_FLhet_nonEggs,WT_FLhetRow_nonEggs))
        YaxisEmbryonicLethal_nonEggs = np.vstack((YaxisEmbryonicLethal_nonEggs,EmbryonicLethalRow_nonEggs))
        YaxisSSIMSRes_nonEggs = np.vstack((YaxisSSIMSRes_nonEggs,SSIMSResRow_nonEggs))
        YaxisSS_FLhet_nonEggs = np.vstack((YaxisSS_FLhet_nonEggs,SS_FLhetRow_nonEggs))
        YaxisOtherGenotype_nonEggs = np.vstack((YaxisOtherGenotype_nonEggs,OtherGenotypeRow_nonEggs))
        
        YaxisTotal_Adults = np.vstack((YaxisTotal_Adults,totalRow_Adults))
        YaxisWT_Adults = np.vstack((YaxisWT_Adults,WTRow_Adults))
        YaxisSSIMS_Adults = np.vstack((YaxisSSIMS_Adults,SSIMSRow_Adults))
        YaxisSS_Adults = np.vstack((YaxisSS_Adults,SSRow_Adults))
        YaxisFL_Adults = np.vstack((YaxisFL_Adults,FLRow_Adults))
        YaxisWT_FLhet_Adults = np.vstack((YaxisWT_FLhet_Adults,WT_FLhetRow_Adults))
        YaxisEmbryonicLethal_Adults = np.vstack((YaxisEmbryonicLethal_Adults,EmbryonicLethalRow_Adults))
        YaxisSSIMSRes_Adults = np.vstack((YaxisSSIMSRes_Adults,SSIMSResRow_Adults))
        YaxisSS_FLhet_Adults = np.vstack((YaxisSS_FLhet_Adults,SS_FLhetRow_Adults))
        YaxisOtherGenotype_Adults = np.vstack((YaxisOtherGenotype_Adults,OtherGenotypeRow_Adults))

    filenameb = filename.split('.')[0] +'_allMos_Total.csv'
    np.savetxt(filenameb,YaxisTotal_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_WT.csv'
    np.savetxt(filenameb,YaxisWT_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_SSIMS.csv'
    np.savetxt(filenameb,YaxisSSIMS_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_SS.csv'
    np.savetxt(filenameb,YaxisSS_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_FL.csv'
    np.savetxt(filenameb,YaxisFL_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_WT_FLhet.csv'
    np.savetxt(filenameb,YaxisWT_FLhet_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_SSIMSRes.csv'
    np.savetxt(filenameb,YaxisSSIMSRes_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_SS_FLhet.csv'
    np.savetxt(filenameb,YaxisSS_FLhet_allMos, delimiter=',')
    filenameb = filename.split('.')[0] +'_allMos_OtherGenotype.csv'
    np.savetxt(filenameb,YaxisOtherGenotype_allMos, delimiter=',')

    filenameb = filename.split('.')[0] +'_non_eggs_Total.csv'
    np.savetxt(filenameb,YaxisTotal_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_WT.csv'
    np.savetxt(filenameb,YaxisWT_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_SSIMS.csv'
    np.savetxt(filenameb,YaxisSSIMS_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_SS.csv'
    np.savetxt(filenameb,YaxisSS_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_FL.csv'
    np.savetxt(filenameb,YaxisFL_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_WT_FLhet.csv'
    np.savetxt(filenameb,YaxisWT_FLhet_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_SSIMSRes.csv'
    np.savetxt(filenameb,YaxisSSIMSRes_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_SS_FLhet.csv'
    np.savetxt(filenameb,YaxisSS_FLhet_nonEggs, delimiter=',')
    filenameb = filename.split('.')[0] +'_non_eggs_OtherGenotype.csv'
    np.savetxt(filenameb,YaxisOtherGenotype_nonEggs, delimiter=',')

    filenameb = filename.split('.')[0] +'_Adults_Total.csv'
    np.savetxt(filenameb,YaxisTotal_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_WT.csv'
    np.savetxt(filenameb,YaxisWT_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_SSIMS.csv'
    np.savetxt(filenameb,YaxisSSIMS_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_SS.csv'
    np.savetxt(filenameb,YaxisSS_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_FL.csv'
    np.savetxt(filenameb,YaxisFL_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_WT_FLhet.csv'
    np.savetxt(filenameb,YaxisWT_FLhet_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_EmbryonicLethal.csv'
    np.savetxt(filenameb,YaxisEmbryonicLethal_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_SSIMSRes.csv'
    np.savetxt(filenameb,YaxisSSIMSRes_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_SS_FLhet.csv'
    np.savetxt(filenameb,YaxisSS_FLhet_Adults, delimiter=',')
    filenameb = filename.split('.')[0] +'_Adults_OtherGenotype.csv'
    np.savetxt(filenameb,YaxisOtherGenotype_Adults, delimiter=',')


filename = 'FAMSSResistance__.pkl'
StateGenPosData = pd.read_pickle(filename)
cells = 3                                                              # Number of populations being simulated
height = 1        
processModelData(StateGenPosData, filename, cells, height)