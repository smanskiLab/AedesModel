# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:03:27 2018

@author: smanski

This code will simulate a 7x7 grid of neighboring Aedes aegypti populations.

"""
# Main model parameters
experiment = 0                                      # This is used for Filename generation
typeExp = 'Figure_4aSSIMS_repeatReleaseGridREDO'    # This is used for Filename generation 
startingGMOgenotype = 'PPTTLL'                      # Genotype of GE organisms
wildtypeGenotype = 'ppttll'                         # 'ppttll' or 'WW' for SSIMS vs SGD models
HomingFreq = .98
NHEJfreq = .025
promoterConversionFreq = 0.001
startingWT = 500                                    # wild-type agents placed at time step 0
startingGMO = 4500                                  # fGE agents placed at time step 0
addedGMO = 200                                      # intermittently added GE agents
timeStepAdded = 40                                  # number of time steps between intermittent addition
timeSteps = 322                                     # Duration of simulation in steps
migrationRate = 0.01                                # rate of migration to neighboring cell per adult per timestep

#Secondary Model parameters
sexRatio = 0.50

gridWidth = 7                                       #number of columns of population matrix
gridHeight = 7                                      #number of rows of population matrix
borderCells = 2                                     #with of 'Wild-type only' border

import random
import matplotlib.pyplot as plt
import numpy as np
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
import time

gestating=[[],[],[],[],[],[]]                       #Stores gestating eggs until they are ready to be born into model

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[{}]'.format(self.name))
        print('Elapsed: {}'.format(time.time() - self.tstart))

class AedesMosq(Agent):
    """An Aedes aegypti mosquito"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.lifeStage = 1                          
        self.state = 'egg'
        if random.random() < sexRatio:                  
            self.sex = 'male'
        else:
            self.sex = 'female'
        self.pregnant = False
        self.mated = False
        self.genotype = 'ppttll'
        self.gestating = 0
        self.mateGenotype = ''
        self.model.zzz += 1

        
    def mature(self):                                   # Progression from one life=stage to the next at each time step
        if 0 < self.gestating < 8:                      # prevents females from mating to soon after previous mating
            self.gestating = self.gestating + 1
        elif self.gestating == 8:
            self.gestating = 0
            self.pregnant == False
        else:
            self.gestating = 0
            
        if -7 < self.lifeStage < 0:                     # gestating eggs in mother
            self.lifeStage = self.lifeStage + 1
            self.state = 'gestating'
        elif self.lifeStage == 0:
            self.lifeStage = 1
            self.state = 'egg'
            
        elif self.lifeStage == 1:
            if 'p' in self.genotype and 'T' in self.genotype:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
            elif 'L' in self.genotype and self.sex == 'female':
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
            elif 'GG' in self.genotype:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
            elif 'G' in self.genotype and 'I' in self.genotype:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
            elif random.random() < self.model.Sone[self.pos[0],self.pos[1]]: # Time step updated Parameter = ddm k1 [Dye 1984]
                self.lifeStage = 2
                self.state = 'l2'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        
        elif self.lifeStage == 2:
            if random.random() < self.model.Sone[self.pos[0],self.pos[1]]: # Time step updated Parameter = ddm k1 [Dye 1984]
                self.lifeStage = 3
                self.state = 'l3'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 3:
            if random.random() < self.model.Stwo[self.pos[0],self.pos[1]]:  # Time step updated Parameter = ddm k2 [Dye 1984]
                self.lifeStage = 4
                self.state = 'l4'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 4:
            if random.random() < self.model.Stwo[self.pos[0],self.pos[1]]:  # Time step updated Parameter = ddm k2 [Dye 1984]
                self.lifeStage = 5
                self.state = 'l5'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 5:
            if random.random() < 0.765:                                     # Parameter = dis k3 [Dye 1984]
                self.lifeStage = 6
                self.state = 'l6'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 6:
            if random.random() < 0.765:                                     # Parameter = dis k3 [Dye 1984]
                self.lifeStage = 7
                self.state = 'l7'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 7:
            if random.random() < 0.765:                                     # Parameter = dis k3 [Dye 1984]
                self.lifeStage = 8
                self.state = 'l8'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 8:
            if random.random() < 0.754:                                     # Parameter = dis k4 [Dye 1984]
                self.lifeStage = 9
                self.state = 'l9'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 9:
            if random.random() < 0.754:                                     # Parameter = dis k4 [Dye 1984]
                self.lifeStage = 10
                self.state = 'l10'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 10:
            if random.random() < 0.754:                                     # Parameter = dis k4 [Dye 1984]
                self.lifeStage = 11
                self.state = 'pupa'
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif 10 < self.lifeStage < 37:
            if random.random() < 0.95:                                      # Parameter = adult per TS survival [Goindin]
                self.lifeStage = self.lifeStage + 1
                self.state = 'A' + str(self.lifeStage - 11)
            else:
                self.lifeStage = 38
                self.state = 'dead'
                self.model.grid._remove_agent(self.pos, self)
                self.model.schedule.remove(self)
        elif self.lifeStage == 37:
            self.lifeStage = 38
            self.state = 'dead'
            self.model.grid._remove_agent(self.pos, self)
            self.model.schedule.remove(self)
        
        
            
    
    def move(self):                                                     # Moves mosquitoes between neighboring populations
        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            moore=True,                                                 # includes diagonol moves, otherwise von neumann
            include_center=False)
        new_position = random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)
        self.pos = list(new_position)
        
    def reproduce(self, maleGenotype, female):                          # Determine new genotype of each new egg
        for i in range(40):                                             # Parameter Number of eggs per GC
            egg = AedesMosq(self.model.zzz,self.model)
            if len(maleGenotype) == 6:                                  # Rules for SSIMS
                a,b,c,d,e,f = maleGenotype
                if random.random() < 0.5:
                    malePallele = a
                else:
                    malePallele = b
                if random.random() < 0.5:
                    maleTallele = c
                else:
                    maleTallele = d
                if random.random() < 0.5:
                    maleLallele = e
                else: 
                    maleLallele = f
                femaleGenotype = female.genotype
                t,u,v,x,y,z = femaleGenotype
                if random.random() < 0.5:
                    femalePallele = t
                else:
                    femalePallele = u
                if random.random() < 0.5:
                    femaleTallele = v
                else:
                    femaleTallele = x
                if random.random() < 0.5:
                    femaleLallele = y
                else: 
                    femaleLallele = z
                if malePallele == 'P' and femalePallele == 'p':             # Rules for promoter conversion
                    if random.random() < promoterConversionFreq:
                        femalePallele = 'q'
                    else:
                        femalePallele = femalePallele
                elif malePallele == 'p' and femalePallele == 'P':
                    if random.random() < promoterConversionFreq:
                        malePallele = 'q'
                    else:
                        malePallele = malePallele
                else:
                    malePallele = malePallele
                egg.genotype = malePallele + femalePallele + maleTallele + femaleTallele + maleLallele + femaleLallele
                egg.state = 'gestating'
                egg.lifeStage = 0 # from -6
                self.model.grid.place_agent(egg, self.pos)
                #self.model.schedule.add(egg)                               # New eggs will not be added to scheduler until they pop out of 'gestating' list of lists
                gestating[5].append(egg)                                    # Add all gestating eggs to bottom of gestating list of lists

            elif len(maleGenotype) == 2:                                    # Rules for SGD simulations
                a,b = maleGenotype
                if 'G' in maleGenotype:                                     # Decide on inherited genotype from drive allele
                    if 'W' in maleGenotype:
                        randNum = random.random()
                        if randNum < np.add(0.5,np.multiply(0.5,np.multiply(HomingFreq,np.subtract(1,NHEJfreq)))):
                            maleGDallele = 'G'
                        elif np.add(0.5,np.multiply(0.5,np.multiply(HomingFreq,np.subtract(1,NHEJfreq)))) < randNum < np.add(0.5,np.multiply(0.5,HomingFreq)):
                            if random.random() < 0.33:
                                maleGDallele = 'R'
                            else:
                                maleGDallele = 'I'
                        else:
                            maleGDallele = 'W'
                    elif 'R' in maleGenotype:
                        if random.random() < 0.5:
                            maleGDallele = 'R'
                        else:
                            maleGDallele = 'G'
                    else:
                        raise NameError('invalid Male Genotype')
                else:
                    if random.random() < 0.5:
                        maleGDallele = a
                    else:
                        maleGDallele = b
                femaleGenotype = female.genotype
                a,b = femaleGenotype
                if 'G' in femaleGenotype:                                   # Decide on inherited genotype from drive allele
                    if 'W' in femaleGenotype:
                        randNum = random.random()
                        if randNum < np.add(0.5,np.multiply(0.5,np.multiply(HomingFreq,np.subtract(1,NHEJfreq)))):
                            femaleGDallele = 'G'
                        elif np.add(0.5,np.multiply(0.5,np.multiply(HomingFreq,np.subtract(1,NHEJfreq)))) < randNum < np.add(0.5,np.multiply(0.5,HomingFreq)):
                            if random.random() < 0.33:
                                femaleGDallele = 'R'
                            else:
                                femaleGDallele = 'I'
                        else:
                            femaleGDallele = 'W'
                    elif 'R' in femaleGenotype:
                        if random.random() < 0.5:
                            femaleGDallele = 'R'
                        else:
                            femaleGDallele = 'G'
                    else:
                        raise NameError('invalid feMale Genotype')
                else:
                    if random.random() < 0.5:
                        femaleGDallele = a
                    else:
                        femaleGDallele = b
                egg.genotype = maleGDallele + femaleGDallele
                egg.state = 'gestating'
                egg.lifeStage = 0 # from -6
                self.model.grid.place_agent(egg, self.pos)
                #self.model.schedule.add(egg)                               # New eggs will not be added to scheduler until they pop out of 'gestating' list of lists
                gestating[5].append(egg)                                    # Add all gestating eggs to bottom of gestating list of lists
            else:
                raise NameError('incorrect alleleNumber in male.genotype')
        self.model.num_agents = self.model.num_agents + 40       
        
    def female_mate(self):                                                  # Mating-ready females find random mating partner male
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        matingPartners = []
        if len(cellmates) > 1:
            for i in range(len(cellmates)):
                if cellmates[i].sex == 'male' and 'A' in cellmates[i].state: # Will mate with any adult male
                    matingPartners.append(cellmates[i])
                else:
                    i=i+1
    
            if len(matingPartners) > 0:    
                matingPartner = random.choice(matingPartners)
                self.mateGenotype = matingPartner.genotype
                self.reproduce(self.mateGenotype, self)
                matingPartner.mated = False                                 # males can mate multiple times
                self.mated = True
                self.pregnant = True
                self.gestating = 1
                    
            else:
                return
        else:
            return
        
           
    def step(self):
        if self.state == 'egg':                                             # Reduce code runtime by having agents self report to eggmatrix and L1/L2 matrix.
            self.model.tempEggMatrix[self.pos[0]][self.pos[1]] += 1
        elif self.state in ['l2','l3','l4']:
            self.model.tempL1L2Matrix[self.pos[0]][self.pos[1]] += 1
        
        if self.lifeStage == 13 and self.sex == 'female':                   # Females at life stage 13 mate with random male
            self.female_mate()
            self.mature()
        elif self.lifeStage == 21 and self.mated == True and self.sex == 'female':
            self.reproduce(self.mateGenotype,self)
            self.mature()
        elif self.lifeStage == 29 and self.mated ==True and self.sex == 'female':
            self.reproduce(self.mateGenotype,self)
            self.mature()
        elif self.lifeStage > 13:
            if random.random() < self.model.migration:
                self.move()
                self.mature()
            else:
                self.mature()
        else:
            self.mature()
        
                 
class AedesModel(Model):
    """A model with some number of agents."""
    def __init__(self, N, width, height):
        self.num_agents = N
        self.grid = MultiGrid(width, height, True)
        self.schedule = RandomActivation(self)
        self.eggMatrix = np.zeros((width,height))
        self.L1L2Matrix = np.zeros((width,height))
        self.tempEggMatrix = np.zeros((width,height))
        self.tempL1L2Matrix = np.zeros((width,height))
        self.Kone = np.zeros((width,height))
        self.Sone = np.ones((width,height))
        self.Ktwo = np.zeros((width,height))
        self.Stwo = np.ones((width,height))
        self.migration = 0
        self.eggMatrixList = []
        self.L1L2MatrixList=[]
        self.zzz = 0
        for i in range(0,322):                                              # change with timesteps
            self.eggMatrixList.append(np.zeros((width,height)))
            self.L1L2MatrixList.append(np.zeros((width,height)))

        for i in range(self.num_agents):
            a = AedesMosq(i, self)
            self.schedule.add(a)
            
            # Add the agent to a random grid cell, this step is overridden if agent is placed to specific cell
            x = random.randrange(self.grid.width)
            y = random.randrange(self.grid.height)
            self.grid.place_agent(a,(x,y))
            
        self.datacollector = DataCollector(
            agent_reporters={#"Sex": lambda a: a.sex,                       # Commented-out to reduce RAM requirements
                #"Pregnant": lambda a: a.pregnant,
                #"Mated": lambda a: a.mated,
                "State": lambda a: a.state,
                "Genotype": lambda a: a.genotype,
                "Position": lambda a: a.pos})
 
    def step(self):
        '''Advance the model by one stop.'''
        newGeneration = gestating.pop(0)                                    # Pop gestating list of lists, shiftin each list up one position
        gestating.append([])
        for ag in newGeneration:
            self.schedule.add(ag)                                           # add new generation from gestation to schedule
            
        self.datacollector.collect(self)
        self.eggMatrix = self.tempEggMatrix
        self.L1L2Matrix = self.tempL1L2Matrix
        
        self.tempEggMatrix = self.eggMatrixList.pop()                       
        self.tempL1L2Matrix = self.L1L2MatrixList.pop()
        
        self.Kone = np.multiply(0.031,np.power(self.eggMatrix,0.523))       # Equations for ddm and dds from [Dye 1984]
        self.Sone = np.sqrt(np.divide(1,np.power(10,self.Kone)))
        self.Ktwo = np.multiply(0.0003331,np.power(self.L1L2Matrix,1.177))
        self.Stwo = np.sqrt(np.divide(1,np.power(10,self.Ktwo)))
        self.schedule.step()
        print("Model step # ")
        

#RunParameters

filename = typeExp + 'experiment' + str(experiment) +'.pkl'
model = AedesModel(0,gridHeight,gridWidth)
cells = gridWidth * gridHeight
model.migration = migrationRate
zzz=1                                                                       # zzz tracks number of agents, increases by one with each birth

                                                                            # put WT mosquitoes in each cell
for i in range(gridHeight):
    for j in range(gridWidth):
        for h in range(startingWT):
            new = AedesMosq(model.zzz,model)
            new.genotype = wildtypeGenotype
            new.state = 'new'
            if random.random() < 0.95:                                      # 95% of seed organisms are larval or earlier life-stages
                new.lifeStage = np.random.randint(-6,10+1)
            else:
                new.lifeStage = np.random.randint(11,37+1)                  # 5% of seed organisms are adults
            model.grid.place_agent(new,[i,j])
            model.schedule.add(new)
        print("Seeded ",startingWT," wildtype in [", i, ",", j, "]")

    
                                                                            # Put GMO mosquitoes in center cells   
for i in range(gridHeight-(2*borderCells)):
    for j in range(gridWidth-(2*borderCells)):
        for h in range(startingGMO):
            new = AedesMosq(model.zzz,model)
            new.genotype = startingGMOgenotype
            new.state = 'new'
            if random.random() < 0.95:
                new.lifeStage = np.random.randint(-6,10+1)
            else:
                new.lifeStage = np.random.randint(11,37+1)
            model.grid.place_agent(new, [i+borderCells,j+borderCells])
            model.schedule.add(new)
        print("Seeded ",startingGMO," GMO in [", i+borderCells, ",", j+borderCells, "]")

timeList=[]                                                                 # Model counter for debugging purposes
with Timer('FULL EXECUTION'):            
    for i in range(timeSteps):
        
        print("time step ", i, "of ", timeSteps)                            # Print intermediate model statistics
        print("Number of eggs = ", model.eggMatrix)
        print("Kone = ",model.Kone, " Sone = ", model.Sone)
        print("Number if L1L2instars = ", model.L1L2Matrix)
        print("Ktwo = ",model.Ktwo, " Stwo = ", model.Stwo)
        print("NUMAGENTS={}".format(model.num_agents))
                                                                            # Adding new GE organisms and saving compressed datafile at each multiple of timestep added
        if i > 0 and i % timeStepAdded == 0:
            for j in range(gridHeight-(2*borderCells)):
                for k in range(gridWidth-(2*borderCells)):
                    for l in range(addedGMO):
                        new =  AedesMosq(model.zzz,model)
                        new.genotype = startingGMOgenotype
                        new.state = 'new'
                        new.lifeStage = np.random.randint(4,11+1)
                        model.grid.place_agent(new,[j+borderCells,k+borderCells])
                        model.schedule.add(new)
            print("Added ", addedGMO," GMO to each treatment cell")
            model.datacollector.get_agent_vars_dataframe()[['State','Genotype','Position']].copy().to_pickle(filename)
            with Timer('model step {}'.format(i)):
                model.step()
        else:
            with Timer('model step {}'.format(i)):
                model.step()

gestatingCtr=0
eggCtr=0
elseCtr=0
for a in model.schedule.agents:
    if a.state == 'gestating':
        gestatingCtr+=1
    elif a.state == 'egg':
        eggCtr+=1
    else:
        elseCtr+=1
print('gestating: {}'.format(gestatingCtr))
print('egg: {}'.format(eggCtr))
print('elseCtr: {}'.format(elseCtr))
    
mosquito_counts = np.zeros((model.grid.width, model.grid.height))           # This will generate a heatmap with the final number of mosquitoes in each population; not super important but used as a sanity check that the code is working well.
for cell in model.grid.coord_iter():
    cell_content, x, y = cell
    mosquito_count = len(cell_content)
    mosquito_counts[x][y] = mosquito_count
plt.imshow(mosquito_counts, interpolation='nearest')
plt.colorbar()
plt.show()




