# AedesModel

This code will simulate population control of Aedes aegypti using a variety of genetic population control
methods. Python files provided here are sufficient to run the simulations and to process data stored in the
compressed dataframe (.pkl file). Final results output as .csv files show population numbers for each genotype
in each sub-population.

### Prerequisites
Python 3

Mesa: Agent-based modeling in Python 3+(https://github.com/projectmesa/mesa)

matplotlib

pandas

numpy

## Files
1. GridSimulation.py - Use to run the model simulation.
2. GridProcessorSSIMS.py - Performs data processing, and figure generation for self-sorting incompatible male system.
3. GridProcessorSGD.py - Performs data processing, and figure generation for suppression gene drive.
4. SinglePopulationTriplicate.py - Runs the simulation using a 3 x 1 grid with no migration. May be treated as three independent experiments.

## Running an experiment
Execution of the model may be carried out within an IDE such as Spyder.
Before running the model, modify the below variables within the GridSimulation.py or SinglePopulationTriplicate.py files.
Input an exeriment # and typeExp for unique filename generation.
Input an appropriate starting GMO Genotype and wildtype genotype: For SSIMS/SGI/FL, the wild-type genotype
must be 'ppttll', and the starting GMO genotype will be: 'PPTTLL'/'PPTTll'/'ppttLL' (respectively). For SGD
simulations, the GMO genotype is 'GW' or 'WG' (both are equivalent), and the wildtype genotype is 'WW'.
Input values for Homing frequency and NHEJ frequency (for SGD simulations). For no resistance, input 1 and 0, respectively.
Input values for promoter conversion frequency (SSIMS/SGI simulations only). For no resistance, input 0.
Input the number of wild-type and GMO mosquitos to be seeded at time step 0 ("starting...")
For re-release control strategies, input the number of agents to be released at each intermittent timestep.
For re-release control strategies, input the time-step for spacing of intermittent releases.
NOTE: intermediate and final .pkl files will be generated at each re-release. If no intermittent released
is desired, input the desired number of timesteps here, and input the desired number of timesteps +1 for the
'timeSteps' variable.
Input the total number of timesteps to be tracked
Input a migration rate. This migration rate corresponds to the likelihood that a given adult will migrate
to a neighboring cell at each timestep.

Once simulation variables have been input, run script in a Python3 IDE. Updates are printed at each timestep
to track the total number of agents as will as population-dependent variables such as the density-dependent mortality/
survival rates for eggs and early larva. When the model has finished, a new .pkl file should be generated.
This file can be processed to produce .csv files showing population numbers for each genotype at each
time-step.


## Authors

* **Michael J. Smanski**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Mosquito simulation model was based largely on a similar model described by Christopher 
Dye [Journal of Animal Ecology 1984], with additional run parameters taken from a paper by Goindin D, 
et al. [PLoS ONE 2015]
* Many thanks to the developers of the Mesa platform for simulation models.


