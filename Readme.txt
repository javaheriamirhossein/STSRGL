===============================================================
This repository contains a stand-alone MATLAB code for 
(undirected/directed) graph learning and signal recovery 
from corrupted spatio-temporal data. The folder 'Utils' 
contains some utility functions (metrics and auxiliary 
functions) and the 'Main' folder contains the principal functions 
including the implementation of the 'STSRGL' method proposed in 
the following paper:


Amirhossein Javaheri, Arash Amini, Farokh Marvasti, Daniel P. Palomar, 
"Learning Spatio-Temporal Graphical Models From Incomplete Observations", 
IEEE Transactions on Signal Processing, 2024.

The 'graph_types.mat' in the 'Data' folder, contains several 
types of graphs to be used for genration of random synthetic data.

--------------------------------------------------------------
To compare with the benchmark algorithms one can edit the 
following files by adding some code line for implentation of 
the benchmark methods

graph_learning_algorithms.m:   for undirected graph learning
X_recovery_algorithms.m    :   for signal recovery
AVAR_learning_algorithms.m :   for directed graph learning or VAR model 
                               estimation


===============================================================

You can run the following demos:
-----------------------------------
graph_learn_demo/visual_demo:   evaluate/visualize the graph learning performance versus 
                        	different parameters at different scenarios


signal_recovery_demo:   	evaluate the signal recovery performance versus 
                        	different parameters at different scenarios


AVAR_learn_demo/visual_demo: 	evaluate/visualize the performance of the VAR model estimation  
                       		versus different parameters at different scenarios


converhence_demo:       	examine the rate of convergence of the signal recovery or  
                        	or the graph learning method versus (iteration) time