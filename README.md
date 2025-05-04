# STSRGL

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

To compare with the benchmark algorithms one can edit the 
following files by adding some code line for implentation of 
the benchmark methods

*graph_learning_algorithms.m*:    &emsp;&emsp;&emsp;  for undirected graph learning

*X_recovery_algorithms.m*    :    &emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp; for signal recovery

*AVAR_learning_algorithms.m* :    &emsp;&emsp;&emsp; for directed graph learning or VAR model 
                                estimation


===============================================================

You can run the following demos:
-----------------------------------
* graph_learn_demo/visual_demo:
  
&nbsp;&nbsp;&nbsp; evaluate/visualize the graph learning performance versus different parameters at different scenarios



* signal_recovery_demo:
    
&nbsp;&nbsp;&nbsp; evaluate the signal recovery performance versus different parameters at different scenarios


* AVAR_learn_demo/visual_demo:
  
&nbsp;&nbsp;&nbsp; evaluate/visualize the performance of the VAR model estimation versus different parameters at different scenarios


* converhence_demo:
     	
&nbsp;&nbsp;&nbsp; examine the rate of convergence of the signal recovery or the graph learning methods versus (iteration) time


===============================================================

Please give a star and cite:

A. Javaheri, A. Amini, F. Marvasti and D. P. Palomar, "Learning Spatiotemporal Graphical Models From Incomplete Observations," in IEEE Transactions on Signal Processing, vol. 72, pp. 1361-1374, 2024, doi: 10.1109/TSP.2024.3354572.

A. Javaheri, A. Amini, F. Marvasti and D. P. Palomar, "Joint Signal Recovery and Graph Learning from Incomplete Time-Series," ICASSP 2024 - 2024 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Seoul, Korea, Republic of, 2024, pp. 13511-13515, doi: 10.1109/ICASSP48485.2024.10448021.

