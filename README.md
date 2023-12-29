# Channel-Modeling-for-60-GHz-Fixed-mmWave-O2I-and-O2O-Uplink-with-Angular-Misalignment

## Abstract: 
In this paper, we examine the effect of misalignment angle on cluster-based power delay profile modeling for a 60 GHz millimeter-wave (mmWave) uplink. The analysis uses real-world data, where fixed uplink scenarios are realized by placing the transmitter at ground level and the receiver at the building level. Both outdoor-to-indoor (O2I) and outdoor-to-outdoor (O2O) scenarios are studied. Using the misalignment angle and the scenario as inputs, we propose a statistical power delay profile (PDP) simulation algorithm based on the Saleh-Valenzuela (SV) model. Different goodness-of-fit metrics reveal that our proposed algorithm is robust to both O2I and O2O scenarios and can approximate the PDPs fairly well, even in case of misalignment.  

## File descriptions:
1. main_generator: Generate the PDPs for O2I scenario
2. padp_plot: Generates the Angular PDPs for O2I scenario
3. dataset3_analysis_cesa: Generates the PDPs and Angular PDPs for the O2O scenario
   The files above must be run first to obtain the required variables in the workspace
4. svparams: Calculates the SV parameters for all PDPs in O2I and O2O
5. dataset1pdpsim: Calculates the SV parameters for each misalignment range in O2I and simulates a PDP
6. dataset3pdpsim: Calculates the SV parameters for each misalignment range in O2O and simulates a PDP

