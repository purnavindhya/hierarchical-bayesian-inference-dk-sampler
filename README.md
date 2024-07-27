# Hierarchical Bayesian Inference with Durbin-Koopman Sampler
This repository contains the code for the paper "Matsuoka, Kodai, Kiyoyuki Kaito, and Masamichi Sogabe. "Bayesian time–frequency analysis of the vehicle–bridge dynamic interaction effect on simple-supported resonant railway bridges." Mechanical Systems and Signal Processing 135 (2020): 106373."
The code for data generation in given in the 'SSI_Data_Generation.m' file and the Kalman filter with Durbin-Koopman smoother code is in 'Filter_Smoother_DK.m' file. 
## Instructions:
1. Run 'SSI_Data_Generation.m' to generate the data. This generates a data file 'approx_z.xls' in the same folder.
2. Run 'Filter_Smoother_DK.m' that reads 'approx_z.xls' file and generates the response of the structure.
