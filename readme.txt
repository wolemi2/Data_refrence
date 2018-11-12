Title: Bayesian emulation and calibration of an individual-based model of microbial communities (2018) by O.K. Oyebamiji, D.J. Wilkinson, B. Li, P.G. Jayathilake, P. Zuliani, T.P. Curtis 

This folder contains data and R routines to reproduce the results in the manuscripts.
(a) The codes m1.R, m2.R, m3.R and m4.R - R functions for fitting the four multivariate DLMGP emulators and for producing cross-validation results.
 (i) The codes for plotting all the emulation results in the manuscript and supplementary material are also given in m4.R.

 (b) The folder named "data" - contains all the training (input and output) data for building the emulators.
(i) There are 97 training points (3 non-valid simulation data have been removed). 
(ii) Each training output eg "output_14.csv" is a 96 X 4 matrix and denotes the 14th simulation data. The training input eg "input_14.csv" is a 96 X7 matrix of scaled input parameters.
(iii) There are 11 test data points, each denoted as eg "ytest_1.csv" and "inptest_1.csv" with varying dimension depending on when the simulation terminates.
(iii) The outputs are Biomass concentration, total number of particles, biofilm height and surface roughness.
(iv) There are 7 input parameters namely: carbon substrate, o2, nh4, muhet, muaob, Yhet and Yaob respectively.
(v) There are 11 test points (10 varied simulations plus (1 standard simulation where parameters in Table 1 are kept constant)

(c) The folder named "data2" - contains all the iDynoMiCS (both test and training output) data for model calibration.

(d) The "calibration.R" is an R routine for fitting the calibration model using the iDynoMiCS data.

(e) The "rwish.R" other miscellaneous R functions.
(f) The text file "max2min2.txt" denotes the scales (maximum and minimum) values of input parameters.
NOTE: All output data have been logarithm-transformed and input data scaled between 0 and 1.


##################################################################
Dr. O.K. Oyebamiji
Mathematics and Statistics
Lancaster University
Lancaster
LA1 4YF
UK
E-mails: wolemi2@yahoo.com; o.oyebamiji@lancaster.ac.uk