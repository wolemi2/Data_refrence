Title: Bayesian emulation and calibration of an individual-based model of microbial communities
O.K. Oyebamiji, D.J. Wilkinson, B. Li, P.G. Jayathilake, P. Zuliani, T.P. Curtis 

This folder contains data and R routines to reproduce the results in the manuscripts.
(a) codes m1.R m2.R, m3.R and m4.R - R functions for fitting the four multivariate DLMGP emulators
 and for producing crossvalidation results.
 (i) The codes for plotting all the Figures in the manuscript and supplementary material are in m4.R.

 (b) Folder "data" - contain all the training (input and output) data for building the emulators.
(i) There are 97 training points (3 non valid simulation data have been removed). 
(ii) Each output eg "output_14.csv" denotes the 14th simulation data and is a 96 X 4 matrix.
(iii) The outputs are Biomass concentration, total number of particles, biofilm height and surface roughness.
(iv) The inputs are 7 parameters namely: carbon substrate, o2, nh4, muhet, muaob, Yhet and Yaob
(v) There are 11 test points (10 varied simulations plus 1 standard simulation (parameters in Table 1 are kept constant)

(c) Folder "data2" - contain all the Idynomics (input and output) data for model calibration.



###################################################################
Dr. O.K. Oyebamiji
Mathematics and Statistics
Lancaster University
Lancaster
LA1 4YF
UK
E-mails: wolemi2@yahoo.com; o.oyebamiji@lancaster.ac.uk