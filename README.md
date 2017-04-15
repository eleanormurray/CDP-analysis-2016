# CDP-analysis-2016

The repository CDP-analysis-2016 archives the SAS code for adherence-adjusted analyses in the placebo arm of the Coronary Drug Project trial, published in Clinical Trials in 2016. [1] This code is included in the SAS appendix of the Clinical Trials paper. 

This repository contains the SAS programs used for the main analyses to replicate the 1980 CDP findings and update these analyses using inverse probability weighting. SAS 9.4 was used for all analyses. The code appendix contains the following programs:
1.	Program 1: Data management
This program reads in the CDP dataset and defines the variables for the original and updated analyses. It creates two datasets for the placebo arm. Changing the value of ITR to 3 produces the same datasets for the clofibrate arm. 
o	Binary.sas7bdat for crude and baseline adjusted analyses
o	Censwt_ag.sas7bdat for post-randomization adjusted analyses
2.	Program 2: Replication of 1980 analysis
This program estimates the crude and baseline adjustment using the original 1980 definition of adherence. The output gives the first 2 columns of row 1 of Table 2 when using the original adherence coding and the second column of row 2 when using the updated adherence coding. It also outputs the results for Table 1.
3.	Program 3: Updated analyses with baseline variables
This program estimates the crude and baseline adjustment using updated definition of adherence and outputs columns 1 and 3 of row 2 of Table 2. Changing the adherence variable input allows estimation of the third column of row 1 in Table 2. 
4.	Program 4: Updated analysis with post-randomization  covariates
This program requires the restricted cubic spline macro developed by F.E. Harrell[2, 3] in order to run. It estimates the association between adherence and death using adjustment for post-randomization covariates using inverse probability weighting and outputs final column for row 2 of Table 2.

The datasets for the Coronary Drug Project have been submitted to the National Heart, Lung, and Blood Institute (NHLBI), and will be available through application to the NHLBI data center. 


Reference:
1.	Murray EJ, Hernán MA. Adherence adjustment in the Coronary Drug Project: A call for better per-protocol effect estimates in randomized trials. Clinical Trials. 2016; 13(4): 372 – 378.
2.	Harrell FE. %rcspline macro. Clinical Biostatistics, . Duke University Medical Center, 1988.
3.	Devlin TF, Weeks BJ. Spline functions for logistic regression modeling. In:  Proc Eleventh Annual SAS Users Group International,  Atlanta, Georgia,  February 9-12 1986: Cary NC: SAS Institute.
