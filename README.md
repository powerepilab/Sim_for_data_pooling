# Sim_for_data_pooling

README


DIRECTORY OF FILES:

ARIC_sim_data_code_postGitHub_20200728.R
	Reads in the ARIC data and creates/saves the simulation rules for ARIC

HRS_sim_data_code_postGitHub_20200728.R
	Reads in the HRS data and creates/saves the simulation rules for HRS

RCTSim1_postGitHub_20200728.R
	Puts the HRS and ARIC simulation rules into the format for reading into the
	next piece of code

RCTSim2_postGitHub_20200728.R
	Creates the simulated data and runs the desired regression models

RCTSim3_postGitHub_20200728.R
	Summarizes the result over the specified number of simulation runs


R SESSION INFO

This code runs with the following version of R or R packages.  Use of other versions may break or otherwise impact the output of this code.

R version 4.0.1 (2020-06-06)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows Server >= 2012 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] mice_3.9.0    mvtnorm_1.1-1 nnet_7.3-14   purrr_0.3.4   lme4_1.1-23   Matrix_1.2-18
 [7] ggplot2_3.3.2 tibble_3.0.1  tidyr_1.1.0   dplyr_1.0.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.4.6     pillar_1.4.4     compiler_4.0.1   nloptr_1.2.2.1   tools_4.0.1     
 [6] boot_1.3-25      statmod_1.4.34   lifecycle_0.2.0  gtable_0.3.0     nlme_3.1-148    
[11] lattice_0.20-41  pkgconfig_2.0.3  rlang_0.4.6      cli_2.0.2        rstudioapi_0.11 
[16] yaml_2.2.1       withr_2.2.0      generics_0.0.2   vctrs_0.3.1      grid_4.0.1      
[21] tidyselect_1.1.0 glue_1.4.1       R6_2.4.1         fansi_0.4.1      minqa_1.2.4     
[26] magrittr_1.5     backports_1.1.7  scales_1.1.1     ellipsis_0.3.1   splines_4.0.1   
[31] MASS_7.3-51.6    assertthat_0.2.1 colorspace_1.4-1 munsell_0.5.0    broom_0.5.6     
[36] crayon_1.3.4     




