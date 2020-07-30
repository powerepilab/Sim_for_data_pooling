##################################
# ~ Proof of concept example for use of simulation to allow data pooling despite privacy restrictions
# ~ RCTSim1
# Author: Teresa Filshtein
##################################

#### Orientation to the Code####

# Note: Step 1 - Determine the causal structure - occurs prior to any coding
# HRS_sim_data_code and ARIC_sim_data_code: Implements Step 2 - estimates non parametric/parametric models from real data, used for data generating process

# RCTSim1: Collates "rules" generated in Step 2 - run this script to pull in data/model fits/etc from RCTSim0, to be used in RCTSim2
# RCTSim2: Implements Steps 3 through 5 - data generation process, using model fits, summary statistics from RCTSim0, generate data, fit models

###################################

### RCTSim1 ###
# RCTSim1 consists of a single function that puts all the data into the right format to be used in RCTSim2


## Set working directory
setwd("FILEPATH")



# Load packages
packages = c('dplyr','tidyr','tibble','ggplot2','lme4',
             'purrr','nnet','mvtnorm','mice')
lapply(packages, require,ch = T)


# Function preamble_sim
preamble_sim = function(source, # cohort: 'aric', 'hrs'
                        dRDA, # basic information for the simulation
                        memRDA, # var cor matrix, memory model
                        demogRDA, # age, race, gender demographics, proportions
                        educRDA, # model fit for education
                        csesRDA, # model fit for cses
                        hba1cRDA, # model fit for hba1c
                        dietPRDA, # model fit for prudent diet
                        dietWRDA, # model fit for western diet
                        nDF, # # of participants
                        n.obs = 4, # # of observations
                        time.int = 2 # time between observations (yrs)
                        ){

load(dRDA)
  
# Basics
n = get(nDF)
n.obs_sim   = n.obs   # maximum number of longitudinal observations
time.int_sim = time.int    # time interval between observations. 


#############################
##### Common Variables ######
#############################


assign("dist_dem", get(demogRDA))
assign("dist_mem", get(memRDA))
assign("dist_educ", get(educRDA))
assign("dist_hba1c", get(hba1cRDA))



###########################
######## Memory ###########
###########################

varMem = as.data.frame(dist_mem$reVar)

sdRanInt= varMem$sdcor[1]
sdRanSlope = varMem$sdcor[2]
covIntSlope = sdRanInt*sdRanSlope*varMem$sdcor[3]
sdError = varMem$sdcor[which(varMem$grp=='Residual')]

if(source == "hrs"){
  ## unique to HRS
  assign("dist_cses", get(csesRDA))
list(
  n = n, 
  n_obs = n.obs_sim,
  time.int = time.int_sim,
  dist_mem = dist_mem,
  dist_dem = dist_dem,
  dist_educ = dist_educ,
  dist_cses = dist_cses,
  dist_hba1c = dist_hba1c,
  sdRanInt = sdRanInt,
  sdRanSlope = sdRanSlope, 
  covIntSlope = covIntSlope,
  sdError = sdError
)
}else if (source=='aric'){
  ## unique to ARIC
  assign("dist_diet_p", get(dietPRDA))
  assign("dist_diet_w", get(dietWRDA))
list(
    n = n, 
    n_obs = n.obs_sim,
    time.int = time.int_sim,
    dist_mem = dist_mem,
    dist_dem = dist_dem,
    dist_educ = dist_educ,
    dist_diet_p = dist_diet_p,
    dist_diet_w = dist_diet_w,
    dist_hba1c = dist_hba1c,
    sdRanInt = sdRanInt,
    sdRanSlope = sdRanSlope, 
    covIntSlope = covIntSlope,
    sdError = sdError
    
  )  
}

}

## Run function for each data set


RCT2_hrs = preamble_sim(source = "hrs",
                        dRDA = 'HRS_sim_data.RDA',
                        memRDA = "jd_mem", 
                        demogRDA = 'desc.tab0',
                        educRDA = 'jd_educ',
                        csesRDA = 'jd_cses', 
                        hba1cRDA = 'jd_hba1c',
                        dietPRDA = 'jd_diet_p',
                        dietWRDA = 'jd_diet_w',
                        nDF = 'n_hrs',
                        n.obs = 4, time.int = 2)

RCT2_aric = preamble_sim(source = 'aric',
                         dRDA = 'ARIC_sim_data.RDA',
                         memRDA = 'jd_mem',
                        demogRDA = 'desc_tab0',
                        educRDA = 'jd_educ',
                        csesRDA = 'jd_cses', 
                        hba1cRDA = 'jd_hba1c',
                        dietPRDA = 'jd_diet_p',
                        dietWRDA = 'jd_diet_w',
                        nDF = 'n_aric',
                        n.obs = 4, time.int = 2)


