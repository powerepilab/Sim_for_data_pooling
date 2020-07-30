##################################
# ~ RCT Why Not, Project 1
# ~ RCTSim0
# Author: Teresa Filshtein
# Date: May 15, 2019
##################################

#################################

#### Orientation to the Code####

# Note: Step 1 - Determine the causal structure - occurs prior to any coding

# HRS_sim_data_code and ARIC_sim_data_code: Implements Step 2 - estimates non parametric/parametric models from real data, used for data generating process
# RCTSim1: Collates "rules" generated in Step 2 - preamble - run this script to pull in data/model fits/etc from RCTSim0, to be used in RCTSim2
# RCTSim2: Implements Steps 3 through 5 - data generation process, using model fits, summary statistics from RCTSim0, generate data, fit models

###################################

### RCTSim0 ###

## Set working directory
setwd("FILEPATH")

# Load packages
packages = c('tidyverse','ggplot2','lme4','purrr','nnet','stargazer','dplyr')
lapply(packages, require,ch = T)

## Load data
load("hrs_cross.rda") #baseline data only, wide dataset
load("hrs_long.rda") #long dataset with structure for linear mixed model of cognitive change


#################
#Table 1
#################

#create variable for HbA1c above and below that is missing if NA
hrs_cross = hrs_cross %>% mutate(
  hba1c_O65only = ifelse(hba1c>=6.5, hba1c, NA),
  hba1c_U65only = ifelse(hba1c<6.5, hba1c, NA)
)

#Table 1 statistics
hrstab1=t(cbind(hrs_cross %>% ungroup() %>% summarise(
  Count=n(),
  Age=paste0(round(mean(base.age),2),"(",round(sd(base.age),2),")"),
  hba1c=paste0(round(mean(hba1c),2),"(",round(sd(hba1c),2),")"),
  HbA1cBelow=paste0(round(mean(hba1c_U65only, na.rm = TRUE),2),"(",round(sd(hba1c_U65only, na.rm = TRUE),2),")"),
  HbA1cAbove=paste0(round(mean(hba1c_O65only, na.rm = TRUE),2),"(",round(sd(hba1c_O65only, na.rm = TRUE),2),")"),
  cses=paste0(round(mean(cses),2),"(",round(sd(cses),2),")"),
  Female=paste0(length(ID[which(sex=="Female")]),"(",round(100*length(ID[which(sex=="Female")])/n(),2),")"),
  White=paste0(length(ID[which(race=="White")]),"(",round(100*length(ID[which(race=="White")])/n(),2),")"),
  Black=paste0(length(ID[which(race=="Black")]),"(",round(100*length(ID[which(race=="Black")])/n(),2),")"),
  LessThanHS=paste0(length(ID[which(educ=="LessThanHS")]),"(",round(100*length(ID[which(educ=="LessThanHS")])/n(),2),")"),
  HS=paste0(length(ID[which(educ=="HS")]),"(",round(100*length(ID[which(educ=="HS")])/n(),2),")"),
  College=paste0(length(ID[which(educ=="College")]),"(",round(100*length(ID[which(educ=="College")])/n(),2),")")
)))

#cognitive assessment stats for 1st assessment
cog1=hrs_long %>% ungroup() %>% filter(assessNo<=1) %>% summarize(average = mean(mem.s), std = sd(mem.s))



write.table(hrstab1, "FILEPATH/HRSTable1.txt", sep="\t")
write.table(cog1, "FILEPATH/HRSTable1_cog.txt", sep="\t")





####################################
# Estimate data Generating rules by generating summary stats/model fits
####################################
n_hrs = nrow(hrs_cross)
n_hrs #11509
##################################
# ~ age, sex, race, educ
##################################

desc.tab0 = hrs_cross %>% group_by(agecat,sex, race) %>% summarise(
  N = n()) %>% ungroup()  %>% 
  mutate(
    prop = prop.table(N)
  )



####################################
# ~ cses (Equation 1)
####################################

jd_cses <- lm(cses~ base.age60 + race, data = hrs_cross)


####################################
# ~ educ (Equation 2 & 3)
####################################

jd_educ <- multinom(educ ~ base.age60 + sex + race + cses, 
                    data = hrs_cross)


####################################
# ~ diet (Equation 4 & 5)
####################################

#Omitted because missing from HRS data

####################################
# ~ hba1c (Equation 6)
####################################
#library(truncnorm)


jd_hba1c <- glm(
  formula = hba1c ~  base.age60 + sex + race + educ + cses,
  family = Gamma(link = "log"), 
data = hrs_cross)


####################################
# ~ memory (Equation 7 and 8)
####################################


jd_mem <-   lmer(I(mem.s) ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex +  cses) + 
         (1 + age60 | ID),
       data = hrs_long)

#save output for tables
sink("HRS_memorymodel.txt")
summary(jd_mem)
sink()

###

#remove input data files
rm(hrs_cross, hrs_long)

####################################################################
## remove individual level data from model fits objects

## NOTE: R OBJECTS INCLUDE EMBEDDED GLOBAL EVIRONMENT AND OTHER INFORMATION THAT MAY HOLD
## RAW DATA.  WE RECOMMEND USING A DIFFERENT APPROACH, MORE EXPLICITLY SAVING OUT WHAT YOU NEED, RATHER THAN REMOVING
## THINGS FROM THE R OBJECT TO ENSURE INDIVIDUAL LEVEL DATA ISN'T INADVERTANTLY SHARED IN THE FUTURE.


jd_cses = list(
  form = formula(jd_cses)[-2],
  fixEF = coefficients(jd_cses),
  sigma = summary(jd_cses)[['sigma']]
)

## create these for aric diet_p, diet_w

# jd_diet_p = list(
#   form = formula(jd_diet_p)[-2],
#   fixEF = coefficients(jd_diet_p),
#   sigma = summary(jd_diet_p)[['sigma']]
# )

# jd_diet_w = list(
#   form = formula(jd_diet_w)[-2],
#   fixEF = coefficients(jd_diet_w),
#   sigma = summary(jd_diet_w)[['sigma']]
# )

rm_educ = c(1:length(jd_educ))[names(jd_educ) %in% c('model', 'residuals','effects','fitted.values', 'weights')]
jd_educ[rm_educ] = NULL

rm_hba1c = c(1:length(jd_hba1c))[names(jd_hba1c) %in% c('model', 'residuals','effects','fitted.values', 'weights', 'linear.predictors','y','data','prior.weights')]
jd_hba1c$disp.param = summary(jd_hba1c)$dispersion
jd_hba1c[rm_hba1c] = NULL
jd_hba1c$qr$qr = c()


jd_mem = list(final.form = formula(jd_mem, fixed.only  = TRUE)[-2],
              final.fixef = fixef(jd_mem),
              reVar = VarCorr(jd_mem)
)


#remove all attributes that contain contents of global environment just in case raw data is accidentally
# still loaded into global environment when you save files
environment(jd_educ[["terms"]]) = NULL
environment(jd_mem[["final.form"]]) = NULL
environment(jd_cses[["form"]])=NULL
environment(jd_hba1c[["terms"]]) = NULL
environment(jd_hba1c[["family"]]) = NULL
environment(jd_hba1c[["formula"]]) = NULL


#remove unnecessary things from the global environment
rm(rm_educ, rm_hba1c)

####################################
# Save Data
####################################

######################################################################################################
save(n_hrs, desc.tab0, jd_educ, jd_cses,jd_hba1c,jd_mem,
     file = 'HRS_sim_data.RDA')
######################################################################################################
