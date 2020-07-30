########################################################################################
## 
########################################################################################
## This file is used to create summary
## data for use by UCSF
########################################################################################

########################################################################################



## Set working directory
setwd("FILEPATH")

## Load data
load("Atable_0.rda")  ##wide dataset - values at baseline - age, sex, race/ethnicity, education, hba1c, diet scores
load("memData2.rda")  ##Long dataset one row per visit with cognitive data, which also includes all the relevant variables from baseline

## Libraries
library(nnet) ## multinom function
library(lme4)
library(tidyverse)

#################
# Table 1 - ARIC
#################

#set HbA1c above and below missing if NA
Atable_0 = Atable_0 %>% mutate(
  hba1c_O65only = ifelse(hba1c>=6.5, hba1c, NA),
  hba1c_U65only = ifelse(hba1c<6.5, hba1c, NA)
)

#Table 1 statistics
tab3=t(cbind(Atable_0 %>% ungroup() %>% summarise(
  Count=n(),
  Age=paste0(round(mean(V2AGE),2),"(",round(sd(V2AGE),2),")"),
  hba1c=paste0(round(mean(hba1c),2),"(",round(sd(hba1c),2),")"),
  HbA1cBelow=paste0(round(mean(hba1c_U65only, na.rm = TRUE),2),"(",round(sd(hba1c_U65only, na.rm = TRUE),2),")"),
  HbA1cAbove=paste0(round(mean(hba1c_O65only, na.rm = TRUE),2),"(",round(sd(hba1c_O65only, na.rm = TRUE),2),")"),
  diet_w=paste0(round(mean(diet_w),2),"(",round(sd(diet_w),2),")"),
  diet_p=paste0(round(mean(diet_p),2),"(",round(sd(diet_p),2),")"),
  Female=paste0(length(ID[which(sex=="Female")]),"(",round(100*length(ID[which(sex=="Female")])/n(),2),")"),
  White=paste0(length(ID[which(race=="White")]),"(",round(100*length(ID[which(race=="White")])/n(),2),")"),
  Black=paste0(length(ID[which(race=="Black")]),"(",round(100*length(ID[which(race=="Black")])/n(),2),")"),
  LessThanHS=paste0(length(ID[which(educ=="LessThanHS")]),"(",round(100*length(ID[which(educ=="LessThanHS")])/n(),2),")"),
  HS=paste0(length(ID[which(educ=="HS")]),"(",round(100*length(ID[which(educ=="HS")])/n(),2),")"),
  College=paste0(length(ID[which(educ=="College")]),"(",round(100*length(ID[which(educ=="College")])/n(),2),")"),
  "Cognitive:Assessment 1"=paste0(round(mean(S2DWRT,na.rm=T),2),"(",round(sd(S2DWRT,na.rm = T),2),")"),
  "Cognitive:Assessment 2"=paste0(round(mean(S4DWRT,na.rm=T),2),"(",round(sd(S4DWRT,na.rm = T),2),")"),
  "Cognitive:Assessment 3"=paste0(round(mean(S5DWRT,na.rm=T),2),"(",round(sd(S5DWRT,na.rm = T),2),")"),
  "Cognitive:Assessment 4"=paste0(round(mean(S6DWRT,na.rm=T),2),"(",round(sd(S6DWRT,na.rm = T),2),")")
  
)))

write.table(tab3, "FILENAME/ARICTable1.txt", sep="\t")


#############
## Non-Parametric
#############

## Check N
n_aric=nrow(Atable_0) #12043

## distribution of agecat, sex, race
desc_tab0 <- Atable_0 %>% group_by(agecat, sex, race) %>% summarise(N=n()) %>% ungroup() %>% 
  mutate(prop=prop.table(N))


#############
## Models
#############


## model for educ
jd_educ <- multinom(educ ~ base.age60 + sex + race, data=Atable_0)


## model for diet_p
jd_diet_p <- lm(diet_p ~ base.age60 + sex + race + educ, data=Atable_0)


## model for diet_w
jd_diet_w <- lm(diet_w ~ base.age60 + sex + race + educ, data=Atable_0)


## model for hba1c
jd_hba1c <- glm(hba1c ~ base.age60 + sex + race + educ + diet_p + diet_w,
                family = Gamma(link="log"),data=Atable_0)



## model for memory
jd_mem <- lmer(mem.s ~ age60 * (hba1c_U65 + hba1c_O65 + race + educ + sex + diet_p + diet_w) + (1 + age60 | ID),
                    data = memData2,
                    na.action = na.omit)


#output memory model in original data
sink("ARIC_memorymodel.txt")
summary(jd_mem)
sink()


####################################################################
## remove individual level data from model fits objects

## NOTE:  WE DO NOT RECOMMENT THIS APPROACH IN THE FUTURE.  WE WOULD RECOMMEND SAVING OUT TEXT-BASED LISTS OF RELEVANT INPUTS,
## AS R OBJECTS OFTEN HOLD MORE THAN IS IMMEDIATELY OBVIOUS TO THE USER.

# remove raw data so it doesn't end up in the R objects
rm(memData2, Atable_0)

#jd_cses = list(
#  form = formula(jd_cses)[-2],
#  fixEF = coefficients(jd_cses),
#  sigma = summary(jd_cses)[['sigma']]
#)

 jd_diet_p = list(
   form = formula(jd_diet_p)[-2],
   fixEF = coefficients(jd_diet_p),
   sigma = summary(jd_diet_p)[['sigma']]
 )

 jd_diet_w = list(
   form = formula(jd_diet_w)[-2],
   fixEF = coefficients(jd_diet_w),
   sigma = summary(jd_diet_w)[['sigma']]
 )

 
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


#remove global environment from attributes just in case
environment(jd_educ[["terms"]]) = NULL
environment(jd_mem[["final.form"]]) = NULL
environment(jd_diet_p[["form"]])=NULL
environment(jd_diet_w[["form"]])=NULL
environment(jd_hba1c[["terms"]]) = NULL
environment(jd_hba1c[["family"]]) = NULL
environment(jd_hba1c[["formula"]]) = NULL

#remove unnecessary things from the global environment
rm(rm_educ, rm_hba1c)

################

#save files for use by UCSF
save(n_aric, desc_tab0, jd_educ, jd_diet_p, jd_diet_w, jd_hba1c, jd_mem,  file = "ARIC_sim_data.RDA")
