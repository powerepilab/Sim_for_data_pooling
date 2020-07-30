##################################
# ~ Proof of concept example for use of simulation to allow data pooling despite privacy restrictions
# ~ RCTSim2
# Author: Teresa Filshtein
##################################


#### Orientation to the Code####

# Note: Step 1 - Determine the causal structure - occurs prior to any coding

# HRS_sim_data_code and ARIC_sim_data_code: Implements Step 2 - estimates non parametric/parametric models from real data, used for data generating process
# RCTSim1: Collates "rules" generated in Step 2 - run this script to pull in data/model fits/etc from RCTSim0, to be used in RCTSim2
# RCTSim2: Implements Steps 3 through 5 - data generation process, using model fits, summary statistics from RCTSim0, generate data, fit models

###################################

#NOTE:  NEED TO RUN RCTSim1 and have the resulting objects loaded into memory for this to run.

### RCTSim2 consists of 3 functions ###

## Functions

# 1. dataGen(cohort, dfSource) 
# - cohort specifies which cohort to generate data for, 
# - dfSource are the summary stastistics, model fits, etc for cohort. 

# 2. resFunc(modres)
# - simple function to pull results from a model fit

# 3. modFitPool(cohort1 = 'hrs', cohort2 = 'aric')
# - calls dataGen for each cohort listed
# - data is generated for each cohort and then pooled together
# - necessary variables are imputed
# - models are fit for the simulated data and pooled data


#----------------------------------------------------------------------------------------------
### source RCTSim1  - to pull in data for use in modFitPool()

#source('/RCTSim1_publication.R')
# RCTSim1_publication.R is a script that pulls necessary information from summary statistics, models fits, etc.
# It is not necessary, but used to gather information for the data simulation process.
# The output is a list:
#list(
#   n = number of participants, 
#   n_obs = number of observations for each participant,
#   time.int = time between observations,
#   dist_mem = model fit for memory~hba1c,
#   dist_dem = demographics (age category, race, gender proportions),
#   dist_educ = model fit for education,
#   dist_cses = model fit for cses,
#   dist_diet_p = model fit for diet prudent,
#   dist_diet_w = model fit for diet western,
#   dist_hba1c = model fit for hba1c,
#   sdRanInt = estimated sd for intercept random effect,
#   sdRanSlope = estimated sd for slope random effect, 
#   covIntSlope = estimated covariance between int and slope random effects,
#   sdError = estimated sd for error
# ) 

#----------------------------------------------------------------------------------------------

# Function 1: Generates cohort specific data based off of summary statistics/models fits provided. 
dataGen = function(cohort,dfSource){

#######################################################
### Age Category, Gender, Race,  Education and CSES ###
#######################################################

df_sim = get(dfSource)  # pull in data sources (from RCTSim1_publication.R)


#######################################################
### Age Category, Gender, Race,  Education and CSES ###
#######################################################

# df_sim$dist_dem is a data frame where each row is an age category gender and race combindation 
d_dem = sample(x = 1:nrow(df_sim$dist_dem),size = df_sim$n,prob = df_sim$dist_dem$prop,replace = TRUE)

df1 = data.frame(
  ID = 1:df_sim$n, # id for each individual 
  df_sim$dist_dem[d_dem,] # age, gender, race
) %>% select(-c(N, prop)) %>%

### Generate Age, uniform within each category
 mutate(
  Age_Low = 
  case_when(
    agecat == "Agecat0" ~ 50,
    agecat == "Agecat1" ~ 55,
    agecat == "Agecat2" ~ 60,
    agecat == "Agecat3" ~ 65,
    agecat == "Agecat4" ~ 70,
    agecat == "Agecat5" ~ 75,
    agecat == "Agecat6" ~ 80
  ),
  Age_High = 
    case_when(
      agecat == "Agecat0" ~ 55,
      agecat == "Agecat1" ~ 60,
      agecat == "Agecat2" ~ 65,
      agecat == "Agecat3" ~ 70,
      agecat == "Agecat4" ~ 75,
      agecat == "Agecat5" ~ 80,
      agecat == "Agecat6" ~ 100
    )) %>% group_by(ID) %>% mutate(
  base.age = runif(n = 1,Age_Low,Age_High), # generate a random uniform number for each individual within their age category
 base.age60 = (base.age-60)/10 # center at 60, divide by 10 for decades

) %>% add_column(
  ### Generate CSES - function of age category, race
    cses = if(cohort =='hrs'){
      as.numeric(model.matrix(df_sim$dist_cses$form, .) %*% df_sim$dist_cses$fixEF + 
      rnorm(df_sim$n,mean = 0,sd = df_sim$dist_cses$sigma))# noise
    }else{NA}
    )

#################
### Education ###
#################

### compute education probabilities for each individual

educ_prob = data.frame(
  predict(df_sim$dist_educ,newdata = df1,type = 'probs'), 
  U = runif(df_sim$n)
  ) %>% 
  mutate(
    educ = ifelse(U<X1,'College',ifelse(U>=X1 & U<(X1+X2),'HS',"LessThanHS"))
  ) #confirmed X1=College, X2=HS, X3=LtHs

df1$educ = educ_prob$educ



#####################
### Diet         ####
#####################

## Generate diet (Equation 10)


df1 = df1 %>% 
  add_column(
  diet_p = if(cohort == 'aric'){
    as.numeric(model.matrix(df_sim$dist_diet_p$form, .) %*% df_sim$dist_diet_p$fixEF + 
        rnorm(df_sim$n,mean = 0,sd = df_sim$dist_diet_p$sigma))# noise
      
  }else{NA},
  diet_w = if(cohort == 'aric'){
    as.numeric(model.matrix(df_sim$dist_diet_w$form, .) %*% df_sim$dist_diet_w$fixEF + 
      rnorm(df_sim$n,mean = 0,sd = df_sim$dist_diet_w$sigma))# noise
    
  }else{NA}
)

#####################
### HBA1c        ####
#####################

## Generate HbA1C

# Gamma distribution 
# E[hba1c|X] = u =  exp(b0+ b1Age + b2Diet + b3Educ + b5race + b6sex)
# shape = 1/disp.param, scale = u*disp.param
# parents =  "Age"    "Diet"   "Educ"   "race"   "sex" 

disp.param = df_sim$dist_hba1c$disp.param #dispersion parameter
hba1c_pred = data.frame(mu_i = predict(df_sim$dist_hba1c,newdata = df1,type = 'response'),
                       shape_i = 1/disp.param) %>%
  mutate(
   scale_i = mu_i*disp.param
 ) %>% rowwise() %>%
 mutate(
   hba1c = rgamma(n = 1,shape = shape_i, scale = scale_i)
 )

df1$hba1c  = hba1c_pred$hba1c


####################################################################################
## Longitudinal Measurements
####################################################################################

### Expand data - make long 
### Merge with baseline measurements. 
ranef = 
    data.frame(
              rmvnorm(n = df_sim$n, 
                      mean = c(0,0),
                      sigma = matrix(c(df_sim$sdRanInt^2,rep(df_sim$covIntSlope,2),df_sim$sdRanSlope^2),nrow=2)
                      )
              ) %>% 
    rename(int = X1,
           slope  = X2)
    
df1 = df1 %>% add_column(
  int = ranef$int,
  slope = ranef$slope
)


df2 = merge(
  expand.grid(ID = df1$ID,time = seq(from=0,to = (df_sim$n_obs*df_sim$time.int - df_sim$time.int),by = df_sim$time.int)) %>% 
              arrange(ID),
            df1,by = 'ID', all.x = TRUE) %>%
  mutate(
      error = rnorm(n = n(), mean = 0, sd = df_sim$sdError),
      e.age = rnorm(n()),
      age70 = (base.age+time+e.age-70)/10,
      age60 = (base.age+time-60)/10,
      hba1c.c = hba1c-6.5
      ) %>%
  
  ### above and below Hba1c 6.5
    mutate(
      hba1c_O65 = (hba1c.c>=0)*(hba1c.c),
      hba1c_U65 = (hba1c.c<0)*(hba1c.c),
      hba1c_O65only = ifelse(hba1c>=6.5,hba1c,NA),
      hba1c_U65only = ifelse(hba1c<6.5,hba1c,NA)
  ) %>%

  #####################
  #####################
  
  mutate(
      mem = model.matrix(df_sim$dist_mem$final.form, data = .)%*%df_sim$dist_mem$final.fixef + int + (age60)*slope + error  # error changing every time point. 
    )  ###########

df2  = df2 %>% mutate(
    mem.s = as.numeric(mem) 
)

list(df = df2)

}

#-----------------------------------------------------------------------
### Function 2: simple function to pull results from a model fit

## function for pulling results from model fits
resFunc = function(modres){
  data.frame(Est = summary(modres)[["coefficients"]][,'Estimate'],
             SE =summary(modres)[["coefficients"]][,'Std. Error'])
}

#-----------------------------------------------------------------------
### Function 3: implements Steps 4 and 5 described in the mansucript

## Pool data fit (Step 4)
# We recommend setting a seed so that all values can be replicated.
set.seed(19770827)

#cohort1='hrs'
#cohort2='aric'

modFitPool = function(cohort1 = 'hrs', cohort2 = 'aric'){
  
  #specify data source for dataGen()
  dfS1 = paste('RCT2',cohort1, sep = "_")
  dfS2 = paste('RCT2',cohort2, sep = "_")

  #generate data for each cohort
  df1 = dataGen(cohort = cohort1,dfSource = dfS1)
  df2 = dataGen(cohort = cohort2,dfSource = dfS2)

  ### Pool both generated data sets
  
  dfPool = rbind(
    df1$df %>% add_column(
      group = cohort1) %>% 
      mutate(
      groupID = paste(cohort1,ID,sep= "")
      ),
    df2$df %>% add_column(
      group = cohort2) %>% mutate(
      groupID = paste(cohort2,ID,sep= "")
    )
  ) %>% mutate(educ=as.factor(educ))  %>% select(
    group, groupID, ID, 
    sex, race, educ, 
    cses, diet_p, diet_w, 
    hba1c_U65, hba1c_O65,
    age60,time,mem.s
  )
  
  # Impute 'missing' values - impute cses for ARIC and western, prudent for HRS
  
  dfPool_imp = mice(data = dfPool %>% distinct(groupID,.keep_all = TRUE) %>% 
                      select(-c(group, groupID, ID)),
                    m = 1, 
                    method='norm',
                    printFlag=FALSE)
  
  # append groupID to the imputed data
  
  dfPool_complete = cbind(dfPool %>% distinct(groupID,.keep_all = TRUE) %>% 
                            select(groupID), 
                          mice::complete(dfPool_imp)
                          )
  
  # bind imputed values with pooled data set
  
  dfPool_long = left_join(dfPool %>% select(-c(cses, diet_p, diet_w)), 
                          dfPool_complete %>% select(groupID, cses, diet_p, diet_w), by = 'groupID')
 
  #dfPool_long = dfPool_long %>% mutate(educ=as.character(educ))
    
### Analysis (Step 5)
  
## 1. Fit separate models for each cohort with non-imputed data
  
  ## cohort 1  
  mod_HRS = dfPool_long %>% filter(group==cohort1) %>% 
    lmer(mem.s ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex + cses) +
           (1+age60|ID),
         data = ., 
         control = lmerControl(optCtrl = list(maxfun=2e8),
          optimizer="bobyqa"))
  
  ## cohort 2
  mod_ARIC = dfPool_long %>% filter(group==cohort2) %>% 
    lmer(mem.s ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex + diet_p + diet_w) +
           (1+age60|ID),
         data = ., 
         control = lmerControl(optCtrl = list(maxfun=2e8),
                               optimizer="bobyqa"))
  

## 2. Fit separate models for each cohort with imputed data
  
  ## cohort 1
  mod_imp_HRS = dfPool_long %>% filter(group==cohort1) %>% 
    lmer(mem.s ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex + cses + diet_p + diet_w) +
           (1+age60|groupID),
         data = .,
         control = lmerControl(optCtrl = list(maxfun=2e8),
                               optimizer="bobyqa")) 
  
  ## cohort 2
  mod_imp_ARIC = dfPool_long %>% filter(group==cohort2) %>% 
    lmer(mem.s ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex + cses + diet_p + diet_w) +
           (1+age60|groupID),
         data = .,  
         control = lmerControl(optCtrl = list(maxfun=2e8),
                               optimizer="bobyqa"))

## 3. Fit a model for pooled data (with imputed variables)
  
  # 3.a cohort specific intercept and slopes
  mod_imp_POOL = dfPool_long %>%
    lmer(mem.s ~  group*(hba1c_U65 + hba1c_O65) + 
           group*age60*(race + educ + sex + cses + diet_p + diet_w) +
           age60*(hba1c_U65 + hba1c_O65) +
           (1+age60|groupID),
         data = ., 
         control = lmerControl(optCtrl = list(maxfun=2e9),
                               optimizer="bobyqa")) 

  # 3.b no cohort specific effects
  mod_imp_POOL_nogroup1 = dfPool_long %>%
    lmer(mem.s ~  age60*(hba1c_U65 + hba1c_O65 + race + educ + sex + cses + diet_p + diet_w) +
           (1+age60|groupID),
         data = ., 
         control = lmerControl(optCtrl = list(maxfun=2e8),
                               optimizer="bobyqa"))

   
  ## 4. Fit a model for pooled data (without the imputed variables)
  
  mod_imp_POOL_red = dfPool_long %>%
    lmer(mem.s ~  group*(hba1c_U65 + hba1c_O65) + 
           group*age60*(race + educ + sex) + age60*(hba1c_U65 + hba1c_O65) +
           (1+age60|groupID),
         data = ., 
         control = lmerControl(optCtrl = list(maxfun=2e8),
                               optimizer="bobyqa"))

## collect baseline information on simulated data for each cohort
  
  # cohort 1
   tabhrs = df1$df %>% distinct(ID,.keep_all = TRUE) %>% summarise(
    Count = n(),
    Age_mean = mean(base.age),
    Age_sd = sd(base.age),
    Hba1c_mean = mean(hba1c), 
    Hba1c_sd = sd(hba1c),
    hba1c_U65_mean = mean(hba1c_U65only, na.rm=TRUE),
    hba1c_U65_sd = sd(hba1c_U65only, na.rm=TRUE),
    hba1c_O65_mean = mean(hba1c_O65only, na.rm=TRUE),
    hba1c_O65_sd = sd(hba1c_O65only, na.rm=TRUE),
    SexMale = length(ID[sex=='Male']),
    SexFemale = length(ID[sex=='Female']),
    RaceWhite = length(ID[race=='White']),
    RaceBlack = length(ID[race=='Black']),
    EducHS = length(ID[educ=='HS']),
    EducLessHS = length(ID[educ=='LessThanHS']),
    EducCollege = length(ID[educ=='College']),
    cses_mean = mean(cses),
    cses_sd = sd(cses),
    mem_base_mean = mean(mem.s),
    mem_base_sd = sd(mem.s)
  )
  
   # cohort 2
  tabaric = df2$df %>% distinct(ID,.keep_all = TRUE) %>% summarise(
    Count = n(),
    Age_mean = mean(base.age),
    Age_sd = sd(base.age),
    Hba1c_mean = mean(hba1c), 
    Hba1c_sd = sd(hba1c),
    hba1c_U65_mean = mean(hba1c_U65only, na.rm=TRUE),
    hba1c_U65_sd = sd(hba1c_U65only, na.rm=TRUE),
    hba1c_O65_mean = mean(hba1c_O65only, na.rm=TRUE),
    hba1c_O65_sd = sd(hba1c_O65only, na.rm=TRUE),
    SexMale = length(ID[sex=='Male']),
    SexFemale = length(ID[sex=='Female']),
    RaceWhite = length(ID[race=='White']),
    RaceBlack = length(ID[race=='Black']),
    EducHS = length(ID[educ=='HS']),
    EducLessHS = length(ID[educ=='LessThanHS']),
    EducCollege = length(ID[educ=='College']),
    prud_mean = mean(diet_p),
    prud_sd = sd(diet_p),
    west_mean = mean(diet_w),
    west_sd = sd(diet_w),
    mem_base_mean = mean(mem.s),
    mem_base_sd = sd(mem.s)
  )
  
# collected longitudinal data for each simulated data set
tabmem = dfPool_long %>% group_by(group, time) %>% summarise(
  Mean_mem =  mean(mem.s),
  SD_mem = sd(mem.s)
)  

# results
  list(
       HRS_sim = resFunc(mod_HRS),
       ARIC_sim = resFunc(mod_ARIC), 
       HRS_imp = resFunc(mod_imp_HRS),
       ARIC_imp = resFunc(mod_imp_ARIC),
       resPOOL = resFunc(mod_imp_POOL),
       resPOOL_n1 = resFunc(mod_imp_POOL_nogroup1),
       resPOOL_red = resFunc(mod_imp_POOL_red),
       tabhrs = tabhrs,
       tabaric = tabaric,
       tabmem = tabmem

 ) 
}




### Run the simulation B times

#B = 5000 #number of simulations
B = 5000 #number of simulations

full_results = lapply(1:B,function(x){
  print(paste("Simulation ",x))
  modFitPool()
})

save(full_results, file="FILEPAT/full_results_sim5000.rda")

