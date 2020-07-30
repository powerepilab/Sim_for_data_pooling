### Summarise Results ###
library(dplyr)
library(tidyr)
library(purrr)

## load in your full_results file

load("FILEPATH/full_results_sim5000.rda")


#pull HRS model results
resHRS = bind_cols(map(map(full_results,1),1)) %>%
  mutate(Var  = row.names(map(full_results,1)[[1]])) %>%
  select(Var, everything())

#pull ARIC model results
resARIC = bind_cols(map(map(full_results,2),1)) %>%
  mutate(Var  = row.names(map(full_results,2)[[1]])) %>%
  select(Var, everything())

#pull HRS imputed model results
resHRSimp = bind_cols(map(map(full_results,3),1)) %>%
  mutate(Var  = row.names(map(full_results,3)[[1]])) %>%
  select(Var, everything())

#pull ARIC imputed model results
resARICimp = bind_cols(map(map(full_results,4),1)) %>%
  mutate(Var  = row.names(map(full_results,4)[[1]])) %>%
  select(Var, everything())

#pull full pooled model results
resPOOL = bind_cols(map(map(full_results,5),1)) %>%
  mutate(Var  = row.names(map(full_results,5)[[1]])) %>%
  select(Var, everything())

#pull pooled without group interaction results
resPOOLn1 =  bind_cols(map(map(full_results,6),1)) %>%
  mutate(Var  = row.names(map(full_results,6)[[1]])) %>%
  select(Var, everything())

#pool reduced model, no imputation results
resPOOLred =  bind_cols(map(map(full_results,7),1)) %>%
  mutate(Var  = row.names(map(full_results,7)[[1]])) %>%
  select(Var, everything())



## just changing the variable name to be consistent with the non pooled models. 

resPOOL$Var[which(resPOOL$Var=='hba1c_U65:age60')] = 'age60:hba1c_U65'
resPOOL$Var[which(resPOOL$Var=='hba1c_O65:age60')] = 'age60:hba1c_O65'


resPOOLn1$Var[which(resPOOLn1$Var=='hba1c_U65:age60')] = 'age60:hba1c_U65'
resPOOLn1$Var[which(resPOOLn1$Var=='hba1c_O65:age60')] = 'age60:hba1c_O65'
 

resPOOLred$Var[which(resPOOLred$Var=='hba1c_U65:age60')] = 'age60:hba1c_U65'
resPOOLred$Var[which(resPOOLred$Var=='hba1c_O65:age60')] = 'age60:hba1c_O65'




rPool = left_join(data.frame(Var = resPOOL$Var,
                             resPOOL_Est  = apply(resPOOL %>% select(-Var),1,mean, na.rm = TRUE),
                             resPOOL_SE = apply(resPOOL %>% select(-Var),1,sd, na.rm = TRUE)),

    data.frame(Var = resPOOLn1$Var,
               resPOOLn1_Est  = apply(resPOOLn1 %>% select(-Var),1,mean, na.rm = TRUE),
               resPOOLn1_SE = apply(resPOOLn1 %>% select(-Var),1,sd, na.rm = TRUE)),
               
    by = 'Var') %>% 
  left_join(data.frame( Var = resPOOLred$Var,  
                        resPOOLred_Est = apply(resPOOLred%>% select(-Var),1,mean, na.rm = TRUE),
                        resPOOLred_SE = apply(resPOOLred%>% select(-Var),1,sd, na.rm = TRUE)), by = 'Var'
  )%>% 
  left_join(data.frame( Var = resHRS$Var,  
                        resHRS_Est = apply(resHRS %>% select(-Var),1,mean, na.rm = TRUE),
                        resHRS_SE = apply(resHRS %>% select(-Var),1,sd, na.rm = TRUE)), by = 'Var'
  )%>% 
  left_join(data.frame( Var = resARIC$Var,  
                        resARIC_Est = apply(resARIC %>% select(-Var),1,mean, na.rm = TRUE),
                        resARIC_SE = apply(resARIC %>% select(-Var),1,sd, na.rm = TRUE)), by = 'Var'
  ) %>% 
  left_join(data.frame( Var = resHRSimp$Var,  
                        resHRSimp_Est = apply(resHRSimp %>% select(-Var),1,mean, na.rm = TRUE),
                        resHRSimp_SE = apply(resHRSimp %>% select(-Var),1,sd, na.rm = TRUE)), by = 'Var'
  )%>% 
  left_join(data.frame( Var = resARICimp$Var,  
                        resARICimp_Est = apply(resARICimp %>% select(-Var),1,mean, na.rm = TRUE),
                        resARICimp_SE = apply(resARICimp %>% select(-Var),1,sd, na.rm = TRUE)), by = 'Var'
  )



write.csv(rPool, file = 'rPool_final_sim5000.csv', row.names = FALSE)

### create summary tables for participant statistics

options(scipen=999)

hrs_sim_df = do.call(rbind,map(full_results, 8))
aric_sim_df = do.call(rbind,map(full_results, 9))
hrs_sim_tab = data.frame(HRS = apply(hrs_sim_df, 2, mean))
aric_sim_tab = data.frame(ARIC= apply(aric_sim_df, 2, mean))

write.csv(hrs_sim_tab, file = 'hrs_sim_tab_sim5000.csv', row.names = FALSE)
write.csv(aric_sim_tab, file = 'aric_sim_tab_sim5000.csv', row.names = FALSE)

sink("sessionInfo_sim5000.txt")
sessionInfo()
sink()

