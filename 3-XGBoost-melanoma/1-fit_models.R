# title: "XGBoost for ligands trial"
# author: "BP Nicolet"
# date: "14/12/2021"

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(reshape)
library(circlize)
library(stringr)
library(doMC)
library(caret)
library(randomForest)
library(e1071)
library(xgboost)
library(Matrix)

set.seed(12345)

setwd("~/Dropbox/Dropbox_Kaspar/Ligandome_Analyses/HLA_ligandome_project_revision/")

read_tsv("3-XGBoost-analyses/Data/HLA_train_set_complete.tsv") %>% 
  count(tumor, ligand)

read_tsv("3-XGBoost-analyses/Data/HLA_train_set_complete.tsv") %>% 
  na.omit() %>% 
  count(tumor, ligand)

### importing Full_data ###


read_tsv("3-XGBoost-analyses/Data/HLA_train_set_complete.tsv") %>% 
  na.omit() %>% 
  group_by(ligand) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map2(data, c(40000,10000), sample_n)) %>% 
  unnest(data) %>% 
  select(ligand, rna, ribo, rank, swissprot_id) -> Full_data

read_tsv("2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip") %>% 
  mutate(across(everything(), replace_na, 0)) %>% 
  inner_join(Full_data, by = c("Entry" = "swissprot_id")) %>% 
  na.omit() -> Full_data

Full_data <- select(Full_data, !contains("miR"))



Full_data$swissprot_id <- NULL
Full_data$Entry <- NULL



###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

### Modeling ###

###---------###
## Set up ##
###---------###
gc()
registerDoMC(4)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        search="grid",
                        verboseIter = TRUE,
                        sampling = "down",
                        allowParallel = T, 
                        returnData = F)

xgbGrid <- expand.grid(nrounds = 1000,
                       max_depth = 1,
                       colsample_bytree = 0.5,
                       eta = 0.3,
                       gamma=1,
                       min_child_weight = 0.9,
                       subsample = 1)

set.seed(12345)


#objective = "binary:logistic"



###-----------------------###
## Aff only model ##
###-----------------------###
print("Aff model")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, rank)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_only|08-12-2023.RDS")



rm(XGB.model)
gc()


###-----------------------###
## Aff + lib
###-----------------------###
print("Aff + lib model")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(-c( rna, ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib|08-12-2023.RDS")






###-----------------------###
## Aff + lib SMALL 
###-----------------------###

XGB.model <- readRDS("3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib|08-12-2023.RDS")


importances <- varImp(XGB.model)

rm(XGB.model)

gc()

importances$importance %>% 
  as_tibble(rownames = "feature") %>% 
  filter(Overall > 0 ) %>% 
  pull(feature) -> features

Full_data <- select(Full_data, all_of(c("ligand", "rank", "rna","ribo", features)))

print("Aff + lib model")

start_time <- Sys.time()

Full_data %>% 
  select( all_of(c("ligand", "rank",  features)))
train(as.factor(ligand)~. ,
      data=.,
      method="xgbTree",
      trControl=control,
      metric="Accuracy",
      tuneGrid= xgbGrid,
      na.action = na.omit,
      nthread = 12,
      verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib_small|08-12-2023.RDS")









###-----------------------###
## Aff + lib SMALL + RNA
###-----------------------###

print("Aff + libsmall + rna model")

start_time <- Sys.time()

Full_data %>% 
  select(all_of(c("ligand", "rank","rna",  features))) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib_small_rna|08-12-2023.RDS")




###-----------------------###
## Aff + lib SMALL + RIBO
###-----------------------###


rm(XGB.model)

gc()

print("Aff + lib + ribo model")

start_time <- Sys.time()

Full_data %>% 
  select( all_of(c("ligand", "rank", "ribo", features))) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib_small_ribo|08-12-2023.RDS")




###-----------------------###
## Aff + lib SMALL + RNA + RIBO
###-----------------------###


rm(XGB.model)

gc()

print("Aff + lib + ribo + RNA model")

start_time <- Sys.time()

Full_data %>% 
  select( all_of(c("ligand", "rank", "rna", "ribo", features))) %>% 
train(as.factor(ligand)~. ,
      data=.,
      method="xgbTree",
      trControl=control,
      metric="Accuracy",
      tuneGrid= xgbGrid,
      na.action = na.omit,
      nthread = 12,
      verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib_small_rna_ribo|08-12-2023.RDS")





###-----------------------###
## Aff + RNA 
###-----------------------###


rm(XGB.model)

gc()

print("Aff + RNA model")

start_time <- Sys.time()

Full_data %>% 
  select( all_of(c("ligand", "rank", "rna"))) %>% 
train(as.factor(ligand)~. ,
      data=.,
      method="xgbTree",
      trControl=control,
      metric="Accuracy",
      tuneGrid= xgbGrid,
      na.action = na.omit,
      nthread = 12,
      verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_rna|08-12-2023.RDS")



###-----------------------###
## Aff + Ribo 
###-----------------------###


rm(XGB.model)

gc()

print("Aff + Ribo model")

start_time <- Sys.time()

Full_data %>% 
  select(all_of(c("ligand", "rank", "ribo"))) %>% 
train(as.factor(ligand)~. ,
      data=.,
      method="xgbTree",
      trControl=control,
      metric="Accuracy",
      tuneGrid= xgbGrid,
      na.action = na.omit,
      nthread = 12,
      verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_ribo|08-12-2023.RDS")


###-----------------------###
## Aff + RNA + Ribo
###-----------------------###


rm(XGB.model)

gc()

print("Aff + RNA + Ribo model")

start_time <- Sys.time()

Full_data %>% 
  select( all_of(c("ligand", "rank", "rna", "ribo"))) %>% 
train(as.factor(ligand)~. ,
      data=.,
      method="xgbTree",
      trControl=control,
      metric="Accuracy",
      tuneGrid= xgbGrid,
      na.action = na.omit,
      nthread = 12,
      verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_rna_ribo|08-12-2023.RDS")



