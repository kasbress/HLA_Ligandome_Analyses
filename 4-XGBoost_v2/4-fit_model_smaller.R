# title: "XGBoost for ligands trial"
# author: K.Bresser & B. Nicolet
# date: "14/03/2024"

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

## set work directory
setwd("~/PATH/TO/DIR")


## Import melanoma data
read_tsv("/Output/Melanoma/Melanoma_Train_Table_all.tsv") %>% 
  transmute(swissprot.id = swissprot_id, ligand = ligand, rank = rank) %>% 
  group_by(ligand) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map2(data, c(120000,12000), sample_n)) %>% 
  unnest(data) -> train.data

## Import mono-allelic data and bind rows
read_tsv("7-New_Train_set/Output/Abelin/Broad_Train_Table_All.tsv") %>% 
  select(rank, swissprot.id, ligand) %>% 
  bind_rows(train.data) -> train.data

## Import Lung adenocarcinoma data and bind rows
read_tsv("7-New_Train_set/Output/Lung/Lung_Train_Table_all.tsv") %>% 
  select(rank, swissprot.id, ligand) %>% 
  bind_rows(train.data) -> train.data



## Import the v1 sequence feature model
lib.model <- read_rds("3-XGBoost-analyses/XGB_models_final/xgb_ligands|aff_lib|08-12-2023.RDS")

## Extract important features
importances <- varImp(lib.model)
importances$importance %>% 
  as_tibble(rownames = "feature") %>% 
  filter(Overall > 0 ) %>% 
  pull(feature) -> features

## Filter the sequence feature library, to only keep important features
read_tsv("2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip") %>% 
  mutate(across(everything(), replace_na, 0)) %>% 
  select(Entry, any_of(features)) %>% 
  inner_join(train.data, by = c("Entry" = "swissprot.id")) %>% 
  na.omit() -> Full_data


## Remove stuff that's not needed
Full_data$swissprot.id <- NULL
Full_data$Entry <- NULL
rm(importances)
rm(lib.model)

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

### Modeling ###

###---------###
## Set up ##
###---------###
gc()
registerDoMC(6)
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



###-----------------------###
## Aff + lib
###-----------------------###
print("Aff + lib model")

start_time <- Sys.time()

Full_data %>% 
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


## Output the sequence feature model 2.0
saveRDS(XGB.model,"7-New_Train_set/Output/xgb_ligands|aff_lib_2|08-12-2023.RDS")


