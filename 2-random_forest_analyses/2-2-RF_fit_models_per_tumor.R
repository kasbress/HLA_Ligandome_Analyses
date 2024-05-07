library(caTools)
library(AUC)
library(caret)
library(randomForest)
library(doMC)
library(e1071)
library(tidyverse)
set.seed(1)

# This script is to obtain a feature ranking for each of the 4 subsets of
# features. the features are stored as separate files, one for each subset. I
# will use a partition of the train data as an input of the feature ranking.

# Import data -------------------------------------------------------------

# Set working directory
setwd("PATH/TO/DIR/2-random_forest_analyses")

# Read in train set
train.peptides <- read_tsv("./Output/test_train_sets/HLA_train_set_complete.tsv")

# Drop columns not needed for training and switch ligand column to factor
train.peptides %>% 
  mutate(ligand = factor(case_when(ligand == TRUE ~ "yes", TRUE ~ "no"), 
                         levels = c("yes", "no") )) %>% 
  dplyr::select(tumor, ligand, swissprot_id) -> train.peptides

# Sub-sample train set
train.peptides %>% 
  group_by(tumor, ligand) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(sample.size = rep(c(8000, 2000), 3) ) %>% 
  mutate(data = map2(data, sample.size, sample_n)) %>% 
  select(-sample.size) %>% 
  unnest(data) -> train.peptides


# Import features tables
paths <- paste0("./Output/feature_files/", list.files("./Output/feature_files/"))
#paths <- paste0("./Output/feature_files/count_features.tsv")

list.files("./Output/feature_files/") %>% 
  str_split("_") %>% 
  map_chr(1) -> table.names

#table.names <- "count"

paths %>% 
  map(read_tsv) %>% 
  set_names(table.names) -> feature.tables


feature.tables %>% 
  enframe("classes", "values") %>% 
  mutate(values = map(values, inner_join, train.peptides, c("Entry" = "swissprot_id")),
         values = map(values, group_by, tumor),
         values = map(values, nest)) %>% 
  unnest(values) -> features.per.tumor




# random forest models ----------------------------------------------------


# Set up two class summary function
twoClassSummaryCustom <- function(data, lev = NULL, model = NULL) 
{
  lvls <- levels(data$obs)
  if(length(lvls) > 2)
    stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
  if(!all(levels(data[, "pred"]) == lvls)) 
    stop("levels of observed and predicted data do not match")
  rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0, 
                                     1), data[, lvls[1]])
  #  pred <- factor(ifelse(data[, "no"] > 0.5, "no", "yes"))
  #  pred <- factor(ifelse(data[, "yes"] > 0.99, "yes", "no"))
  
  out <- c(rocAUC,
           sensitivity(data[, "pred"], data[, "obs"], lev[1]),
           specificity(data[, "pred"], data[, "obs"], lev[2]),
           posPredValue(data[, "pred"], data[, "obs"], lev[1]))
  names(out) <- c("ROC", "Sens","Spec", "Prec")
  out
}


# Register cores for parallel processing
registerDoMC(6)

# Set the train control
trctrl <-  trainControl(method = "repeatedcv", number = 6, repeats = 1, 
                        summaryFunction = twoClassSummaryCustom, classProbs = TRUE, 
                        returnData = FALSE, sampling = "down", allowParallel = T)

# Define function to build formula from column names
get_formulas <- function(df){
  df %>% 
    select(-c(Entry, ligand)) %>% 
    names() %>% 
    paste( collapse = " + ") -> feats
  
  form <- paste("ligand", feats, sep = " ~ ")
  form <- as.formula(form)
  form
}


# Get the formula's and add to tibble
features.per.tumor %>% 
  mutate(formulas = map(data, get_formulas)) -> features.per.tumor
  


features.per.tumor %>% 
  mutate(rf.models = map2(formulas, 
                          data, 
                          ~train(.x , data = .y, 
                                 method = "rf",
                                 maximize = TRUE,
                                 trControl=trctrl,
                                 tuneGrid= data.frame(mtry = c(round(sqrt(ncol(.y)-1))))*1.5,
                                 importance = T, 
                                 nodesize = 2, 
                                 ntree = 5000,
                                 metric = "ROC"))) %>% 
  select(tumor, classes, rf.models) -> rf.models



write_rds(rf.models, "./Output/RF_per_tumor_new.rds")

