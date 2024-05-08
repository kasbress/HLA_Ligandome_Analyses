Test XGB models
================
Kaspar Bresser
08/01/2024

``` r
library(caret)
library(xgboost)
library(tidyverse)
```

Import the feature library

``` r
feature.table <- read_tsv("../2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip") %>% 
  mutate(across(everything(), replace_na, 0))

feature.table <- rename(feature.table, swissprot_id = Entry)
```

Import feature library model, we’ll start with those predictions, as
these are the most memory intensive.

``` r
afflib.model <- read_rds("./XGB_models_final/xgb_ligands|aff_lib|08-12-2023.RDS")
afflib.model.small <- read_rds("./XGB_models_final/xgb_ligands|aff_lib_small|08-12-2023.RDS")
afflib.model.small.rna <- read_rds("./XGB_models_final/xgb_ligands|aff_lib_small_rna|08-12-2023.RDS")
afflib.model.small.ribo <- read_rds("./XGB_models_final/xgb_ligands|aff_lib_small_ribo|08-12-2023.RDS")
afflib.model.small.rnaribo <- read_rds("./XGB_models_final/xgb_ligands|aff_lib_small_rna_ribo|08-12-2023.RDS")
```

Import the test set, remove NAs

``` r
test.set <- read_tsv("Data/HLA_test_set_complete.tsv") %>% na.omit()
```

Define a function to perform the predictions for the library. I was
using chuncked input to `predict()` before, and used `cbind()` to
aggregate the results, but annoyingly the output was not always the same
length, couldn’t find a way to fix this.

The function below is very slow, but it works. It takes a chunk of the
test table, joins it with the feature table, nests the data into a
column, and uses `map()` to output prediction results into a new column.
Finally the dataframe is formatted and returned.

``` r
predictions_lib <- function(dat){
  gc()
  print("predicting...")
  dat %>% 
    left_join(feature.table, by = "swissprot_id") %>% 
    na.omit() %>% 
    nest() %>% 
    mutate(model = list(afflib.model),
           model.small = list(afflib.model.small),
           model.small.rna = list(afflib.model.small.rna),
           model.small.ribo = list(afflib.model.small.ribo),
           model.small.rnaribo = list(afflib.model.small.rnaribo)) %>% 
    mutate(Pred = map2(model, data, predict, "prob"),
           Pred.small = map2(model.small, data, predict, "prob"),
           Pred.small.rna = map2(model.small.rna, data, predict, "prob"),
           Pred.small.ribo = map2(model.small.ribo, data, predict, "prob"),
           Pred.small.rnaribo = map2(model.small.rnaribo, data, predict, "prob")) %>% 
    mutate(data = map(data, select, ligand, Peptide, swissprot_id, tumor)) %>% 
    select(!c(model, model.small, model.small.rna, model.small.ribo,model.small.rnaribo)) %>% 
    mutate(Pred = map(Pred, transmute, aff_lib = `TRUE`),
           Pred.small = map(Pred.small, transmute, aff_lib_small = `TRUE`),
           Pred.small.rna = map(Pred.small.rna, transmute, aff_lib_small_rna = `TRUE`),
           Pred.small.ribo = map(Pred.small.ribo, transmute, aff_lib_small_ribo = `TRUE`),
           Pred.small.rnaribo = map(Pred.small.rnaribo, transmute, aff_lib_small_rna_ribo = `TRUE`)) %>% 
    unnest(c(Pred, Pred.small,Pred.small.rna, Pred.small.ribo, Pred.small.rnaribo, data)) 
}
```

Divide the test set into chunks and `map()` the function. Aggregate all
results together.

``` r
#### Library predictions
num_groups = 75

test.set %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest() %>% 
  pull(data) %>%
  map(predictions_lib) %>% 
  reduce(bind_rows) -> preds1
```

Get the remaining models. Get the file list, extract the names that will
be used as column names

``` r
files <- list.files("XGB_models_final", pattern = "xgb")
files <- files[!grepl('lib',files)]

files %>% 
  str_extract("\\|.*\\|") %>% 
  str_remove_all("\\|") -> names

files <- set_names(files, names)

files %>% 
  map(~read_rds(paste0("XGB_models_final/", .))) %>% 
  set_names(names) -> XGB.models

#rm(c(afflib.model. afflib.model.small, feature.table))
```

These predictions don’t need all the features, so they can be
immediately performed on the table in 1 go.

``` r
XGB.models %>% 
  map2_dfc(list(test.set) , ~predict(object = .x, newdata = .y, type = "prob")$`TRUE`) %>% 
  bind_cols(test.set) -> preds
```

Output data.

``` r
preds1 %>% 
  left_join(preds) %>% 
  write_tsv("Output/all_predictions_new.tsv")
```
