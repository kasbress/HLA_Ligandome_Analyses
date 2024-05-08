Test XGB models on T B data
================
Kaspar Bresser
28/03/24

``` r
library(caret)
library(xgboost)
library(tidyverse)
```

Import the feature library

``` r
feature.table <- read_tsv("../2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip") %>% 
  mutate(across(everything(), replace_na, 0))

feature.table <- dplyr::rename(feature.table, swissprot.id = Entry)
```

Import feature library model, we’ll start with those predictions, as
these are the most memory intensive.

``` r
afflib.model <- read_rds("models/xgb_ligands|aff_lib_small|08-12-2023.RDS")

feature.table %>% 
  select( one_of(afflib.model$coefnames), swissprot.id) -> feature.table
```

Import the test set, remove NAs

``` r
test.set <- read_tsv("Output/T-B_Table_all.tsv") %>% na.omit()
```

``` r
afflib.model <- read_rds("models/xgb_ligands|aff_lib_2|08-12-2023.RDS")
#affExplib.model <- read_rds("models/xgb_ligands|affExp_lib_MORE_smaller|08-12-2023.RDS")
afflib.model.small <- read_rds("models/xgb_ligands|aff_lib_small|08-12-2023.RDS")
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
    left_join(feature.table, by = "swissprot.id") %>% 
    na.omit() %>% 
    nest() %>% 
    mutate(model = list(afflib.model),
           model.small = list(afflib.model.small)) %>% 
    mutate(Pred = map2(model, data, predict, "prob"),
           Pred.small = map2(model.small, data, predict, "prob")) %>% 
    mutate(data = map(data, select, ligand, sequence, swissprot.id, tumor)) %>% 
    select(!c(model, model.small)) %>% 
    mutate(Pred = map(Pred, transmute, aff_lib = `TRUE`),
           Pred.small = map(Pred.small, transmute, aff_lib_small = `TRUE`)) %>% 
    unnest(c(Pred , Pred.small, data)) 
}
```

Divide the test set into chunks and `map()` the function. Aggregate all
results together.

``` r
#### Library predictions
num_groups = 30

test.set %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest() %>% 
  pull(data) %>%
  map(predictions_lib) %>% 
  reduce(bind_rows) -> preds1
```

``` r
preds1 %>% 
  select(sequence, tumor, aff_lib) %>% 
  right_join(test.set, by = c("tumor", "sequence")) %>% 
  write_tsv("Output/all_predictions_Mening.tsv")
```

Get the remaining models. Get the file list, extract the names that will
be used as column names

``` r
XGB.models <- list(aff_only = read_rds("models/xgb_ligands|aff_only|08-12-2023.RDS"))
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
  write_tsv("Output/all_predictions_T-B.tsv")
```
