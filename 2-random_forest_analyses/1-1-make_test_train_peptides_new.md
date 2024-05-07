Make test and train sets
================
Kaspar Bresser
04/10/2023


Used the code below to generate the train and test sets used for the
random forests.

``` r
library(tidyverse)
library(furrr)
```


Read the MS detected peptides into a nested tibble, then unnest and
select the columns we want to continue with. Also, only keep unique
cases.

``` r
file.names <- list.files(("Data/MS_peptides")) 

HLA.ligands <- tibble(
  
  ligand.table = map( paste0("Data/MS_peptides/", file.names), read_tsv ),
  tumor = map_chr(str_split(file.names, "_"), 1)

  )

#HLA.ligands <- mutate(HLA.ligands, case_when(tumor))

HLA.ligands
```

    ## # A tibble: 3 × 2
    ##   ligand.table            tumor
    ##   <list>                  <chr>
    ## 1 <spc_tbl_ [4,606 × 37]> M026 
    ## 2 <spc_tbl_ [6,024 × 30]> mel95
    ## 3 <spc_tbl_ [7,198 × 28]> RTE

``` r
HLA.ligands <- mutate(HLA.ligands, tumor = c("M026X1", "SKMEL95", "NKIRTIL006"))

# Define Columns to keep from the ligand tables
keep <- c("rna","ribo",  "swissprot_id", "ligand", "peptide" )

HLA.ligands %>% 
  unnest( cols = ligand.table) %>% 
  select(one_of(keep, "tumor")) %>% 
  distinct() -> HLA.ligands

HLA.ligands
```

    ## # A tibble: 16,619 × 6
    ##      rna   ribo swissprot_id ligand peptide   tumor 
    ##    <dbl>  <dbl> <chr>        <lgl>  <chr>     <chr> 
    ##  1  29.5  2556. Q02952       TRUE   AAAEEEKVL M026X1
    ##  2 468.    962. O00567       TRUE   AAAEITRKL M026X1
    ##  3  11.8  1044. Q9P2E5       TRUE   AAAELERRY M026X1
    ##  4 134.    840. Q13425       TRUE   AAAELIKEV M026X1
    ##  5  28.6   562. Q9C0D5       TRUE   AAAGHMKLV M026X1
    ##  6  44.5   575. P54802       TRUE   AAAGLHRYL M026X1
    ##  7  15.1   202  Q9NRK6       TRUE   AAANAIRVY M026X1
    ##  8 583.  17890. P24821       TRUE   AAAPDVKEL M026X1
    ##  9  13.1    82  Q8IW45       TRUE   AAAPVIKAY M026X1
    ## 10  25.9   136. Q68EM7       TRUE   AAATSVHVV M026X1
    ## # ℹ 16,609 more rows

Split test, train. 80/20

``` r
HLA.ligands %>% 
  group_by(tumor) %>% 
  sample_frac(.8) -> train.ligands

HLA.ligands %>% 
  anti_join(train.ligands) -> test.ligands
```

import expression info

``` r
file.names <- list.files("Data/protein_data/") 

expression.data <- tibble(
  
  expression.info = map( paste0("Data/protein_data/", file.names), read_tsv ),
  tumor = map_chr(str_split(file.names, "_"), 1)

  )

(expression.data <- unnest(expression.data, cols = expression.info) %>% select(!length))
```

    ## # A tibble: 47,958 × 4
    ##    swissprot_id     rna   ribo tumor 
    ##    <chr>          <dbl>  <dbl> <chr> 
    ##  1 Q8IW70        0.0377   0    M026X1
    ##  2 Q8NAP8        0.0614   2.28 M026X1
    ##  3 Q8IV76        7.78    84    M026X1
    ##  4 O94876       21.8     19    M026X1
    ##  5 Q6P1R3        1.59    13.0  M026X1
    ##  6 Q15714       19.6    357.   M026X1
    ##  7 Q13207        8.92   544.   M026X1
    ##  8 Q9H4H8       29.0    133    M026X1
    ##  9 P31277        0.0632   0    M026X1
    ## 10 Q96KQ7       21.1    305.   M026X1
    ## # ℹ 47,948 more rows

import swissprotIDs that are in the library, import swissprot sequences,
and combine with expression data

``` r
"../2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv" %>% 
  read_tsv() %>% 
  pull(Entry) -> IDs.in.library

uniprot <- read_csv("UniProt_reviewed_input.tsv") %>% rename(swissprot_id = "sequence_id")

uniprot %>% 
  filter(swissprot_id %in% IDs.in.library)  -> uniprot

expression.data %>% 
  nest(expr = !tumor) %>% 
  mutate(expr = map(expr, ~right_join(., uniprot))) %>% 
  unnest(expr) %>% 
  replace_na(list(rna = 0, ribo = 0)) -> expression.data
```

sample peptides for train data.

``` r
get_peptide <- function(tum, pep.len){
  train.ligands %>% 
    filter(tumor == tum) %>% 
    filter(nchar(peptide) == pep.len) %>% 
    pull(peptide) %>% 
    toupper() -> s
  
  num <- nrow(filter(train.ligands, tumor == tum & nchar(peptide) == pep.len))*10
  
  expression.data %>% 
    filter(tumor == tum) %>% 
    mutate(peptide = sequence) %>% 
    mutate(len = nchar(peptide)) %>% 
    slice_sample(n = num*2, weight_by = len, replace = T) %>% 
    rowwise() %>% 
    mutate(num = sample(1:(nchar(peptide)-pep.len), 1),
         peptide = str_sub(peptide, num, num+(pep.len-1))) %>% 
    ungroup() %>% 
    filter(!(peptide %in% s)) %>% 
    select(peptide, swissprot_id, tumor, rna, ribo) %>% 
    mutate(ligand = FALSE) %>% 
    filter(swissprot_id %in% IDs.in.library) %>% 
    distinct(peptide, .keep_all = T) %>% 
    slice_sample(n = num)
}



decoys <- map2_dfr(rep(c("M026X1", "SKMEL95", "NKIRTIL006"), 3), rep(c(9,10,11), each = 3), get_peptide)

train.set <- bind_rows(decoys, train.ligands)

train.set %>% 
  count(tumor, ligand, nchar(peptide))
```

    ## # A tibble: 18 × 4
    ##    tumor      ligand `nchar(peptide)`     n
    ##    <chr>      <lgl>             <int> <int>
    ##  1 M026X1     FALSE                 9 29120
    ##  2 M026X1     FALSE                10  3180
    ##  3 M026X1     FALSE                11  1920
    ##  4 M026X1     TRUE                  9  2912
    ##  5 M026X1     TRUE                 10   318
    ##  6 M026X1     TRUE                 11   192
    ##  7 NKIRTIL006 FALSE                 9 40770
    ##  8 NKIRTIL006 FALSE                10  9780
    ##  9 NKIRTIL006 FALSE                11  4110
    ## 10 NKIRTIL006 TRUE                  9  4077
    ## 11 NKIRTIL006 TRUE                 10   978
    ## 12 NKIRTIL006 TRUE                 11   411
    ## 13 SKMEL95    FALSE                 9 29090
    ## 14 SKMEL95    FALSE                10  9340
    ## 15 SKMEL95    FALSE                11  5640
    ## 16 SKMEL95    TRUE                  9  2909
    ## 17 SKMEL95    TRUE                 10   934
    ## 18 SKMEL95    TRUE                 11   564

sample peptides for test data.

``` r
get_peptide <- function(tum, pep.len){
  test.ligands %>% 
    filter(tumor == tum) %>% 
    filter(nchar(peptide) == pep.len) %>% 
    pull(peptide) %>% 
    toupper() -> s
  
  num <- nrow(filter(test.ligands, tumor == tum & nchar(peptide) == pep.len))*1000
  
  expression.data %>% 
    filter(tumor == tum) %>% 
    mutate(peptide = sequence) %>% 
    mutate(len = nchar(peptide)) %>% 
    slice_sample(n = num*4, weight_by = len, replace = T) %>% 
    rowwise() %>% 
    mutate(num = sample(1:(nchar(peptide)-pep.len), 1),
         peptide = str_sub(peptide, num, num+(pep.len-1))) %>% 
    ungroup() %>% 
    filter(!(peptide %in% s)) %>% 
    select(peptide, swissprot_id, tumor, rna, ribo) %>% 
    mutate(ligand = FALSE) %>% 
    filter(swissprot_id %in% IDs.in.library) %>% 
    distinct(peptide, .keep_all = T) %>% 
    slice_sample(n = num)
}


decoys <- map2_dfr(rep(c("M026X1", "SKMEL95", "NKIRTIL006"), 3), rep(c(9,10,11), each = 3), get_peptide)

test.set <- bind_rows(decoys, test.ligands)

test.set %>% 
  count(tumor, ligand, nchar(peptide))
```

    ## # A tibble: 18 × 4
    ##    tumor      ligand `nchar(peptide)`       n
    ##    <chr>      <lgl>             <int>   <int>
    ##  1 M026X1     FALSE                 9  719000
    ##  2 M026X1     FALSE                10   75000
    ##  3 M026X1     FALSE                11   61000
    ##  4 M026X1     TRUE                  9     719
    ##  5 M026X1     TRUE                 10      75
    ##  6 M026X1     TRUE                 11      61
    ##  7 NKIRTIL006 FALSE                 9 1021000
    ##  8 NKIRTIL006 FALSE                10  250000
    ##  9 NKIRTIL006 FALSE                11   96000
    ## 10 NKIRTIL006 TRUE                  9    1021
    ## 11 NKIRTIL006 TRUE                 10     250
    ## 12 NKIRTIL006 TRUE                 11      96
    ## 13 SKMEL95    FALSE                 9  700000
    ## 14 SKMEL95    FALSE                10  237000
    ## 15 SKMEL95    FALSE                11  165000
    ## 16 SKMEL95    TRUE                  9     700
    ## 17 SKMEL95    TRUE                 10     237
    ## 18 SKMEL95    TRUE                 11     165

Write out the table and peptide files

``` r
write_tsv(train.set, "./Output/test_train_sets/HLA_train_set.tsv")

train.set %>% 
  transmute(pep.len = paste0(nchar(peptide), "AA"), peptide = peptide) %>% 
  nest(data = peptide) %>% 
  mutate(peptides = map(data, pull, peptide),
         peptides = set_names(peptides, pep.len)) %>% 
  pull(peptides) %>% 
  map2(names(.), ~write_lines(.x, paste0("./Output/test_train_sets/Broad_Train_Peptides_", .y, ".tsv")))
```
