Generate test set for XGB models
================
Kaspar Bresser
25/03/2024

- [Finish test table](#finish-test-table)

``` r
library(tidyverse)
```

## Finish test table

Import the tables and combine into single table.

``` r
files <- list.files("./Output/netMHCpan_predictions_clean/")

files %>% 
  str_split("\\.") %>% 
  map(1) %>% 
  as.character() -> alleles

tibble(allele = alleles, 
       dat = map(files, ~read_table(paste0("./Output/netMHCpan_predictions_clean/", .)))) %>% 
  unnest(dat) -> test.table

test.table
```

    ## # A tibble: 10,010,000 × 14
    ##    allele   Pos MHC   Peptide Core     Of    Gp    Gl    Ip    Il Icore Identity
    ##    <chr>  <dbl> <chr> <chr>   <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <chr> <chr>   
    ##  1 A0201      1 HLA-… LLLPPP… LLLP…     0     7     1     0     0 LLLP… PEPLIST 
    ##  2 A0201      1 HLA-… ALYEYQ… ALYE…     0     4     1     0     0 ALYE… PEPLIST 
    ##  3 A0201      1 HLA-… LLYDFQ… LLYD…     0     4     1     0     0 LLYD… PEPLIST 
    ##  4 A0201      1 HLA-… SLFAGQ… SLFG…     0     3     1     0     0 SLFA… PEPLIST 
    ##  5 A0201      1 HLA-… YMLPDG… YMLP…     0     8     1     0     0 YMLP… PEPLIST 
    ##  6 A0201      1 HLA-… VLFDVS… VLFD…     0     4     1     0     0 VLFD… PEPLIST 
    ##  7 A0201      1 HLA-… KLLDEA… KLLD…     0     5     1     0     0 KLLD… PEPLIST 
    ##  8 A0201      1 HLA-… KIWDVS… KIWD…     0     7     1     0     0 KIWD… PEPLIST 
    ##  9 A0201      1 HLA-… YIMDNK… YIMD…     0     4     1     0     0 YIMD… PEPLIST 
    ## 10 A0201      1 HLA-… LLYEDI… LLYE…     0     7     1     0     0 LLYE… PEPLIST 
    ## # ℹ 10,009,990 more rows
    ## # ℹ 2 more variables: Score_EL <dbl>, `%Rank_EL` <dbl>

Add to test table

``` r
test.table %>% 
  transmute(sequence = Peptide, allele = allele, rank = `%Rank_EL`) %>% 
  nest(dat = -allele) -> test.table

read_tsv("./Output/Broad_Test_Table.tsv") %>% 
  nest(dat2 = -allele) %>% 
  inner_join(test.table) %>% 
  mutate(dat = map2(dat, dat2, ~inner_join(.x, .y, by = c("sequence")))) %>% 
  select(-dat2) %>% 
  unnest(dat) -> test.table

test.table %>% 
  count(allele, ligand)
```

    ## # A tibble: 40 × 3
    ##    allele ligand      n
    ##    <chr>  <lgl>   <int>
    ##  1 A0201  FALSE  499985
    ##  2 A0201  TRUE      500
    ##  3 A1101  FALSE  499989
    ##  4 A1101  TRUE      500
    ##  5 A2301  FALSE  499988
    ##  6 A2301  TRUE      500
    ##  7 A2402  FALSE  499991
    ##  8 A2402  TRUE      500
    ##  9 A3101  FALSE  499987
    ## 10 A3101  TRUE      500
    ## # ℹ 30 more rows

``` r
write_tsv(test.table, "Output/Test_Table_Broad.tsv")
```
