Generate peptide database
================
Kaspar Bresser
22/03/2024

- [Import ligands](#import-ligands)
- [Select decoys](#select-decoys)
- [Run netMHCpan](#run-netmhcpan)
- [Run netMHCpanExp](#run-netmhcpanexp)

``` r
library(babelgene)
library(tidyverse)
library(readxl)
library(furrr)
```

## Import ligands

``` r
excel_sheets("Data/Abelin_Immunity_data.xlsx")[6:19]
```

    ##  [1] "A0203" "A0204" "A0207" "A0301" "A2402" "A2902" "A3101" "A6802" "B3501"
    ## [10] "B4402" "B4403" "B5101" "B5401" "B5701"

``` r
alleles <- excel_sheets("Data/Abelin_Immunity_data.xlsx")[6:19]

excel.file <- "Data/Abelin_Immunity_data.xlsx"

tibble(data = map(alleles, ~read_excel(excel.file, sheet = . )), allele = alleles) %>% 
  mutate(data = map(data, ~select(., sequence, entry_name))) %>% 
  unnest(data) %>% 
  rename(UCSC.id = entry_name) %>% 
  mutate(UCSC.id = str_remove(UCSC.id, "\\..+"),
         sequence = toupper(sequence)) -> broad.HLA.data
  
broad.HLA.data
```

    ## # A tibble: 23,492 × 3
    ##    sequence     UCSC.id  allele
    ##    <chr>        <chr>    <chr> 
    ##  1 KLYDIDVAKV   uc002hci A0203 
    ##  2 VLKDLVKVI    uc003ymt A0203 
    ##  3 VLYPHEPTAVF  uc002ysi A0203 
    ##  4 KMYEEFLSKV   uc002jhg A0203 
    ##  5 GLADLERAQAMI uc002amp A0203 
    ##  6 ILMEHIHKL    uc002hrq A0203 
    ##  7 ILMEHIHKL    uc002hrq A0203 
    ##  8 RVLDFDPMAV   uc002anr A0203 
    ##  9 ALFPHLLQPVLW uc002xty A0203 
    ## 10 KMYEEFLSKV   uc002jhg A0203 
    ## # ℹ 23,482 more rows

Import master sheet

``` r
master.sheets <- read_tsv("Data/Abelin_immunity_data.tsv")
```

Prep master sheet

``` r
master.sheets %>% 
  pivot_longer(cols = everything(), names_to = "allele", values_to = "sequence", values_drop_na = T) %>% 
  mutate(sequence = toupper(sequence),
         allele = str_remove(allele, ":")) -> master.sheets
```

Get conversions for IDs to swissprot

``` r
conversion.table <- read_tsv("Data/swissprot_conversion_2023_12_21.tsv.gz")

conversion.table %>% 
  transmute(swissprot.id = Entry, UCSC.id = str_remove(UCSC, "\\..+")) -> conversion.table
```

Add swissprot IDs and only retain peptides in master sheet

``` r
broad.HLA.data %>% 
  inner_join(conversion.table) %>% 
  semi_join(master.sheets, by = c("allele", "sequence")) %>% 
  distinct(sequence, allele, .keep_all = T) -> broad.HLA.data

broad.HLA.data
```

    ## # A tibble: 13,783 × 4
    ##    sequence         UCSC.id  allele swissprot.id
    ##    <chr>            <chr>    <chr>  <chr>       
    ##  1 KLYDIDVAKV       uc002hci A0203  P62750      
    ##  2 KMYEEFLSKV       uc002jhg A0203  P10644      
    ##  3 GLADLERAQAMI     uc002amp A0203  Q15751      
    ##  4 ILMEHIHKL        uc002hrq A0203  P84098      
    ##  5 ALFPHLLQPVLW     uc002xty A0203  P55060      
    ##  6 KVKDILSKV        uc001zkb A0203  Q52LJ0      
    ##  7 FTFPNRLLTT       uc001bza A0203  Q9UBB6      
    ##  8 SLKEMVSKL        uc001arc A0203  P52209      
    ##  9 SVISHLLRV        uc003eib A0203  O95219      
    ## 10 AFGGSGGRGSSSGGGY uc001sat A0203  P35908      
    ## # ℹ 13,773 more rows

Subsample the peptide pools to get a more workable number, 500 is
slightly under the smallest allele-set.

``` r
"../2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip" %>% 
  read_tsv() %>% 
  pull(Entry) -> IDs.in.library

broad.HLA.data %>% 
  filter(nchar(sequence) %in% 9:11) %>%  
  filter(swissprot.id %in% IDs.in.library) %>% 
  group_by(allele) %>% 
  slice_sample(n = 500) %>% 
  mutate(ligand = TRUE) %>% 
  ungroup() -> filtered.peptides

filtered.peptides
```

    ## # A tibble: 6,971 × 5
    ##    sequence   UCSC.id  allele swissprot.id ligand
    ##    <chr>      <chr>    <chr>  <chr>        <lgl> 
    ##  1 FLHARLRTA  uc002htn A0203  O43242       TRUE  
    ##  2 YLGIVELLV  uc001wtf A0203  P25963       TRUE  
    ##  3 LLAEKVEQL  uc003wwt A0203  Q13454       TRUE  
    ##  4 YLTEVFLHV  uc003hvx A0203  Q8NDB2       TRUE  
    ##  5 RLAALGRQV  uc002mhi A0203  Q8IUR0       TRUE  
    ##  6 ALAEIAKAEL uc001bys A0203  P23246       TRUE  
    ##  7 LLAELPASV  uc002hdc A0203  Q14254       TRUE  
    ##  8 ILAQQPLSV  uc002azp A0203  Q6ZRI6       TRUE  
    ##  9 YIYDKDMEII uc002gym A0203  Q9UPT9       TRUE  
    ## 10 FMASQMLKV  uc002mgx A0203  Q9HCS7       TRUE  
    ## # ℹ 6,961 more rows

``` r
write_tsv(filtered.peptides, "Output/Abelin/Broad_peptides_subsample.tsv")
```

## Select decoys

Import UniProt sequences, use the rna expression table to make an
ensembl/UniProt matching table.

Define a function that, for each allele: (1) Checks which sequences were
detected for that allele, (2) samples an excess of proteins from the
uniprot table weighted by their length, (3) then for each sample
extracts a 9mer, (4) Remove detected peptides, filter for uniqueness and
down sample to 350,000.

``` r
uniprot <- read_csv("Data/UniProt_reviewed_input.tsv") %>% rename(swissprot_id = "sequence_id")


get_peptide <- function(a, pep.len){
  broad.HLA.data %>% 
    filter(allele == a) %>% 
    filter(nchar(sequence) == pep.len) %>% 
    pull(sequence) %>% 
    toupper() -> s
  
    amount.needed <- nrow(filter(filtered.peptides, allele == a & nchar(sequence) == pep.len))*10
  

    uniprot %>% 
      mutate(len = nchar(sequence)) %>% 
      slice_sample(n = amount.needed*2, weight_by = len, replace = T) %>% 
      rowwise() %>% 
      mutate(number = sample(1:(nchar(sequence)-pep.len), 1),
           sequence = str_sub(sequence, number, number+(pep.len-1))) %>% 
      ungroup() %>% 
      filter(!(sequence %in% s)) %>% 
      transmute(sequence = sequence, swissprot.id = swissprot_id, allele = a, ligand = FALSE) %>% 
      filter(swissprot.id %in% IDs.in.library) %>% 
      distinct(sequence, .keep_all = T) %>% 
      slice_sample(n = amount.needed)
}

decoys <- map2_dfr(rep(alleles, 3), rep(c(9,10,11), each = length(alleles)), get_peptide)

filtered.peptides %>% 
  select(!UCSC.id) %>% 
  bind_rows(decoys) -> test.set


test.set %>% 
  count(allele, ligand, nchar(sequence))
```

Write out the table and peptide files

## Run netMHCpan

Run system command for netMHCpan, for each allele.

Define function to read in predictions for each peptide length, cleanup
the netMHC output, and combine them in a single tsv file

``` r
combine_and_clean <- function(allele){

  file.list <- list.files("./Output/Abelin/netMHCpan_predictions/", pattern = allele)
  file.list <- paste0("./Output/Abelin/netMHCpan_predictions/", file.list)
  
  file.list %>% 
    map(read_lines) %>% 
    map(~.[grepl("   1 HLA", .)]) %>% 
    unlist() %>% 
    c(" Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL", .) %>% 
    write_lines(paste0("./Output/Abelin/netMHCpan_predictions_clean/", allele,".tsv"))
  
  gc()
}
```

Apply the function for each allele

``` r
map(alleles, combine_and_clean)
```

## Run netMHCpanExp

Run system command for netMHCpan, for each allele.

Define function to read in predictions for each peptide length, cleanup
the netMHC output, and combine them in a single tsv file

``` r
combine_and_clean <- function(allele){

  file.list <- list.files("./Output/Abelin/netMHCpanExp_predictions/", pattern = allele)
  file.list <- paste0("./Output/Abelin/netMHCpanExp_predictions/", file.list)
  
  file.list %>% 
    map(read_lines) %>% 
    map(~.[grepl("   1 HLA", .)]) %>% 
    unlist() %>% 
    c(" Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL", .) %>% 
    write_lines(paste0("./Output/Abelin/netMHCpanExp_predictions_clean/", allele,".tsv"))
  
  gc()
}
```

Apply the function for each allele

``` r
map(alleles, combine_and_clean)
```

Import the tables and combine into single table.

``` r
files <- list.files("./Output/Abelin/netMHCpan_predictions_clean/")

alleles <- str_split_i(files, "\\.", 1)

tibble(allele = alleles, 
       dat = map(files, ~read_table(paste0("./Output/Abelin/netMHCpan_predictions_clean/", .)))) %>% 
  unnest(dat) -> train.table.NetMHC

train.table.NetMHC
```

``` r
files <- list.files("./Output/Abelin/netMHCpanExp_predictions_clean/")

alleles <- str_split_i(files, "\\.", 1)

tibble(allele = alleles, 
       dat = map(files, ~read_table(paste0("./Output/Abelin/netMHCpanExp_predictions_clean/", .)))) %>% 
  unnest(dat) -> train.table.NetMHCexp

train.table.NetMHCexp
```

Add to test table

``` r
train.table.NetMHC %>% 
  transmute(sequence = Peptide, allele = allele, rank = `%Rank_EL`) %>% 
  nest(netMHC = -allele) -> train.table.NetMHC

train.table.NetMHCexp %>% 
  transmute(sequence = Peptide, allele = allele, rankExp = `%Rank_EL`) %>% 
  nest(netMHCexp = -allele) -> train.table.NetMHCexp

read_tsv("./Output/Abelin/Broad_Train_Table.tsv") %>% 
  nest(train = -allele) %>% 
  inner_join(train.table.NetMHC) %>% 
  inner_join(train.table.NetMHCexp) %>% 
  mutate(train = map2(train, netMHC, ~inner_join(.x, .y, by = c("sequence")))) %>% 
  mutate(train = map2(train, netMHCexp, ~inner_join(.x, .y, by = c("sequence")))) %>% 
  select(-netMHC, -netMHCexp) %>% 
  unnest(train) %>% 
  distinct(allele, sequence, .keep_all = T) %>% 
  na.omit() -> train.table

train.table %>% 
  count(allele, ligand)
```

``` r
write_tsv(train.table, "Output/Abelin/Broad_Train_Table_All.tsv")
```
