nuORF data prep
================
Kaspar Bresser
28/03/2024

- [Import ligands](#import-ligands)
- [Select decoys](#select-decoys)
- [Run netMHCpan](#run-netmhcpan)
- [Finish test table](#finish-test-table)

``` r
library(babelgene)
library(tidyverse)
library(furrr)
library(readxl)
```

## Import ligands

``` r
MS.data <- read_excel("nuORFs.xlsx")
```

``` r
MS.data %>% 
  count(plotType)
```

    ## # A tibble: 8 × 2
    ##   plotType            n
    ##   <chr>           <int>
    ## 1 3' Overlap dORF   311
    ## 2 3' dORF           497
    ## 3 5' Overlap uORF  1542
    ## 4 5' uORF          2050
    ## 5 Other             163
    ## 6 Out-of-Frame     1619
    ## 7 Pseudogene        341
    ## 8 lncRNA           2044

``` r
internal.ORFs <- c("3' Overlap dORF", "3' dORF", "5' Overlap uORF", "5' uORF", "Out-of-Frame")

MS.data %>% 
  filter(plotType %in% internal.ORFs) %>% 
  filter(nchar(sequence) < 15 & nchar(sequence) > 7) %>% 
  select(sequence, allele, `Mapped protein`, plotType ) %>% 
  group_by(allele) %>% 
  filter(n() > 100) %>% 
  ungroup() -> MS.data.filtered
```

Tidy up, and add peptide lengths.

``` r
MS.data.filtered %>% 
  mutate(pep.length = nchar(sequence))  -> MS.data.filtered

MS.data.filtered
```

    ## # A tibble: 3,027 × 5
    ##    sequence    allele `Mapped protein`                       plotType pep.length
    ##    <chr>       <chr>  <chr>                                  <chr>         <int>
    ##  1 ALWELSAALLR A0301  ENST00000571146.1_3_17:7160220-716028… 3' dORF          11
    ##  2 RVIQIHLKK   A0301  ENST00000393415.7_1_12:69229612-69230… Out-of-…          9
    ##  3 AVIVWSFMK   A3101  ENST00000514133.1_1_4:52895917-528959… 3' dORF           9
    ##  4 GTPKLIWQQR  A3101  ENST00000360051.7_1_2:242256941-24225… 5' uORF          10
    ##  5 ALGPPRPLR   A3101  ENST00000341068.7_1_2:112641540-11264… 5' uORF           9
    ##  6 AVLKPGAWSRK A0301  ENST00000308521.9_1_22:39436844-39436… 5' uORF          11
    ##  7 TTRPWEVTLGR A3101  ENST00000407997.3_1_22:39473245-39473… 5' uORF          11
    ##  8 AVLKPEAWSRK A0301  ENST00000407997.3_1_22:39473245-39473… 5' uORF          11
    ##  9 MLLLMLLYK   A0301  ENST00000359142.7_1_12:90049779-90049… 5' uORF           9
    ## 10 AAMPARLLGWR A3101  ENST00000622395.4_1_1:161147236-16114… 5' uORF          11
    ## # ℹ 3,017 more rows

``` r
convert <- read_tsv("conversion.txt") %>% select(`Ensembl Transcript ID`, `UniProt/SwissProt Accession`) %>% na.omit()

MS.data.filtered %>% 
  mutate(ens.transcript = str_extract(`Mapped protein`, "ENST\\d*")) %>% 
  inner_join(convert, by = c("ens.transcript" = "Ensembl Transcript ID")) %>% 
  distinct() -> MS.data.filtered
```

    ## Warning in inner_join(., convert, by = c(ens.transcript = "Ensembl Transcript ID")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 2564 of `x` matches multiple rows in `y`.
    ## ℹ Row 24375 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
MS.data.filtered
```

    ## # A tibble: 1,910 × 7
    ##    sequence       allele `Mapped protein`     plotType pep.length ens.transcript
    ##    <chr>          <chr>  <chr>                <chr>         <int> <chr>         
    ##  1 GTPKLIWQQR     A3101  ENST00000360051.7_1… 5' uORF          10 ENST000003600…
    ##  2 ALGPPRPLR      A3101  ENST00000341068.7_1… 5' uORF           9 ENST000003410…
    ##  3 AVLKPGAWSRK    A0301  ENST00000308521.9_1… 5' uORF          11 ENST000003085…
    ##  4 TTRPWEVTLGR    A3101  ENST00000407997.3_1… 5' uORF          11 ENST000004079…
    ##  5 AVLKPEAWSRK    A0301  ENST00000407997.3_1… 5' uORF          11 ENST000004079…
    ##  6 MLLLMLLYK      A0301  ENST00000359142.7_1… 5' uORF           9 ENST000003591…
    ##  7 ALFSMLTTK      A0301  ENST00000378119.8_1… 5' uORF           9 ENST000003781…
    ##  8 TSLSRPLIGK     A0301  ENST00000336557.9_1… 5' uORF          10 ENST000003365…
    ##  9 RVPFSHPPR      A3101  ENST00000528542.6_2… 5' uORF           9 ENST000005285…
    ## 10 KIKSFEKSNSWSWK A0301  ENST00000382952.7_2… 5' uORF          14 ENST000003829…
    ## # ℹ 1,900 more rows
    ## # ℹ 1 more variable: `UniProt/SwissProt Accession` <chr>

Subsample the peptide pools to get a more workable number, 1000 should
suffice

``` r
"../2-random_forest_analyses/Data/Protein_per_Uniprot_entry_library_v3.csv.zip" %>% 
  read_tsv() %>% 
  pull(Entry) -> IDs.in.library



MS.data.filtered %>% 
  rename(swissprot.id = `UniProt/SwissProt Accession`) %>% 
  filter(swissprot.id %in% IDs.in.library) %>% 
  mutate(ligand = TRUE) %>% 
  ungroup() %>% 
  select(sequence, allele, swissprot.id, ligand) -> filtered.peptides

filtered.peptides
```

    ## # A tibble: 1,892 × 4
    ##    sequence       allele swissprot.id ligand
    ##    <chr>          <chr>  <chr>        <lgl> 
    ##  1 GTPKLIWQQR     A3101  Q15019       TRUE  
    ##  2 ALGPPRPLR      A3101  Q9H1A4       TRUE  
    ##  3 AVLKPGAWSRK    A0301  Q8IUX4       TRUE  
    ##  4 TTRPWEVTLGR    A3101  Q9HC16       TRUE  
    ##  5 AVLKPEAWSRK    A0301  Q9HC16       TRUE  
    ##  6 MLLLMLLYK      A0301  P20020       TRUE  
    ##  7 ALFSMLTTK      A0301  Q9GZU0       TRUE  
    ##  8 TSLSRPLIGK     A0301  Q9UHD4       TRUE  
    ##  9 RVPFSHPPR      A3101  Q8N4Y2       TRUE  
    ## 10 KIKSFEKSNSWSWK A0301  Q13363       TRUE  
    ## # ℹ 1,882 more rows

## Select decoys

Import UniProt sequences, use the rna expression table to make an
ensembl/UniProt matching table.

Define a function that, for each allele: (1) Checks which sequences were
detected for that allele, (2) samples an excess of proteins from the
uniprot table weighted by their length, (3) then for each sample
extracts a 9mer, (4) Remove detected peptides, filter for uniqueness and
down sample to 350,000.

``` r
canon <- read_excel("nuORFs.xlsx", sheet = "Canonical")

uniprot <- read_csv("UniProt_reviewed_input.tsv") %>% rename(swissprot.id = "sequence_id")


get_peptide <- function(al, pep.len){
  MS.data %>% 
    bind_rows(canon) %>% 
    mutate(sequence = str_to_upper(sequence)) %>% 
    filter(allele == al) %>% 
    filter(nchar(sequence) == pep.len) %>% 
    pull(sequence) -> s
  
    amount.needed <- nrow(filter(filtered.peptides, allele == al & nchar(sequence) == pep.len))*1000
    
    if (amount.needed == 0) {
      tibble(sequence = NA, swissprot.id = NA, allele = NA, ligand = NA)
      return()
    }
    
    else{
      uniprot %>% 
    filter(swissprot.id %in% IDs.in.library) %>%
    mutate(len = nchar(sequence)) %>% 
    slice_sample(n = amount.needed*2, weight_by = len, replace = T) %>% 
    rowwise() %>% 
    mutate(number = sample(1:(nchar(sequence)-pep.len), 1),
         sequence = str_sub(sequence, number, number+(pep.len-1))) %>% 
    ungroup() %>% 
    filter(!(sequence %in% s)) %>% 
    transmute(sequence = sequence, swissprot.id = swissprot.id, allele = al, ligand = FALSE) %>% 
    distinct(sequence, .keep_all = T) %>% 
    slice_sample(n = amount.needed)
    }
  

}

alleles <- unique(filtered.peptides$allele)

decoys <- map2_dfr(rep(alleles, 7), rep(8:14, each = 19), get_peptide)


filtered.peptides %>% 
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

  file.list <- list.files("./Output/netMHCpan_predictions/", pattern = allele)
  file.list <- paste0("./Output/netMHCpan_predictions/", file.list)
  
  file.list %>% 
    map(read_lines) %>% 
    map(~.[grepl("   1 HLA", .)]) %>% 
    unlist() %>% 
    c(" Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL", .) %>% 
    write_lines(paste0("./Output/netMHCpan_predictions_clean/", allele,".tsv"))
  
  gc()
}
```

Apply the function for each allele

``` r
map(alleles, combine_and_clean)
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

Add to test table

``` r
test.table %>% 
  transmute(sequence = Peptide, allele = allele, rank = `%Rank_EL`) %>% 
  nest(dat = -allele) -> test.table

read_tsv("./Output/nuORF_Test_Table.tsv") %>% 
  nest(dat2 = -allele) %>% 
  inner_join(test.table) %>% 
  mutate(dat = map2(dat, dat2, ~inner_join(.x, .y, by = c("sequence")))) %>% 
  select(-dat2) %>% 
  unnest(dat) -> test.table

test.table %>% 
  count(allele, ligand)
```

``` r
write_tsv(test.table, "Output/Test_Table_nuORF_complete.tsv")
```
