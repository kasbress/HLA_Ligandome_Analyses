XGB results top 0.1% Broad data
================
Kaspar Bresser
28/03/2024

- [Import and tidy data](#import-and-tidy-data)
- [Positive predictive value](#positive-predictive-value)
- [Area under the curve](#area-under-the-curve)
- [Sensitivity ect](#sensitivity-ect)

Used the analysis below to compare the performance of the XGB models.

``` r
library(bayestestR)
library(tidyverse)
library(pROC)
library(rstatix)
library(ggpubr)
```

## Import and tidy data

Import the results table.

``` r
XGB.results <- read_tsv("./Output/all_predictions_Lung.tsv")
```

Convert the ligand column to a binary so it can be used for cumsum
calculation, pivot to longer data by putting the models in a grouped
format.

``` r
XGB.results %>% 
  mutate(rank = -rank) %>% 
  mutate(detected = case_when(ligand == TRUE ~ 1,
                              TRUE ~ 0)) %>% 
  pivot_longer(cols = !c("detected","ligand", "tumor", "sequence", "swissprot.id"), 
               names_to =  "model", 
               values_to =  "score") -> XGB.results

XGB.results
```

    ## # A tibble: 20,019,550 × 7
    ##    ligand sequence    swissprot.id tumor     detected model          score
    ##    <lgl>  <chr>       <chr>        <chr>        <dbl> <chr>          <dbl>
    ##  1 TRUE   VAYPHDGKIFF Q9UKN5       C3N-02287        1 aff_lib        0.808
    ##  2 TRUE   VAYPHDGKIFF Q9UKN5       C3N-02287        1 aff_lib_small  0.539
    ##  3 TRUE   VAYPHDGKIFF Q9UKN5       C3N-02287        1 aff_only       0.906
    ##  4 TRUE   VAYPHDGKIFF Q9UKN5       C3N-02287        1 rank          -0.277
    ##  5 TRUE   VAYPHDGKIFF Q9UKN5       C3N-02287        1 rankExp        0.821
    ##  6 TRUE   SVLRSFRVAK  Q6UXN9       C3N-02287        1 aff_lib        0.918
    ##  7 TRUE   SVLRSFRVAK  Q6UXN9       C3N-02287        1 aff_lib_small  0.920
    ##  8 TRUE   SVLRSFRVAK  Q6UXN9       C3N-02287        1 aff_only       0.906
    ##  9 TRUE   SVLRSFRVAK  Q6UXN9       C3N-02287        1 rank          -0.259
    ## 10 TRUE   SVLRSFRVAK  Q6UXN9       C3N-02287        1 rankExp        0.255
    ## # ℹ 20,019,540 more rows

For these analysis we’ll focus on the top 0.1% scoring peptides for each
model. Let’s subset on those. Note that I arrange on both both model
scores and a random number for tie breaking.

Next we’ll calculate the cumulative sums by grouping by allele and
model.

Note that `cumsum()` takes an ordered dataframe.

``` r
XGB.results %>% 
  group_by(tumor, model) %>% 
  mutate(random = sample(1:n())) %>% 
  arrange(desc(score), random, .by_group = T) %>% 
    mutate(peptides = 1:n()/n(), 
         detected_ligands = cumsum(detected)/sum(detected)) %>% 
  slice_head(prop = .001) -> XGB.results
 
XGB.results
```

    ## # A tibble: 20,000 × 10
    ## # Groups:   tumor, model [40]
    ##    ligand sequence   swissprot.id tumor     detected model score random peptides
    ##    <lgl>  <chr>      <chr>        <chr>        <dbl> <chr> <dbl>  <int>    <dbl>
    ##  1 TRUE   ALKNPPINTK O15511       C3N-02287        1 aff_… 1.00  343167  2.00e-6
    ##  2 TRUE   TENDIANFF  P31942       C3N-02287        1 aff_… 1.00  344309  4.00e-6
    ##  3 FALSE  KAQYEDIAQK P04264       C3N-02287        0 aff_… 1.00  463885  5.99e-6
    ##  4 TRUE   EEQKKLEEL  Q9Y2S6       C3N-02287        1 aff_… 1.00   59372  7.99e-6
    ##  5 FALSE  SHVDPIHIF  Q6P6C2       C3N-02287        0 aff_… 0.999 225595  9.99e-6
    ##  6 FALSE  FTKEAYKKY  P13693       C3N-02287        0 aff_… 0.999 383926  1.20e-5
    ##  7 FALSE  GLFNRIIRK  Q9H0A0       C3N-02287        0 aff_… 0.999 159105  1.40e-5
    ##  8 FALSE  KMMDKVEFV  P49591       C3N-02287        0 aff_… 0.999  19908  1.60e-5
    ##  9 FALSE  KMNPPKFSK  P35580       C3N-02287        0 aff_… 0.999 158628  1.80e-5
    ## 10 TRUE   SESNLRKAF  P55072       C3N-02287        1 aff_… 0.999 402101  2.00e-5
    ## # ℹ 19,990 more rows
    ## # ℹ 1 more variable: detected_ligands <dbl>

Lastly, set the ordering in which we’d like the models to appear in
plots, by converting the models variable to a factor

``` r
unique(XGB.results$model)
```

    ## [1] "aff_lib"       "aff_lib_small" "aff_only"      "rank"         
    ## [5] "rankExp"

``` r
model.order <- c("aff_only", "aff_lib_small","aff_lib")

tumor.order <- as.character(unique(XGB.results$tumor))

colors <- c("#6BAED6", "#756CB1", "#5F3D8C")

XGB.results %>% 
  filter(model %in% model.order) %>% 
  mutate(model = as_factor(model),
         model = fct_relevel(model, model.order),
         tumor = as_factor(tumor),
         tumor = fct_relevel(tumor, tumor.order)) -> XGB.results

XGB.results
```

    ## # A tibble: 12,000 × 10
    ## # Groups:   tumor, model [24]
    ##    ligand sequence   swissprot.id tumor     detected model score random peptides
    ##    <lgl>  <chr>      <chr>        <fct>        <dbl> <fct> <dbl>  <int>    <dbl>
    ##  1 TRUE   ALKNPPINTK O15511       C3N-02287        1 aff_… 1.00  343167  2.00e-6
    ##  2 TRUE   TENDIANFF  P31942       C3N-02287        1 aff_… 1.00  344309  4.00e-6
    ##  3 FALSE  KAQYEDIAQK P04264       C3N-02287        0 aff_… 1.00  463885  5.99e-6
    ##  4 TRUE   EEQKKLEEL  Q9Y2S6       C3N-02287        1 aff_… 1.00   59372  7.99e-6
    ##  5 FALSE  SHVDPIHIF  Q6P6C2       C3N-02287        0 aff_… 0.999 225595  9.99e-6
    ##  6 FALSE  FTKEAYKKY  P13693       C3N-02287        0 aff_… 0.999 383926  1.20e-5
    ##  7 FALSE  GLFNRIIRK  Q9H0A0       C3N-02287        0 aff_… 0.999 159105  1.40e-5
    ##  8 FALSE  KMMDKVEFV  P49591       C3N-02287        0 aff_… 0.999  19908  1.60e-5
    ##  9 FALSE  KMNPPKFSK  P35580       C3N-02287        0 aff_… 0.999 158628  1.80e-5
    ## 10 TRUE   SESNLRKAF  P55072       C3N-02287        1 aff_… 0.999 402101  2.00e-5
    ## # ℹ 11,990 more rows
    ## # ℹ 1 more variable: detected_ligands <dbl>

## Positive predictive value

Lets first plot the positive predictive value (PPV), i.e. the number of
true positives within a certain threshold. Calculate the PPV as the
number of true positives divided by the total number of peptides within
the selected pool.

``` r
XGB.results %>%  
  group_by(tumor, model) %>% 
  mutate(model = fct_drop(model)) %>% 
  summarise(PPV = sum(detected)/n()) %>% 
  ungroup() %>% 
  rstatix::t_test(PPV ~ model, paired = T) %>% 
  add_xy_position(x = "model",  fun = "max",step.increase = .04,) %>% 
  mutate(p.adj = round(p.adj, 3))-> stats
```

``` r
XGB.results %>%  
  group_by(tumor, model) %>% 
  summarise(PPV = sum(detected)/n()) %>% 
   ggplot( aes(x = model, y = PPV, group = tumor)) +
    geom_boxplot(aes(group = model), width = 0.45)+
    geom_line()+
    geom_point(aes(fill = model), shape = 21, size = 2.5 )+
    scale_fill_manual(values = colors, labels = c("Aff", "Aff+Lib", "Aff+Lib 2.0"))+
    stat_pvalue_manual(stats, label = "p.adj", 
                     tip.length = 0.00, hide.ns = F, label.size = 2 )+
    ggtitle("PPV in top 0.1%")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          panel.grid.major.y = element_line(), axis.text.x=element_blank(), legend.position = "bottom",)
```

<img src="3-XGB_model_performance_files/figure-gfm/PVV_plot_increment2-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/XGB_PPV_strip.pdf", width = 1.7, height = 3.2, scale = 1.2)
```

## Area under the curve

Or calculate the AUC of these curves as a summary metric and plot that.
Used the `area_under_curve()` function from the
[bayestestR](https://easystats.github.io/bayestestR/index.html) package.

``` r
XGB.results %>% 
  group_by(tumor, model) %>% 
  summarise(auc = area_under_curve(peptides, detected_ligands)*1000) %>% 
  mutate(model = fct_drop(model)) %>%
  ungroup() %>% 
  rstatix::t_test( auc ~ model, paired = T) %>% 
  add_xy_position(x = "model", step.increase = 0.04) %>% 
  mutate(p.adj = round(p.adj, 4)) -> stats

stats
```

    ## # A tibble: 3 × 14
    ##   .y.   group1     group2    n1    n2 statistic    df       p p.adj p.adj.signif
    ##   <chr> <chr>      <chr>  <int> <int>     <dbl> <dbl>   <dbl> <dbl> <chr>       
    ## 1 auc   aff_only   aff_l…     8     8    -0.533     7 6.1 e-1 0.61  ns          
    ## 2 auc   aff_only   aff_l…     8     8    -4.98      7 2   e-3 0.003 **          
    ## 3 auc   aff_lib_s… aff_l…     8     8    -6.11      7 4.84e-4 0.001 **          
    ## # ℹ 4 more variables: y.position <dbl>, groups <named list>, xmin <dbl>,
    ## #   xmax <dbl>

``` r
XGB.results %>% 
  group_by(tumor, model) %>% 
  summarise(auc = area_under_curve(peptides, detected_ligands)*1000) %>% 
    ggplot( aes(x = model, y = auc, group = tumor)) +
    geom_boxplot(aes(group = model), width = 0.45)+
    geom_line()+
    geom_point(aes(fill = model), shape = 21, size = 2.5 )+
    scale_fill_manual(values = colors, labels = c("Aff", "Aff+Lib", "Aff+Lib 2.0"))+
    stat_pvalue_manual(stats, label = "p.adj", 
                     tip.length = 0.00, hide.ns = F, label.size = 2 )+
    ggtitle("AUC in top 0.1%")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          panel.grid.major.y = element_line(), axis.text.x=element_blank(), legend.position = "bottom",)
```

<img src="3-XGB_model_performance_files/figure-gfm/PVV_plot_increment3-1.png" style="display: block; margin: auto;" />

``` r
ggsave("Figs/XGB_AUC_strip.pdf", width = 1.7, height = 3.2, scale = 1.2)
```

## Sensitivity ect

Read in the data to get metrics of model performance

``` r
library(pROC)

XGB.results <- read_tsv("./Output/all_predictions_Lung.tsv")

XGB.results <- mutate(XGB.results, rank = -rank)
```

Define a function to retrieve the performance metrics, and a function to
get the cutoffs at a certain threshold of sensitivity

``` r
get_sens_spec <- function(lig, preds, mods){
 
  map(preds, ~roc(lig, .)) %>% 
     map(~coords(., x = "all", input = "sensitivity", ret = "all") ) -> tmp

    
  tibble(model = mods, dat = tmp) %>% 
    unnest()
    
}

get_cutoffs <- function(lig, preds, mods, thresh){
  map(preds, ~roc(lig, .)) %>%
    map(~coords(., x = thresh, input = "sensitivity", ret = "all") ) %>% 
  map(~select(., !c(threshold))) %>% 
  bind_rows() %>% 
    mutate(model = mods)
}
```

Get the information

``` r
model.results.info <- get_sens_spec(lig = XGB.results$ligand,
                                    preds = map2(list(XGB.results), model.order, ~pull(.x, .y)),
                                    mods =  model.order)
```

    ## Warning: `cols` is now required when using `unnest()`.
    ## ℹ Please use `cols = c(dat)`.

``` r
coordinates <- get_cutoffs(lig = XGB.results$ligand,
                           preds = map2(list(XGB.results), model.order, ~pull(.x, .y)),
                           mods =  model.order,
                           thresh = 0.1)

model.results.info
```

    ## # A tibble: 792,034 × 25
    ##    model    threshold specificity sensitivity accuracy     tn    tp    fn     fp
    ##    <chr>        <dbl>       <dbl>       <dbl>    <dbl>  <dbl> <dbl> <dbl>  <dbl>
    ##  1 aff_on… -Inf             0           1     0.000999 0       4000     0 4.00e6
    ##  2 aff_on…    0.00565       0.138       0.994 0.139    5.52e5  3977    23 3.45e6
    ##  3 aff_on…    0.00693       0.159       0.992 0.160    6.38e5  3970    30 3.36e6
    ##  4 aff_on…    0.00819       0.417       0.978 0.418    1.67e6  3914    86 2.33e6
    ##  5 aff_on…    0.0112        0.439       0.977 0.440    1.76e6  3907    93 2.24e6
    ##  6 aff_on…    0.0158        0.498       0.973 0.499    1.99e6  3891   109 2.01e6
    ##  7 aff_on…    0.0278        0.531       0.971 0.532    2.13e6  3885   115 1.87e6
    ##  8 aff_on…    0.0447        0.584       0.968 0.584    2.33e6  3871   129 1.66e6
    ##  9 aff_on…    0.0647        0.608       0.966 0.608    2.43e6  3866   134 1.57e6
    ## 10 aff_on…    0.0950        0.685       0.958 0.686    2.74e6  3831   169 1.26e6
    ## # ℹ 792,024 more rows
    ## # ℹ 16 more variables: npv <dbl>, ppv <dbl>, fdr <dbl>, fpr <dbl>, tpr <dbl>,
    ## #   tnr <dbl>, fnr <dbl>, `1-specificity` <dbl>, `1-sensitivity` <dbl>,
    ## #   `1-accuracy` <dbl>, `1-npv` <dbl>, `1-ppv` <dbl>, precision <dbl>,
    ## #   recall <dbl>, youden <dbl>, closest.topleft <dbl>

And plot at 0.1 sensitivity

``` r
ggplot(model.results.info, aes(x = `1-specificity`, y = sensitivity, color = model))+
  geom_line( linewidth = 1.5)+
  coord_cartesian(ylim = c(0.0, 0.1), xlim = c(0.0, 0.0005))+
  geom_segment(aes(x = 0.0, xend = max(coordinates$`1-specificity`), y = .1, yend = .1), color = "black", linetype = "dashed")+
  geom_segment(data = coordinates, aes(x = `1-specificity`, xend = `1-specificity`, y = 0.0, yend = .1, color = model), 
               linewidth = 0.8, linetype = "dashed")+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major = element_line())
```

    ## Warning in geom_segment(aes(x = 0, xend = max(coordinates$`1-specificity`), : All aesthetics have length 1, but the data has 792034 rows.
    ## ℹ Did you mean to use `annotate()`?

<img src="3-XGB_model_performance_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

``` r
ggsave(("Figs/XGB_SensSpec_10.pdf"), width = 2.8, height = 2, scale = 2)
```

    ## Warning in geom_segment(aes(x = 0, xend = max(coordinates$`1-specificity`), : All aesthetics have length 1, but the data has 792034 rows.
    ## ℹ Did you mean to use `annotate()`?

And plot at 0.5 sensitivity

``` r
coordinates <- get_cutoffs(lig = XGB.results$ligand,
                           preds = map2(list(XGB.results), model.order, ~pull(.x, .y)),
                           mods =  model.order,
                           thresh = 0.5)

ggplot(model.results.info, aes(x = `1-specificity`, y = sensitivity, color = model))+
  geom_line( linewidth = 1.5)+
  coord_cartesian(ylim = c(0.0, 0.5), xlim = c(0.0, 0.005))+
  geom_segment(aes(x = 0.0, xend = max(coordinates$`1-specificity`), y = .5, yend = .5), color = "black", linetype = "dashed")+
  geom_segment(data = coordinates, aes(x = `1-specificity`, xend = `1-specificity`, y = 0.0, yend = .5, color = model), 
               linewidth = 0.8, linetype = "dashed")+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major = element_line())
```

    ## Warning in geom_segment(aes(x = 0, xend = max(coordinates$`1-specificity`), : All aesthetics have length 1, but the data has 792034 rows.
    ## ℹ Did you mean to use `annotate()`?

<img src="3-XGB_model_performance_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
ggsave(("Figs/XGB_SensSpec_50.pdf"), width = 2.8, height = 2, scale = 2)
```

    ## Warning in geom_segment(aes(x = 0, xend = max(coordinates$`1-specificity`), : All aesthetics have length 1, but the data has 792034 rows.
    ## ℹ Did you mean to use `annotate()`?

``` r
ggplot(model.results.info, aes(x = `1-specificity`, y = sensitivity, color = model))+
  geom_line( linewidth = 1.5)+
  coord_cartesian(ylim = c(0.0, 1), xlim = c(0.0, 0.05))+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major = element_line())
```

<img src="3-XGB_model_performance_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />
