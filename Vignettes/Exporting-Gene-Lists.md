Exporting gene lists with R
================

------------------------------------------------------------------------

``` r
library(tidyverse)
```

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ✓ readr   2.1.2     ✓ forcats 0.5.1

    Warning: package 'tidyr' was built under R version 4.1.2

    Warning: package 'readr' was built under R version 4.1.2

    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    x dplyr::filter() masks stats::filter()
    x dplyr::lag()    masks stats::lag()

``` r
tb1 <-
  readxl::read_excel(path = "../Data/Source/1-s2.0-S107476132030159X-mmc3.xlsx") %>%
  dplyr::select(1:6) %>%
  dplyr::rename(
    accession = 1,
    gene = 2,
    kch_rep1 = 3,
    kch_rep2 = 4,
    kn_rep1 = 5,
    kn_rep2 = 6
  ) %>%
  mutate(across(.cols = 3:6,
                .fns = ~ log2(. + 1))) %>%
  mutate(kch_mean = (kch_rep1 + kch_rep2) / 2,
         kn_mean = (kn_rep1 + kn_rep2) / 2) %>%
  mutate(logFC = kch_mean - kn_mean) %>%
  separate(
    col = 2,
    into = "gene",
    sep = "\\|",
    remove = TRUE,
    extra = "drop"
  )
```

``` r
nash <-
  tb1 %>% 
  filter(logFC < -1) %>%
  distinct(gene) %>%
  write_csv(file = "../Data/Processed/nash.csv")

control <-
  tb1 %>% 
  filter(logFC > 1) %>%
  distinct(gene) %>%
  write_csv(file = "../Data/Processed/control.csv")
```
