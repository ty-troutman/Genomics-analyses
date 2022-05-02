Integrating RNA expression data with differential statistics from DeSeq2
================

Purpose:

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

Read in normalized TPM expression data

``` r
read_tsv("../Data/Source/rawT.txt") %>% colnames()
```

    New names:
    * `` -> ...1

    Rows: 23608 Columns: 27
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): ...1
    dbl (26): BALBcJ_F0_Control_BALB01C, BALBcJ_F0_Control_BALB01D, C57BL6J_F0_C...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

     [1] "...1"                       "BALBcJ_F0_Control_BALB01C" 
     [3] "BALBcJ_F0_Control_BALB01D"  "C57BL6J_F0_Control_C5701C" 
     [5] "C57BL6J_F0_Control_C5701D"  "BALBcJ_F0_LPS_20210513_1a" 
     [7] "BALBcJ_F0_LPS_20210513_1b"  "C57BL6J_F0_LPS_20210513_2b"
     [9] "C57BL6J_F0_LPS_20210513_2c" "BALBcJ_F1_Control_cb61a"   
    [11] "BALBcJ_F1_Control_cb61b"    "BALBcJ_F1_Control_cb61c"   
    [13] "BALBcJ_F1_Control_cb61d"    "C57BL6J_F1_Control_cb61a"  
    [15] "C57BL6J_F1_Control_cb61b"   "C57BL6J_F1_Control_cb61c"  
    [17] "C57BL6J_F1_Control_cb61d"   "BALBcJ_F1_LPS_cb61a"       
    [19] "BALBcJ_F1_LPS_cb61b"        "BALBcJ_F1_LPS_cb62a"       
    [21] "BALBcJ_F1_LPS_cb62b"        "BALBcJ_F1_LPS_cb62c"       
    [23] "C57BL6J_F1_LPS_cb61a"       "C57BL6J_F1_LPS_cb61b"      
    [25] "C57BL6J_F1_LPS_cb62a"       "C57BL6J_F1_LPS_cb62b"      
    [27] "C57BL6J_F1_LPS_cb62c"      

``` r
tpm <- read_tsv("../Data/Source/rawT.txt") %>% select(1:5) %>% rename(geneID = 1)
```

    New names:
    * `` -> ...1
    Rows: 23608 Columns: 27── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): ...1
    dbl (26): BALBcJ_F0_Control_BALB01C, BALBcJ_F0_Control_BALB01D, C57BL6J_F0_C...
    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Read in DeSeq2 differential expression data for these selected
comparisons

``` r
diff <- read_tsv("../Data/Source/BALBcJ_F0_Control.vs.C57BL6J_F0_Control.scatter.txt") %>% rename(geneID = 1)
```

    New names:
    * `` -> ...1

    Rows: 7664 Columns: 8
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): ...1, contrast
    dbl (6): baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Now merge data together. This will invovle using the join functions that
are part of dplyr. Dplyr has mutate joins (inner_join, full_join,
right_join, left_join) and filter joins (semi_join and anti_join),
including full_join, semi_join, left_join, right_join, inner_join, and I
probably forgot a few.

For this case we want consider genes that meet minimal expression
criteria and were included in the differential expression analysis. So
we will use inner_join

``` r
merge <-  inner_join(x = tpm, y = diff, by = "geneID")
mergeMean <- merge %>%
  mutate(meanBALB = (BALBcJ_F0_Control_BALB01C + BALBcJ_F0_Control_BALB01D) /
           2) %>%
  mutate(meanC57 = (C57BL6J_F0_Control_C5701C + C57BL6J_F0_Control_C5701D) /
           2)
```
