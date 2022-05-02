require(tidyverse)
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

file2 <-
  list.files(path = dir3, pattern = "*scatter.txt", 
             full.names = TRUE) %>% 
  map_df(~read_plus(.)) %>% rename("Gene" = X1) %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(Gene) %>% distinct() %>% 
  semi_join(tpm, ., by = "Gene")
