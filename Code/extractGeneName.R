library(dplyr)
library(tidyr)
read_tsv("homerData") %>% separate(col = `Annotation/Divergence`, into = "Annotation/Divergence", sep = "\\|", remove = TRUE, extra = "drop")
