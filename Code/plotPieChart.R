a <- read_tsv("/gpfs/data01/scratch/glasslab/scratch_storage_usage.txt", col_names = FALSE) %>%
  rename("String" = X1, "User" = X2) %>% 
  mutate("Bytes" = as.numeric(str_extract(String, "[0-9]+"))) %>%
  mutate("Power" = str_extract(String, "[A-Z]+")) %>%
  mutate(Power = recode(Power, 'T' = 1e12, 'G' = 1e9, 'M' = 1e6, 'K' = 1e3, .missing = 1)) %>%
  mutate("Bytes" = Bytes*Power) %>%
  filter(!(User == "."))
a %>%
  ggplot(aes(
    x = "", y = Bytes/1e12, fill = reorder(User, -Bytes))) +
  geom_bar(stat="identity", width = 3 , colour = "white") +
  coord_polar("y", start = 0) +
  scale_fill_viridis(discrete = T) +
  xlab(NULL) + ylab("Terrabytes") + 
  ggtitle(paste(sum(a$Bytes)/1e12, " terrabytes used", sep = ""))
ggsave("/home/ttroutman/scratchUsage.png")

# write_tsv(a, "/home/ttroutman/scripts/exampleData/scratch_storage_usage.tsv")