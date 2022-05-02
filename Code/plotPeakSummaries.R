lapply(X = c("readr", "ggplot2", "dplyr", "tidyr"), 
       FUN = library, character.only = TRUE)
read_tsv("peak_summary.txt") %>%
  select(-"tag directory") %>%
  mutate(across(where(is.character) & !c(sID), parse_number)) %>%
  pivot_longer(cols = -1) %>%
  ggplot(aes(x = sID, y = value, fill = sID)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    axis.text.x  = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key = element_blank(),
    legend.background = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = NA, colour = "black"),
    strip.background = element_rect(fill = NA)
  ) +
  scale_fill_viridis_d(option = "viridis") +
  facet_wrap(. ~ name, scales = "free_y", ncol = 4)
ggsave("peakSummaryPlot.pdf", height = 8.5, width = 11, units = "in", useDingbats = F)
