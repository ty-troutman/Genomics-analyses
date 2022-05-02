require(tidyverse)
require(reshape2)
a <- read.delim("alignStats.txt", header = T, sep = "\t")
a[,1]
b <- melt(a)
colnames(b)
colnames(b) <- c("ID","Analyte","Value")
plot <- 
  b %>%
  ggplot(aes(x = ID, y = Value, fill = ID)) + 
  geom_bar(stat = "identity") + 
  # geom_point(pch = 21, position = position_jitterdodge()) +
  # geom_jitter(width = 0.1) +
  labs(x = NULL, y = NULL) +
  theme(#panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
    text = element_text(size = 12), 
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    axis.text.x  = element_blank(),
    # legend.position = c(0.87,0.2),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key = element_blank(),legend.background = element_rect(fill = NA),
    # axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # panel.background = element_blank(),
    panel.background = element_rect(fill=NA, colour="black"),
    strip.background = element_rect(fill=NA)) + 
  # scale_fill_brewer(palette = "Set1") +
  # scale_fill_distiller(palette = "Reds")
  scale_fill_viridis_d(option = "D") +
  # scale_color_viridis_d() +
  # expand_limits(x = 0, y = 0) +
  facet_wrap(.~Analyte, scales = "free_y", ncol = 4)
 ggsave("alignPlots.pdf", plot = plot, height = 8.5, width = 11, units = "in", useDingbats = F)

# check whether the unwanted file exists and remove it
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")
