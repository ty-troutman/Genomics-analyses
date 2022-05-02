################################################################################
# Author: Ty Troutman, PhD
# March 23, 2021

# Note (2021-12-28): This function should be updated to accept arguments instead
# of the explicit variables in lines 46, 47, and 62.
################################################################################



# A function to intersect genes of interest (GOI) with HOMER RNA-seq TPM data.
# Data is grouped via specified variable(s) and summary statistics calculated.
# Data is then grouped by GOI and nested into a list.
# Finally, barplots are iterated by summary statistics per GOI using map2().
# Output is added to a specified object:
#   .x$data[1]
# print(.x$plot[1])

nestedBarplots <- function(tpmData, geneList) {

  pointSize <- 10
  lineWidth <- 0.5/2.835
  theme_set(theme(text = element_text(size = pointSize, colour = "black"),
                  axis.ticks.x = element_line(size = lineWidth, colour = "black"),
                  axis.ticks.y = element_line(size = lineWidth, colour = "black"),
                  axis.text.x  = element_text(size = pointSize*0.75, colour = "black"),
                  axis.text.y  = element_text(size = pointSize*0.6, colour = "black"),
                  plot.title = element_text(hjust = 0),
                  plot.subtitle = element_text(hjust = 0),
                  legend.position = "right",
                  legend.title = element_text(size = pointSize*0.75, colour = "black"),
                  legend.key = element_blank(),
                  legend.background = element_blank(),
                  legend.text = element_text(size = pointSize*0.6, colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(size = lineWidth, colour = "black"),
                  plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
                  strip.text = element_text(face = "italic", size = pointSize, hjust = 0),
                  strip.placement = "outside",
                  strip.background = element_blank()
  ))

  summarizedData <- tpmData %>% filter(., Gene %in% nativeC57) %>% #recode nativeC57 as argument variable
    pivot_longer(cols = -1, names_to = c("Strain", "Condition", "Rep"), #recode quoted text as argument variables
                 names_pattern = "(.*)_(.*)_(.*)",
                 values_to = "TPM") %>%
    group_by(Gene, Strain, Condition) %>%
    summarise(n=n(), mean=mean(TPM), sd=sd(TPM))

  nested_data <-
    summarizedData %>%
    group_by(Gene) %>%
    nest()

  nested_plots <-
    nested_data %>%
    mutate(plot = map2(data, Gene,
                       ~ ggplot(data = .x,
                                aes(x = Condition, y = mean, fill = Strain)) + #recode
                         geom_bar(position = position_dodge(), stat ="identity", colour="black",
                                  size = lineWidth, width = 0.7) +
                         geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
                                       colour = "black", size = lineWidth, width = 0.2,
                                       position = position_dodge(0.7)) +
                         expand_limits(x = 0, y = 0) +
                         theme(panel.spacing = unit(1, "lines")) +
                         # facet_wrap(~Gene, nrow = 1, scales = "free_y") +
                         labs(title = Gene, x = NULL, y = c("Mean (TPM)")) +
                         scale_fill_manual(values = c("#377EB8", "#4DAF4A")),
                         # scale_fill_brewer(palette = "Set1")
                         ))
}