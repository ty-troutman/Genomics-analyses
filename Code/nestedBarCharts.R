################################################################################
# Author: Ty D. Troutman, PhD
# Credit: https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/
# Credit: https://www.painblogr.org/2020-06-19-purring-through-exploratory-analyses.html
# Updated: December 28, 2021
################################################################################

# A function to intersect a character list of genes of interest (GOI) with HOMER
# RNA-seq TPM data. Data is grouped via specified variable(s) and summary
# statistics are calculated. Data is then grouped by GOI and nested into a list.
# Finally, bar charts are iterated with purrr::map2() per GOI using calculated
# summary statistics. 

# Usage: Output is added to a specified object, which I call output here. It can
# then be printed to an individual file per GOI as described below:

#   output <-
#     data %>%
#     nestedBarCharts(
#       geneList = "geneList",
#       variable1 = "variable1",
#       variable2 = "variable2",
#       replicateID = "replicateID"
#     )
# print(output$data[1])
# print(output$plot[1])

nestedBarCharts <- function(data, geneList, variable1, variable2, replicateID, pointSize, lineWidth) {

  .pointSize <- pointSize
  .lineWidth <- lineWidth/2.835
  
  theme_set(theme(text = element_text(size = .pointSize, colour = "black"),
                  line = element_line(size = .lineWidth/2.835, colour = "black"), 
                  rect = element_rect(size = .lineWidth/2.835),
                  axis.ticks.x = element_line(size = .lineWidth, colour = "black"),
                  axis.ticks.y = element_line(size = .lineWidth, colour = "black"),
                  axis.text.x  = element_text(size = .pointSize*0.8, colour = "black"),
                  axis.text.y  = element_text(size = .pointSize*0.6, colour = "black"),
                  axis.title = element_text(size = .pointSize*0.8, colour = "black"),
                  plot.title = element_text(hjust = 0, face = "italic"),
                  plot.subtitle = element_text(hjust = 0),
                  legend.position = "right",
                  legend.title = element_text(size = .pointSize*0.8, colour = "black"),
                  legend.key = element_blank(),
                  legend.background = element_blank(),
                  legend.text = element_text(size = .pointSize*0.8, colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(size = .lineWidth, colour = "black"),
                  plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
                  strip.text = element_text(face = "italic", size = .pointSize, hjust = 0),
                  strip.placement = "outside",
                  strip.background = element_blank()
  ))

  summarizedData <- 
    data %>% rename("Gene" = "Annotation/Divergence") %>%
    filter(., Gene %in% goi) %>%
    pivot_longer(
      cols = -1,
      names_to = c({{variable1}}, {{variable2}}, {{replicateID}}),
      names_pattern = "(.*)_(.*)_(.*)",
      values_to = "TPM"
    ) %>%
    group_by(Gene, .data[[variable1]], .data[[variable2]]) %>%
    summarise(n = n(),
              mean = mean(TPM),
              sd = sd(TPM))
  
  nested_data <-
    summarizedData %>%
    group_by(Gene) %>%
    nest()
  
  nested_plots <-
    nested_data %>%
    mutate(plot = map2(
      data,
      Gene,
      ~ ggplot(data = .x,
               aes(
                 x = .data[[variable1]], y = mean, fill = .data[[variable2]]
               )) +
        geom_bar(
          position = position_dodge(),
          stat = "identity",
          colour = "black",
          size = .lineWidth,
          width = 0.7
        ) +
        geom_errorbar(
          aes(ymin = mean - sd, ymax = mean + sd),
          colour = "black",
          size = .lineWidth,
          width = 0.2,
          position = position_dodge(0.7)
        ) +
        expand_limits(y = 0) +
        theme(panel.spacing = unit(1, "lines")) +
        labs(
          title = Gene,
          x = NULL,
          y = c("Mean (TPM)")
        )
    ))
}
