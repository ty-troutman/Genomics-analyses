maPlot <- function(tpmVar, diffVar, upperLabel, lowerLabel, mean1, mean2) {
  plot <- inner_join(tpmVar, diffVar) %>% mutate(., DEG = factor(ifelse(
    log2FoldChange < -1 & padj < 0.05, lowerLabel,
    ifelse(log2FoldChange > 1 & padj < 0.05, upperLabel, "same")))) %>%
    mutate(meanTPM = log2(
      rowMeans(select(., contains(mean1) | contains(mean2)))
      + 1)) %>% 
    arrange(., desc(DEG))
  upper1 <- grobTree(textGrob(paste(upperLabel, select(plot, DEG) %>% 
                                      filter(DEG == upperLabel) %>% tally(), 
                                    sep = ": "), x=0.05,  y=0.95, hjust=0, 
                              gp=gpar(col="black", fontsize=pointSize*0.75)))
  lower1 <- grobTree(textGrob(paste(lowerLabel, select(plot, DEG) %>% 
                                      filter(DEG == lowerLabel) %>% tally(), 
                                    sep = ": "), x=0.05,  y=0.05, hjust=0, 
                              gp=gpar(col="black", fontsize=pointSize*0.75)))
  ggplot(plot, aes(x = meanTPM, 
                   y = log2FoldChange, color = DEG, label = X1)) +
    geom_point(size = 0.25) + 
    scale_color_manual(values = c("#2171b5", "#238b45", "#d9d9d9")) +
    # scale_color_viridis_d() +
    geom_abline(slope = 0, size = lineWidth) +
    xlab("Log2(Mean TPM+1)") + 
    ylab("Log2 Fold Change") +
    annotation_custom(upper1) +
    annotation_custom(lower1) +
    scale_y_continuous(expand = c(.25, .05)) +
    # ylim(-10, 10) +
    geom_text_repel(data = filter(plot, abs(log2FoldChange) > 2 & padj < 0.05 &
                                    meanTPM > log2(32+1)),
                    aes(x = meanTPM, y = log2FoldChange, label = X1),
                    fontface = 'italic', size = 8/2.835, colour = "#E41A1C")
}
