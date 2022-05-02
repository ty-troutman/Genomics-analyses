ratioPlotFunction <- function(processedData, textSize, lineWidth, pointSize) {
  
  textSize <- textSize
  lineWidth <- lineWidth / 2.835
  pointSize <- pointSize
  
  tmpVar <- list(stat = c("log2FoldChange", "padj"))
  vars2 <-
    cross2(pluck(tmpVar, "stat"), pluck(vars, "group")) %>% map(lift(str_c)) %>% as_vector()
  fc <- keep(vars2, str_detect(vars2, "log2"))
  padj <- keep(vars2, str_detect(vars2, "padj"))
  labels <- vars$labels
  
  pal_p3 <- c(viridis_pal(option = "viridis")(4))
  grob_p3 <- grobTree(
    grid.points(
      x = 0.05,
      y = 0.95,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[4]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[4], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[4]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.95,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = textSize * 0.6)
    ),
    grid.points(
      x = 0.05,
      y = 0.87,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[3]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[3], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[3]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.87,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = textSize * 0.6)
    ),
    grid.points(
      x = 0.05,
      y = 0.8,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[2]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[2], ": ",
            select(processedData, DEG) %>%
              filter(DEG == labels[2]) %>%
              tally(),
            sep = ""
      ),
      x = 0.05,
      y = 0.8,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = textSize * 0.6)
    ),
    grid.points(
      x = 0.05,
      y = 0.73,
      pch = 19,
      size = unit(pointSize*5, "point"),
      # default.units = "native", name = NULL,
      gp = gpar(col = pal_p3[1]),
      draw = TRUE,
      vp = NULL
    ),
    textGrob(
      paste(" ", labels[1], ": ",
        select(processedData, DEG) %>%
          filter(DEG == labels[1]) %>%
          tally(),
        sep = ""
      ),
      x = 0.05,
      y = 0.73,
      hjust = 0,
      gp = gpar(col = "black",
                fontsize = textSize * 0.6)
    )
  )
  
  .x <- sym(fc[2])
  .y <- sym(fc[[1]])
  
  processedData %>% 
    arrange(match(DEG, rev(labels))) %>%
    ggplot(.,
           aes(
             x = !!(.x),
             y = !!(.y),
             colour = DEG
           )) +
    ggrastr::rasterise(geom_point(size = pointSize), dpi = 600) +
    geom_abline(slope = 0, size = lineWidth) +
    geom_abline(slope = 1, size = lineWidth) +
    geom_vline(xintercept = 0, size = lineWidth) +
    theme(
      text = element_text(size = textSize, colour = "black"),
      plot.title  = element_text(size = textSize, colour = "black"),
      axis.ticks.x = element_line(size = lineWidth, colour = "black"),
      axis.ticks.y = element_line(size = lineWidth, colour = "black"),
      axis.title  = element_text(size = textSize * 0.8, colour = "black"),
      axis.text.x  = element_text(size = textSize * 0.6, colour = "black"),
      axis.text.y  = element_text(size = textSize * 0.6, colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    # coord_equal() +
    scale_color_manual(values = rev(pal_p3)) +
    annotation_custom(grob_p3) 
}