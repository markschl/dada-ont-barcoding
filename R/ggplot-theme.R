
default_theme <- function(base_size=9, ...) {
  library(ggplot2)
  theme_minimal(base_size, ...) + theme(
    text = element_text(),
    panel.background = element_rect(colour = NA),
    plot.title = element_text(size = rel(1.1)),
    panel.border = element_rect(
      color = 'black',
      fill = NA,
      linetype = 1,
      size = 0.4
    ),
    axis.title = element_text(size = rel(1)),
    axis.text = element_text(colour = 'black', size = rel(1)),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(3, 'pt'),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(1)),
    legend.key = element_rect(colour = NA),
    legend.key.height = unit(0.8, 'lines'),
  )
}
