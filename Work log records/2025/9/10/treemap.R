tm <- voronoiTreemap(
  df,
  levels = c("Tissue", "cluster"),
  fun = sum,
  cell_size = "prop",
  shape = "rounded_rect",#rounded_rect
  seed = 123
)

drawTreemap(
  tm,
  color_type = "categorical",
  color_level = 2,
  color_palette = cell_number1$my_colors,
  legend = F,
  label_autoscale=F,
  label_size = 3,
  label_level = 1:2,
  label_color = c('black','white'),
  border_size = c(6,1),
  border_level = 1:2
)
