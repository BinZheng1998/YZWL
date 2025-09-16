ggplot(res, aes(x = Expr, y = celltype_gut, fill = celltype_gut)) +
  geom_violin(scale = 'width', color = 'white', size = 0.45, alpha = 0.8) +
  facet_wrap(~ gene, scales = 'free_x', nrow = 1) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 2)) +
  theme_bw() +
  labs(x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1),
    legend.position = 'none',
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(size = 10, color = 'black'),
    strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5,color = "black",face = 'italic'),
    strip.placement = "outside",
    panel.spacing.x = unit(0.1, "lines"),
    panel.border = element_rect(size = 1, color = "black")
  )
