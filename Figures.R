

#### Diversity ####
design_diversity="
AB
AB
AB
CC"

Phylum.barplot+BC_16S+legionella.barplot+
  plot_layout(design=design_diversity,guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")

ggsave("/home/jesus/Pierre/madalina/Figures/diversity.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)
