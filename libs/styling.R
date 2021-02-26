
# Themes ------------------------------------------------------------------


my_theme <- function(){
  
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white", color = "grey80"),
        strip.text = element_text(family = "CM Sans", size = 12),
        axis.text = element_text(family = "CM Sans", size = 10),
        axis.title = element_text(family = "CM Sans", size = 12),
        title = element_text(family = "CM Sans", size = 16),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border =  element_rect(fill = "NA", color = "grey80"),
        legend.key         = element_rect(fill = "white"))
}


# Colors Benchmark --------------------------------------------------------

my_colors <- c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F", '#af7ac5', '#E09F3E', '#335C67')
names(my_colors) <- c('salmonTE', 'SQuIRE', 'TEtools', 'TEtranscripts', 'Telescope', 'single', 'paired')

my_color_scale <- scale_color_manual(name ="Tool", values = my_colors)

my_fill_scale <- scale_fill_manual(name = 'Tool', values = my_colors)

my_shapes <- c(19, 17)
names(my_shapes) <- c('ETEs', 'with DETEs')

my_shape_scale <- scale_shape_manual(name = "Simulation", values = my_shapes)