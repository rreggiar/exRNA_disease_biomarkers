save.manuscript.panel <- function(figure, 
				       name, 
				       width = unit(89, 'mm'), 
				       height = unit(183, 'mm'), 
				       plot.in = last_plot()) {

	if(!dir.exists(paste0(figure.dir, '/','fig.',figure))) {
		dir.create(paste0(figure.dir, '/','fig.',figure))
	}
  
	  ggsave(paste0(figure.dir,
			'/',
			'fig.',
			figure,
			'/', 
			name, 
			'.',
			round(width, 3), 
			'x', 
			round(height, 3), 
			'.pdf'),
	       width = width, 
	       height = height, 
	       units = 'mm', 
	       plot = plot.in)
}
