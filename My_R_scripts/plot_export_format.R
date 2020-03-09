
# PNG format
png(filename=opt$plot, units="in", width=25, height=25, pointsize=12, res=150)

# SCG format
svg(filename=opt$plot , width=15, height=15, pointsize=12)

# EPS format
cairo_ps(filename=opt$plot, width=15, height=15, pointsize=12)
