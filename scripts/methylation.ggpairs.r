library("ggplot2")
library("dplyr")
library("GGally")
library("RColorBrewer")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

args = commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],col.names=c("chr","start","end","wgbs", "count","model"))
bin2d_fn <- function(data, mapping, ...){
	p <- ggplot(data = data, mapping = mapping) + 
		geom_bin2d(bins=25) +
		scale_fill_gradientn(colors=r, trans="log10") + theme_bw()
	p
}
df_fig<-df[,c(4,5,6)]
graph<-ggpairs(df_fig,lower=list(continuous=bin2d_fn))
ggsave("methylation.ggpairs.pdf", plot = graph, width = 9, height = 9, units = "in")
