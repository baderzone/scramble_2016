#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

source("PlotUtils.R")

PlotSuppFigure2 <- function(inFilename, outFilename){
  message("Plotting Supp Figure 2")
  fwer = 0.05/43
  
  data <- read.table(inFilename, header=T)
  data$significant = data$p_two <= fwer
  
  p <- ggplot(data, aes(x=segment, y=-log10(p_two), color=significant))
  p <- p + geom_point(aes(size=5))
  p <- p + geom_hline(yintercept=-log10(fwer), linetype="dashed")
  p <- p + scale_x_continuous(breaks=c(1:43))
  p <- p + scale_color_manual(values=c("darkgrey", "red"))
  p <- p + annotate("text", label = "IST3", x = 6, y = -log10(data$p_two[6])+.1, size = 6, colour = "red")
  p <- p + xlab("") + ylab(expression(-log[10](p-value)))
  p <- p + GetScatterTheme(GetSegmentColor())
  
  ggsave(outFilename, plot=p, width=14, height=8)
}



