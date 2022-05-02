#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

source("PlotUtils.R")

PlotFigure2A <- function(input_filename, output_filename, strain_sort_decrease=T){
  message("Plotting Figure 2C")
  
  #reading data and assigning col names
  data <- read.csv(input_filename)
  colnames(data) <- c("strain",1:43)
  
  #melting data for histogram
  data.melted <- melt(data, id="strain")
  colnames(data.melted) <- c("strain", "segment", "copies") 
  
  data.melted$copies[data.melted$copies > 3] = 3
  cn1rank = rep(0, nrow(data))
  
  for (i in 1:nrow(data)){
    cn1rank[i] = sum(data[i, ] == 1)
  }
  
  data = data[order(-cn1rank), ]
  data.melted$strain = factor(data.melted$strain, levels=rev(data$strain))
  data.melted$copies = factor(data.melted$copies)
  
  # ggplot for heatmap
  p <- ggplot(data.melted, aes(x=segment, fill=copies, y=strain))
  p <- p + geom_tile(size=0.1, colour="black", position="identity") 
  p <- p + scale_fill_manual(values=c("darkgrey", "white", "orange", "red"), labels=c("0", "1", "2", ">2"))
  p <- p + labs(x = "", y = "", fill = "Copy number") 
  
  p <- p + GetTileTheme(GetSegmentColor(), rev(GetStrainColor()))
  
  #p <- p + theme(axis.text.x=element_text(angle=-90,  vjust = 0.3, colour=GetSegmentColor()), 
  #               axis.text.y = element_text(colour=rev(GetStrainColor())))
  
  # saving plot to file
  ggsave(output_filename, plot=p, width=14, height=20)
}

PlotFigure2B <- function(input_filename, output_filename){
  message("Plotting Figure 2A")
  
  # reading data and assigning 
  data <- read.table(input_filename)
  colnames(data) <- c("segment", "sample.cn0", "sample.cn1", "sample.cn2", "sample.cnmany")
  
  # melting data for plotting
  data.melted <- melt(data, id="segment")
  colnames(data.melted) <- c("segment", "cn", "samples") 
  
  # ggplot code for stacked bars
  p <- ggplot(data.melted, aes(x=factor(segment), fill=cn, y=samples))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black")
  p <- p + scale_y_continuous(expand=c(0,0))
  p <- p + scale_fill_manual(values=c("darkgrey", "white", "orange", "red"), labels=c("0", "1", "2", ">2"))
  p <- p + labs(x = "", y = "Number of strains", fill = "Copy number") 
  p <- p + GetBarChartTheme(GetSegmentColor())
  
  # saving plot to file
  ggsave(output_filename, plot=p, width=18)
}

PlotFigure2C <- function(input_filename, output_filename, strain_sort_decrease=T){
  message("Plotting Figure 2B")
  
  #reading data and assigning col names
  data <- read.table(input_filename)
  colnames(data) <- c("strain", "sample.cn0", "sample.cn1", "sample.cn2", "sample.cnmany")
  data = data[order(-data$sample.cn1), ]
  
  #melting data for histogram
  data.melted <- melt(data, id="strain")
  colnames(data.melted) <- c("strain", "cn", "segments") 
  data.melted$strain <- factor(data.melted$strain, levels=data$strain)
  
  # ggplot code for stacked bars
  p <- ggplot(data.melted, aes(x=strain, fill=cn, y=segments))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") + scale_fill_manual(values=c("darkgrey", "white", "orange", "red"), labels=c("0", "1", "2", ">2"))
  p <- p + scale_y_continuous(expand=c(0,0))
  p <- p + labs(x = "", y = "Number of segments", fill = "Copy Number") 
  
  p <- p + GetBarChartTheme(GetStrainColor())

  # saving plot to file
  ggsave(output_filename, plot=p, width=22)
  
}

