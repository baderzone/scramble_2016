#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

source("PlotUtils.R")

PlotFigure3A <- function(input_filename, output_filename){
  message("Plotting Figure 3A")
  
  #reading data and assigning column names
  table <- read.table(input_filename)
  colnames(table) <- c("loxPSym", "deletion", "deletionlost", "inversion", "complex", "absent", "parental")
  
  #reshaping column order
  #data <- data.frame(loxPSym=table$loxPSym, parental=table$parental, deletion=table$deletion, deletionlost=table$deletionlost, inversion=table$inversion, complex=table$complex, absent=table$absent)
  data <- data.frame(loxPSym=table$loxPSym, parental=table$parental, deletion=table$deletion, inversion=table$inversion, complex=table$complex)
  data$loxPSym = factor(data$loxPSym, levels=data$loxPSym)
  
  #melting data for histogram
  data.melted <- melt(data, id="loxPSym")
  colnames(data.melted) <- c("loxPSym", "event.type", "event.count") 
  
  # ggplot code for stacked bars
  p <- ggplot(data.melted, aes(x=loxPSym, fill=event.type, y=event.count))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") 
  p <- p + scale_y_continuous(expand = c(0,0))
  p <- p + scale_fill_manual(values=c("white","darkgrey", "lightblue", "red"), labels=c("Parental", "Simple Deletion", "Simple Inversion", "Complex Rearrangement"))
  p <- p + labs(x = "", y = "Number of junctions", fill = "Type of junction")
  p <- p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, colour=(GetSegmentColor())), 
                 legend.position = "top",
                 axis.line = element_line(colour="black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  
  ggsave(output_filename, plot=p, width=14)
}

PlotFigure3B <- function(input_filename, output_filename){
  message("Plotting Figure 3B")
  
  #reading data and assigning 
  table <- read.table(input_filename)
  colnames(table) <- c("strain", "deletion", "inversion", "complex", "absent", "parental")
  
  #reshaping column order
  #data <- data.frame(strain=table$strain, parental=table$parental, deletion=table$deletion, inversion=table$inversion, complex=table$complex, absent=table$absent)
  data <- data.frame(strain=table$strain, parental=table$parental, deletion=table$deletion, inversion=table$inversion, complex=table$complex)
  data$strain = factor(data$strain, levels=data$strain)
  
  #melting data for histogram
  data.melted <- melt(data, id="strain")
  colnames(data.melted) <- c("strain", "event.type", "event.count") 
  
  # ggplot source code
  p <- ggplot(data.melted, aes(x=strain, fill=event.type, y=event.count))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") 
  p <- p + scale_y_continuous(expand = c(0,0))
  p <- p + scale_fill_manual(values=c("white","darkgrey","lightblue","red"), labels=c("Parental", "Simple Deletion", "Simple Inversion", "Complex Rearrangment"))
  p <- p + labs(x = "", y = "Number of Junctions", fill = "Junction Type") 
  p <- p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, colour=(GetStrainColor())), 
                 legend.position = "top",
                 axis.line = element_line(colour="black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  
  # saving data file
  ggsave(output_filename, plot=p, width=14)
}

PlotFigure3C <- function(input_filename, output_filename, strain_sort_decrease=T){
  message("Plotting Figure 3C")
  
  #reading data and assigning col names
  data <- read.csv(input_filename)
  colnames(data) <- c("strain", get_junctions())
  #melting data for histogram
  data.melted <- melt(data, id="strain")
  print(head(data.melted))
  
  colnames(data.melted) <- c("strain", "junction", "copies") 
  data.melted$copies[data.melted$copies > 3] = 3
  data.melted$strain = factor(data.melted$strain, levels=rev(data$strain))
  data.melted$copies = factor(data.melted$copies)
  
  # ggplot for heatmap
  p <- ggplot(data.melted, aes(x=junction, fill=copies, y=strain))
  p <- p + geom_tile(size=0.1, colour="black", position="identity") 
  p <- p + scale_fill_manual(values=c("darkgrey", "white", "orange", "red"), labels=c("0", "1", "2", ">2"))
  p <- p + labs(x = "", y = "", fill = "Copy number") 
  p <- p + theme(axis.text.x=element_text(angle=-90,  vjust = 0.3, colour=GetSegmentColor()), axis.text.y = element_text(colour=rev(GetStrainColor())))
  
  # saving plot to file
  ggsave(output_filename, plot=p, width=14, height=14)
}

PlotFigure3ABis <- function(input_filename, output_filename){
  message("Plotting Figure 3ABIS")
  
  #reading data and assigning column names
  data <- read.csv(input_filename, stringsAsFactors=F)
  
  #melting data for histogram
  data.melted <- melt(data, id="segment")
  colnames(data.melted) <- c("loxPSym", "event.type", "event.count") 
  
  # ggplot code for stacked bars
  p <- ggplot(data.melted, aes(x=loxPSym, fill=event.type, y=event.count))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") 
  p <- p + scale_y_continuous(expand = c(0,0))
  p <- p + scale_fill_manual(values=c("white","darkgrey","green","lightblue", "blue", "darkblue", "orange", "red", "darkred"))
  p <- p + labs(x = "", y = "Number of junctions", fill = "Type of junction")
  p <- p + GetBarChartTheme(GetHalfSiteColor())
  ggsave(output_filename, plot=p, width=24)
}

PlotFigure3BBis <- function(input_filename, output_filename){
  message("Plotting Figure 3BBis")
  
  #reading data and assigning 
  data = read.csv(input_filename, stringsAsFactors=F)
  data = data[order(-data$parental),]
  data$strain = factor(data$strain, levels=data$strain)
  
  #melting data for histogram
  data.melted <- melt(data, id="strain")
  colnames(data.melted) <- c("strain", "event.type", "event.count") 
  
  # ggplot source code
  p <- ggplot(data.melted, aes(x=strain, y=event.count, fill=event.type))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") 
  p <- p + scale_y_continuous(expand = c(0,0))
  p <- p + scale_fill_manual(values=c("white","darkgrey","green","lightblue", "blue", "darkblue", "orange", "red", "darkred"))
  p <- p + labs(x = "", y = "Number of junctions", fill = "Type of junction") 
  #GetStrainColor
  p <- p + GetBarChartTheme(GetStrainColor())
  # saving data file
  ggsave(output_filename, plot=p, width=22)
}

