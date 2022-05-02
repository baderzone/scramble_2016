#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

source("PlotUtils.R")

PlotFigure4A <- function(input_filename, output_filename, show_legend=FALSE){
  message("Plotting Figure 4A")
  #reading data and assigning 
  table <- read.table(input_filename)
  colnames(table) <- c("strain", "event", "junction", "distance.segment", "distance.bp", "segment")
  
  #picking just useful columns
  data <- table[, c(2,4,5)]
  data$event <- as.character(data$event)
  
  #relabelling event column for facet wrap
  data$event[data$event == "SI"] = "Simple Inversion"
  data$event[data$event == "SD"] = "Simple Deletion"
  
  # ggplot code for facet histogram
  p <- ggplot(data, aes(x=distance.bp, fill=event, color=event))
  p <- p + geom_freqpoly(aes(y=..count..), binwidth=100)
  p <- p + scale_y_continuous(limits = c(0, 15))
  p <- p + scale_color_manual(name="Type of Variation", values=c("black", "red"))
  p <- p + labs(x = "Variation length (bp)", y = "Number of events", fill = "Count" )
  p = p + theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(colour="black"),
                axis.text.y = element_text(colour="black"),
                legend.position = "top",
                legend.key=element_blank())
  
  # saving to file
  ggsave(output_filename, plot=p, width=12, height=6)
  return(p)
}

PlotFigure4B <- function(input_filename, output_filename, show_legend=FALSE){
  message("Plotting Figure 4B")
  #reading data and assigning 
  table <- read.table(input_filename)
  colnames(table) <- c("strain", "event", "junction", "distance.segment", "distance.bp", "segment")
  
  #picking just useful columns
  data <- table[, c(2,4,5)]
  data$event <- as.character(data$event)
  
  #relabelling event column for facet wrap
  data$event[data$event == "SI"] = "Simple Inversion"
  data$event[data$event == "SD"] = "Simple Deletion"
  data = data[data$distance.bp <= 5000,]
  
  # ggplot code for facet histogram
  p <- ggplot(data, aes(x=distance.bp, fill=event, color=event))
  p <- p + geom_freqpoly(aes(y=..count..), binwidth=100)
  p <- p + scale_x_continuous(limits = c(0, 5000))
  p <- p + scale_y_continuous(limits = c(0, 15))
  p <- p + scale_color_manual(name="Type of Variation", values=c("black", "red"))
  p <- p + labs(x = "Variation length (bp)", y = "Number of events")
  p = p + theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(colour="black"),
                axis.text.y = element_text(colour="black"),
                legend.position = "top",
                legend.key=element_blank())
    
  # saving to file
  ggsave(output_filename, plot=p, width=12, height=6)
  return(p)
}

PlotFigure4C <- function(input_filename, output_filename, show_legend=FALSE){
  message("Plotting Figure 4C")
  #reading data and assigning 
  table <- read.table(input_filename)
  colnames(table) <- c("strain", "event", "junction", "distance.segment", "distance.bp", "segment")
  
  #picking just useful columns
  data <- table[, c(2,4,5)]
  data$event <- as.character(data$event)
  
  #relabelling event column for facet wrap
  data$event[data$event == "SI"] = "Simple Inversion"
  data$event[data$event == "SD"] = "Simple Deletion"
  data = data[data$distance.bp <= 5000 & data$distance.bp >=900,]
  
  # ggplot code for facet histogram
  p <- ggplot(data, aes(x=distance.bp, fill=event, color=event))
  p <- p + geom_freqpoly(aes(y=..count..), binwidth=100)
  p <- p + scale_x_continuous(limits = c(0, 5000))
  p <- p + scale_y_continuous(limits = c(0, 15))
  p <- p + scale_color_manual(values=c("darkgrey", "red"))
  p <- p + labs(x = "Variation length (bp)", y = "Number of events", fill = "Count" )
  p = p + theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(colour="black"),
                axis.text.y = element_text(colour="black"),
                legend.position = "top")
  
  # saving to file
  ggsave(output_filename, plot=p, width=12, height=6)
  return(p)
}

