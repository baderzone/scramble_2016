#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

#loading ggplot2 + reshape libraries
library(ggplot2)
library(grid)
library(reshape)
library(RColorBrewer)


GetSegmentColor <- function(){
  segment_color <- rep("black", 43)
  # centromer
  segment_color[2] <- "steelblue1"
  # outlining markers in purple
  segment_color[c(14,32)] <- "magenta3"
  # essential genes
  segment_color[c(7,9,10,12,20)] <- "red"
  
  return(segment_color)
}

GetHalfSiteColor <- function(){
  segment_color <- rep("black", 86)
  # centromer
  segment_color[3:4] <- "steelblue1"
  # outlining markers in purple
  segment_color[c(27,28,63:64)] <- "magenta4"
  # essential genes
  segment_color[c(13,14,17,18,19,20,23,24,39,40)] <- "red"
  
  return(segment_color)
}

GetStrainColor <- function(){
  #strain_color <- c(rep("darkgreen",41), rep("orange",8), rep("red", 7), rep("purple",8), rep("black", 2))
  strain_color = rep("black", 66)
  return(strain_color)
}

GetJunctionLabel <- function(){
  junctions <- c("1R|2L","2R|3L","3R|4L","4R|5L","5R|6L","6R|7L","7R|8L","8R|9L","9R|10L","10R|11L","11R|12L","12R|13L","13R|14L","14R|15L","15R|16L","16R|17L","17R|18L","18R|19L","19R|20L","20R|21L","21R|22L","22R|23L","23R|24L","24R|25L","25R|26L","26R|27L","27R|28L","28R|29L","29R|30L","30R|31L","31R|32L","32R|33L","33R|34L","34R|35L","35R|36L","36R|37L","37R|38L","38R|39L","39R|40L","40R|41L","41R|42L","42R|43L","43R|1L")
  return(junctions)  
}

GetBarChartTheme <- function(xcolor){
  return (theme(legend.position = "top",
        axis.line = element_line(colour="black"),                
        axis.text.x=element_text(angle=-90,  vjust = 0.3, colour=xcolor, size=20),
        axis.text.y=element_text(size=20, colour="black"),
        axis.title.x=element_text(size=20, colour="black"),
        axis.title.y=element_text(size=25, colour="black"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()))
}

GetTileTheme <- function(xcolor,ycolor){
  return (theme(legend.position = "top",
                axis.line = element_line(colour="black"),                
                axis.text.x=element_text(angle=-90,  vjust = 0.3, colour=xcolor, size=20),
                axis.text.y=element_text(size=20, colour=ycolor),
                axis.title.x=element_text(size=20, colour="black"),
                axis.title.y=element_text(size=25, colour="black"),
                legend.text=element_text(size=20),
                legend.title=element_text(size=20),        
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()))
}

GetScatterTheme <- function(xcolor){
  return (theme(legend.position = "none",
                axis.line = element_line(colour="black"),                
                axis.text.x=element_text(angle=-90,  vjust = 0.3, colour=xcolor, size=20),
                axis.text.y=element_text(size=20, colour="black"),
                axis.title.x=element_text(size=20, colour="black"),
                axis.title.y=element_text(size=25, colour="black"),
                legend.text=element_text(size=20),
                legend.title=element_text(size=20),        
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()))
}



