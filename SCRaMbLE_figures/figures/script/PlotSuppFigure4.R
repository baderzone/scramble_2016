PlotSuppFigure4 <- function(input_filename, output_filename){
  message("Plotting SuppFigure 4")
  
  #reading data and assigning column names
  table <- read.csv(input_filename)
  
  #melting data for histogram
  data.melted <- melt(table, id="junction")
  colnames(data.melted) <- c("junction", "event.type", "event.count") 
  
  # ggplot code for stacked bars
  p <- ggplot(data.melted, aes(x=junction, fill=event.type, y=event.count))
  p <- p + geom_bar(stat="identity",size=0.1, colour="black") 
  p <- p + scale_y_continuous(expand = c(0,0))
  p <- p + scale_fill_manual(values=c("white","lightblue", "darkgrey", "black", "red"))
  p <- p + labs(x = "", y = "Number of junctions", fill = "Junction Type")
  p <- p + GetBarChartTheme(GetSegmentColor())
  
  ggsave(output_filename, plot=p, width=18)
}

