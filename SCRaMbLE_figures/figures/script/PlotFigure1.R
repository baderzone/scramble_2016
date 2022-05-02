PlotFigure1E <- function(input_filename, output_filename, filter=c("parental"), gradient=T, threshold=3){
  data = read.csv(input_filename, stringsAsFactors=T)
  
  for (i in 1:length(filter)){
    data$count[data$event == filter[i]] = NA   
  }
  
  data$event = factor(data$event, levels=c("deletion", "inversion", "complex"), labels=c("Deletion","Inversion", "Complex"))
  data$count[data$count==0] = NA
  data$count[data$count > threshold] = threshold
  print(max(data$count, na.rm=T))
  print(min(data$count, na.rm=T))
  
  p = ggplot(data, aes(x=src,y=dst))
  if (gradient){
    p = p + geom_tile(data=subset(data, !is.na(count)), aes(fill=event, alpha=count))      
  }else{
    p = p + geom_tile(data=subset(data, !is.na(count)), aes(fill=event))  
  }
  
  p = p + geom_tile(data=subset(data, is.na(count)), aes(fill=NA), name="")  
  p = p + scale_fill_manual(values=c("darkgrey", "lightblue", "red"), name="Type of recombination")
  p = p + scale_alpha_continuous(guide=F, limits=c(1.0, 3.0))
  p = p + xlab("") + ylab("")
  p = p + GetTileTheme(GetHalfSiteColor(), GetHalfSiteColor())

  ggsave(output_filename, plot=p, width=24, height=24)
}