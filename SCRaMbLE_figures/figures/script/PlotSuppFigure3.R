PlotSuppFigure3A <- function(inFilename, outFilename){
  message("Plotting SuppFigure 3A")
  
  raw = read.table(inFilename, header=T)
  data = data.frame(segment=raw$segment, 
                    cds.utr=raw$cds.utr + raw$utr.cds, 
                    cds.cds=raw$cds.cds,
                    cds.nc=raw$cds.nc + raw$nc.cds, 
                    utr.utr=raw$utr.utr, 
                    utr.nc=raw$utr.nc + raw$nc.utr, 
                    nc.nc=raw$nc.nc)
  colnames(data) = c("segment", "CDS-UTR", "CDS-CDS", "CDS-NC", "UTR-UTR", "UTR-NC", "NC-NC")
  data.m = melt(data, id="segment")
  
  
  p = ggplot(data.m, aes(x=factor(segment), y=value, fill=variable))
  p = p + geom_histogram(stat="identity", position="stack")
  p = p + scale_fill_brewer(palette="Paired",name="Type of junction")
  p = p + xlab("") + ylab("Number of junctions")
  p = p + scale_y_continuous(expand=c(0,0))
  p = p + GetBarChartTheme(GetHalfSiteColor())

  ggsave(outFilename, plot=p, width=24, height=10)
}

PlotSuppFigure3B <- function(inFilename, outFilename){
  message("Plotting SuppFigure 3B")
  
  raw = read.table(inFilename, header=T, stringsAsFactors=F)
  data = data.frame(strain=raw$strain, 
                    cds.utr=raw$cds.utr + raw$utr.cds, 
                    cds.cds=raw$cds.cds, 
                    cds.nc=raw$cds.nc + raw$nc.cds, 
                    utr.utr=raw$utr.utr, 
                    utr.nc=raw$utr.nc + raw$nc.utr, 
                    nc.nc=raw$nc.nc)
  
  data = data[order(-data$cds.utr),]
  data$strain = factor(data$strain, levels=data$strain)
  colnames(data) = c("strain", "CDS-UTR", "CDS-CDS", "CDS-NC", "UTR-UTR", "UTR-NC", "NC-NC")
  data.m = melt(data, id="strain")
  
  p = ggplot(data.m, aes(x=strain, y=value, fill=variable))
  p = p + geom_histogram(stat="identity", position="stack")
  p = p + scale_fill_brewer(palette="Paired",name="Type of junction")
  p = p + xlab("") + ylab("Number of junctions")
  p = p + scale_y_continuous(expand=c(0,0))
  p = p + GetBarChartTheme("black")
  
  ggsave(outFilename, plot=p, width=18, height=10)
}