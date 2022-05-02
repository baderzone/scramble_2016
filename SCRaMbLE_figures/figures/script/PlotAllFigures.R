#######################################
# Author: Giovanni Stracquadanio
# Email: stracquadanio@jhu.edu
# License: Creative Commons Attribution-NonCommercial-ShareAlike 
#
#######################################

source("PlotUtils.R")

# loading all files
source("PlotFigure1.R")
source("PlotFigure2.R")
source("PlotFigure3.R")
source("PlotFigure4.R")
source("PlotSuppFigure2.R")
source("PlotSuppFigure3.R")
source("PlotSuppFigure4.R")

PlotFigure1E("../data/scramble_contact_map.txt", "../figures/Figure1E.pdf")

PlotFigure2A("../data/Figure2C.csv", "../figures/Figure2A.pdf")
PlotFigure2B("../data//Figure2A.txt", "../figures/Figure2B.pdf")
PlotFigure2C("../data//Figure2B.txt", "../figures/Figure2C.pdf")

PlotFigure3ABis("../data/Figure3ABis.txt", "../figures/Figure3A.pdf")
PlotFigure3BBis("../data/Figure3BTer.txt", "../figures/Figure3B.pdf")

PlotFigure4A("../data/Figure4_Unique.txt", "../figures/Figure4A.pdf")
PlotFigure4B("../data/Figure4_Unique.txt", "../figures/Figure4B.pdf")
PlotFigure4C("../data/Figure4_Unique.txt", "../figures/Figure4C.pdf")
  
PlotSuppFigure2("../data//SuppFigure2Data.txt", "../figures/SuppFigure2.pdf")

PlotSuppFigure3A("../data//SuppFigure3ABis.txt", "../figures/SuppFigure3A.pdf")
PlotSuppFigure3B("../data//SuppFigure3B.txt", "../figures/SuppFigure3B.pdf")
PlotSuppFigure4("../data//SuppFigure4.txt", "../figures/SuppFigure4.pdf")


