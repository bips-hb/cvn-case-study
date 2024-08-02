# ***************************************************
#
#  Analysis code for the data application in the 
#  CVN paper by Dijkstra, Godt, Foraita
#
#  The data is from the KiKme project and not
#  publicly available.
#
# ***************************************************

# load libraries
library(CVN)
library(CVNSim)
library(biomaRt)          # Bioconductor package
library(dplyr)
library(gtools)
library(ggplot2)   
library(igraph)
library(tidyr)
library(viridis)         # for hamming-distance plot
source("R/functions.R")  # loads important functions



# Run files -------------------------------------------------------------------------
source("analysis/01_Preprocessing.R")
source("analysis/02_CVN-analysis.R")


