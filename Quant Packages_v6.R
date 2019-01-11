install.packages('packrat')
packrat::init()
packrat::snapshot()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("vsn", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSnbase", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")


install.packages("stringr")
install.packages("tidyverse", dependencies = TRUE)
install.packages("psych")
install.packages("dplyr")
install.packages("rgl", dependencies = TRUE)
install.packages('pca3d', dependencies = TRUE)

install.packages('readxl', dependencies = TRUE)
install.packages('writexl', dependencies = TRUE)
install.packages('openxlsx', dependencies = TRUE)

install.packages('tidyr')
install.packages('gridExtra')
install.packages('MASS')
install.packages('robustbase')
install.packages('gplots') 
install.packages('ggpubr')

install.packages('tibble')

#---------------------------------------------------------------------------------------------







install.packages("factoextra", dependencies=TRUE)
install.packages("FactoMineR", dependencies = TRUE)
install.packages('ellipse', dependencies = TRUE)




install.packages('rtools', dependencies = TRUE)
install.packages('stats', dependencies = TRUE)
install.packages('effsize')
install.packages("openxlsx")
install.packages("sjstats")
install.packages('gridExtra')
install.packages('robustbase')

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")

install.packages('jsonlite')
#install.packages ("XLConnect")
install.packages("openxlsx", dependencies = TRUE)

install.packages("devtools")

install_github("easyGgplot2", "kassambara")


library(effsize)

install.packages("Rcmdr")
install.packages('PerformanceAnalytics', 'ape', 'raster' )




install.packages('statmod')


#install.packages('bnstruct')
#library(bnstruct)
#remove.packages('bnstruct')



install.packages('hexbin')





