
install.packages("tidyverse")
install.packages("psych")
install.packages("dplyr")
install.packages("edgeR")
install.packages("rgl")
install.packages("factoextra", dependencies=TRUE)
install.packages("FactoMineR", dependencies = TRUE)
install.packages('ellipse', dependencies = TRUE)
install.packages('pca3d', dependencies = TRUE)
install.packages('readxl', dependencies = TRUE)
install.packages('writexl', dependencies = TRUE)
install.packages('openxlsx', dependencies = TRUE)
install.packages('rtools', dependencies = TRUE)
install.packages('stats', dependencies = TRUE)
install.packages('effsize')

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")



library(readxl)
library(tidyverse) # modern R packages for big data analysis
library(limma) # edgeR will load this if we do not
library(edgeR)
library(dplyr)
library(rgl)
library(factoextra)
library(FactoMineR)
library(ellipse)
library(pca3d)



library(writexl)
library(xlsx)
library(effsize)
