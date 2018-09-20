library(stringr)
library(readxl)
library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(gridExtra)



#set to require in function
library(pca3d)
library(openxlsx)
library(rgl)


#not needed anymore?

library(factoextra)
library(FactoMineR)
library(ellipse)





#install section for packages, some not ended up being used at this time

install.packages("stringr")
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
install.packages("openxlsx")
install.packages("sjstats")
install.packages('gridExtra')

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")



install.packages ("XLConnect")
install.packages("openxlsx", dependencies = TRUE)




library(sjstats)
library(writexl)

install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")


library(effsize)
