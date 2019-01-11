rm(list = ls())
library(tcltk)
library(tidyr)
library(stringr)
library(readxl)
library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(gridExtra)
library(MASS)
library(pcaMethods)
library(preprocessCore)
library(vsn)
library(robustbase)
library(gplots) 
library(limma)
library(ggpubr)
library(tibble)
library(rgl)
library(pca3d)

#-------------------------------------------
# User input
#-------------------------------------------

# Design and Sample Info
#study_design <-"5062 121218 Sample Info.xlsx"
study_design <- file.choose()
input_data <- file.choose()
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


#color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black")
color_choices <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
options(digits=3)

source("Quant Functions Misc v6.R")
source("Quant Functions Norm v6.R")
source("Quant Functions Impute v6.R")
source("Quant Functions Stats v6.R")
source("Quant Functions Plots v6.R")

design_info <- read_excel(study_design, sheet="Design")
sample_info <- read_excel(study_design, sheet="SampleList")


#input_data <- as.character(design_info[1,2])
file_prefix <- as.character(design_info[2,2])
psm_input <- as.logical(design_info[3,2])
psm_to_peptide <- as.logical(design_info[4,2])
protein_peptide_input <- as.logical(design_info[5,2])
peptide_to_protein <- as.logical(design_info[6,2])
normalize_to_protein <- as.logical(design_info[7,2])
pair_comp <- as.logical(design_info[9,2])
adh_spike <- as.logical(design_info[10,2])
missings <- as.character(design_info[11,2])
intensity_cutoff <- as.numeric(design_info[12,2])
area_floor <- as.numeric(design_info[13,2])
pvalue_cutoff <- as.numeric(design_info[14,2])
fc_cutoff <- as.numeric(design_info[15,2])

adh_list <- as.character(na.omit(t(design_info[16,2:11])))
bait_list <- as.character(na.omit(t(design_info[17,2:11])))
avidin_list <- as.character(na.omit(t(design_info[18,2:11])))
casein_list <- as.character(na.omit(t(design_info[19,2:11])))
carbox_list <- as.character(na.omit(t(design_info[20,2:11])))
bira_list <- as.character(na.omit(t(design_info[21,2:11])))
ko1_list <- as.character(na.omit(t(design_info[22,2:11])))
ko2_list <- as.character(na.omit(t(design_info[23,2:11])))
ko3_list <- as.character(na.omit(t(design_info[24,2:11])))
ko4_list <- as.character(na.omit(t(design_info[25,2:11])))
protein_norm_list <- as.character(na.omit(t(design_info[8,2:11])))

misaligned_filter <- as.logical(design_info[26,2])
impute_method <- as.character(design_info[27,2])

#read data file
forward_data <- read_excel(input_data, 1)
# number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file
sample_number <- nrow(sample_info)
comp_number <- length(grep(x = colnames(sample_info), pattern = "^Comp"))
#count number of samples in each group - add column
sample_info$Count <- sapply(sample_info$Group, function(string) sum(string==sample_info$Group))

# create unqiue group dataframe with sample number
sample_groups <- sample_info[6:ncol(sample_info)]
sample_groups <- sample_groups %>% distinct(Group, .keep_all = TRUE)
group_number <- nrow(sample_groups)

#assign start and end columns for each group in pd output
sample_groups$start <- 1
sample_groups$end <- sample_groups$Count[1]
for(i in 2:(group_number)) {
  sample_groups$start[i] <- sample_groups$start[i-1] + sample_groups$Count[i-1]
  sample_groups$end[i] <- sample_groups$start[i] + sample_groups$Count[i] - 1
}

#create data frame for comparisons
comp_groups <- data.frame(seq(from=1, to=comp_number))
colnames(comp_groups) <- "CompNumber"
for(i in 1:comp_number){
  comp_N <- grep("N", sample_groups[[i+1]], ignore.case=TRUE)
  comp_D <- grep("D", sample_groups[[i+1]], ignore.case=TRUE)
  comp_groups$comp_N[i] <- comp_N 
  comp_groups$comp_D[i] <- comp_D
  comp_groups$N_start[i] <- sample_groups$start[comp_N]
  comp_groups$N_end[i] <- sample_groups$end[comp_N]
  comp_groups$D_start[i] <- sample_groups$start[comp_D]
  comp_groups$D_end[i] <- sample_groups$end[comp_D]
  comp_groups$comp_name[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D])
  comp_groups$fc[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_FC")
  comp_groups$fc2[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_FC2")
  comp_groups$pval[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_Pval")
  comp_groups$limma_pval[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_LimmaPval")
  comp_groups$Ncount <- sample_groups$end[comp_N] - sample_groups$start[comp_N] +1
  comp_groups$Dcount <- sample_groups$end[comp_D] - sample_groups$start[comp_D] +1
  }

testme <- rep(1,3)

excel_order <- sample_info$PD_Order
#sample_list <- sample_info$ID
group_list <- sample_groups$Group
treatment_groups <- sample_info$Group

group_color <- color_choices[1:group_number]
sample_groups$colorlist <- color_choices[1:group_number]
sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info$Group, Group)])
sample_groups$title <- str_c(sample_groups$Group,"(",sample_groups$colorlist,")")


color_list<- sample_info$colorlist

group_title <- paste(sample_groups$title, collapse=" ")
group_title <- gsub("\\)", "\\),",  group_title)


#organize column headers for final output
sample_info$Header1 <- str_c(sample_info$ID, " ", sample_info$Group)
sample_info$Header2 <- str_c(sample_info$ID, " ", sample_info$Group, " Normalized")


#Create DEP summarized experiment
dep_se <- data.frame(cbind(sample_info$Label, sample_info$Group, sample_info$Replicate))
colnames(dep_se) <- c("label", "condition", "replicate")
dep_se$label <- as.character(dep_se$label)
dep_se$condition <- as.character(dep_se$condition)
dep_se$replicate <- as.numeric(dep_se$replicate)


#create dataframe to hold cv summaries for normalization strategies
summary_cv <- data.frame(sample_groups$Group)

# create subdirectory to store csv and plots
data_dir <- file_prefix
output_dir <- create_dir(data_dir)
file_prefix1 <- str_c(output_dir, file_prefix)

output_dir2 <- create_dir(str_c(data_dir,"//Backup"))

output_dir3 <- create_dir(str_c(data_dir,"//Extra"))
file_prefix3 <- str_c(output_dir3, file_prefix)

file.copy(input_data, str_c(output_dir2, basename(input_data)))
file.copy(study_design, str_c(output_dir2, basename(study_design)))
file.copy("Quant Setup_v6.R", str_c(output_dir2, "Quant Setup_v6.R"))
file.copy("Quant Main_v6.R", str_c(output_dir2, "Quant Main_v6.R"))
file.copy("Quant Functions Misc v6.R", str_c(output_dir2, "Quant Functions Misc v6.R"))
file.copy("Quant Functions Impute v6.R", str_c(output_dir2, "Quant Functions Impute v6.R"))
file.copy("Quant Functions Stats v6.R", str_c(output_dir2, "Quant Functions Stats v6.R"))
file.copy("Quant Functions Plots v6.R", str_c(output_dir2, "Quant Functions Plots v6.R"))
file.copy("Quant Functions Norm v6.R", str_c(output_dir2, "Quant Functions Norm v6.R"))

#rmarkdown::render("Quant Main_v5.R")
