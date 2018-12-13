rm(list = ls())
#-------------------------------------------
# User input
#-------------------------------------------

# excel export from PD
input_data <- "4227_Striatum_TMT_SL_IRS_Norm.xlsx"
input_sample_info <-"4227 TMT Striatum Combined 120618 Sample Info2.xlsx"

# excel export from PD for decoy data
# decoy_data <- read_excel("4983_MS2_decoyPSM_041918.xlsx", 1)
# file prefix for excel outputs
file_prefix <- "4227_Cortex"

psm_input <- FALSE
psm_to_peptide <- FALSE

protein_peptide_input <- FALSE #PD export with nested protein/peptide, need specific columns
peptide_to_protein <- FALSE  #collapse from peptide back to protein for final output

normalize_to_protein <- FALSE  #set accession numbers below for list of target proteins
pair_comp <- FALSE  # pairwise comparison

adh_spike <- FALSE  #add adh plots if spiked into sample

holes <- "Floor"  # "Impute", Average", "Minimum", "Floor"
intensity_cutoff <- 5000000  # if a replicate group has >50% missing values and a measured value is above this cutoff then the measured value is a misalignment, value will be removed
area_floor <- 0.1


pvalue_cutoff <- 0.05  # extra worksheet is created in excel for each comparison, subset with these cutoff values
fc_cutoff <- 1.5

#color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black")
color_choices <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Protein accession numbers for plots and normalization
adh_list <- c("P00330")
bait_list <- c("Q8K3K8")
avidin_list <- c("P02701")
#avidin_list <- c("P22629")  #steptavidin
casein_list <- c("P02662", "P02663")
carbox_list <- c("Q05920", "Q91ZA3", "Q5SWU9") #mouse
#carbox_list <- c("F1QPL7", "A0A0R4IFJ4", "F1QM37") #zebrafish
#bira_list <- c("P06709")
bira_list <- c("O66837") #biotin ligase 
ko1_list <- c("F6SEU4")
ko2_list <- c("O08759")
ko3_list <- c("Q80Z38")
ko4_list <- c("Q4ACU6")
#protein_norm_list <- c("F1QPL7", "A0A0R4IFJ4", "F1QM37") #use these if normalize_to_protein <- TRUE
protein_norm_list <- c("O66837")

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

options(digits=3)

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

source("Quant Functions Misc v6.R")
source("Quant Functions Norm v6.R")
source("Quant Functions Impute v6.R")
source("Quant Functions Stats v6.R")
source("Quant Functions Plots v6.R")

forward_data <- read_excel(input_data, 1)
sample_info <- read_excel(input_sample_info, 1)

# number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file
sample_number <- nrow(sample_info)
comp_number <- length(grep(x = colnames(sample_info), pattern = "^Comp"))
#count number of samples in each group - add column
sample_info$Count <- sapply(sample_info$Group, function(string) sum(string==sample_info$Group))

# create unqiue group dataframe with sample number
sample_groups <- sample_info[5:ncol(sample_info)]
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


#create dataframe to hold cv summaries for normalization strategies
summary_cv <- data.frame(sample_groups$Group)

# create subdirectory to store csv and plots
if(dir.exists(file.path(".", "output_files"))) { unlink(file.path(".", "output_files"), recursive = TRUE, force=TRUE)}
ifelse(!dir.exists(file.path(".", "output_files")), dir.create(file.path(".", "output_files")), FALSE)
output_dir <- ".//output_files//"
file_prefix1 <- str_c(output_dir, file_prefix)

if(dir.exists(file.path(".", "output_files//Backup"))) { unlink(file.path(".", "output_files//Backup"), recursive = TRUE, force=TRUE)}
ifelse(!dir.exists(file.path(".", "output_files//Backup")), dir.create(file.path(".", "output_files//Backup")), FALSE)
output_dir2 <- ".//output_files//Backup//"

if(dir.exists(file.path(".", "output_files//Extra"))) { unlink(file.path(".", "output_files//Extra"), recursive = TRUE, force=TRUE)}
ifelse(!dir.exists(file.path(".", "output_files//Extra")), dir.create(file.path(".", "output_files//Extra")), FALSE)
output_dir3 <- ".//output_files//Extra//"
file_prefix3 <- str_c(output_dir3, file_prefix)

file.copy(input_data, str_c(output_dir2, input_data))
file.copy(input_sample_info, str_c(output_dir2, input_sample_info))
file.copy("Quant Setup_v6.R", str_c(output_dir2, "Quant Setup_v6.R"))
file.copy("Quant Main_v6.R", str_c(output_dir2, "Quant Main_v6.R"))
file.copy("Quant Functions Misc v6.R", str_c(output_dir2, "Quant Functions Misc v6.R"))
file.copy("Quant Functions Impute v6.R", str_c(output_dir2, "Quant Functions Impute v6.R"))
file.copy("Quant Functions Stats v6.R", str_c(output_dir2, "Quant Functions Stats v6.R"))
file.copy("Quant Functions Plots v6.R", str_c(output_dir2, "Quant Functions Plots v6.R"))
file.copy("Quant Functions Norm v6.R", str_c(output_dir2, "Quant Functions Norm v6.R"))

#rmarkdown::render("Quant Main_v5.R")
