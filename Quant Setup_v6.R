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
msgBox <- tkmessageBox(title = "Design",
                       message = "Load study design file", icon = "info", type = "ok")
study_design <- file.choose()


#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


#color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black")
#color_choices <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#color_choices <- c("dodgerblue", "firebrick", "green", "gold", "grey", "bisque", "cyan", "slategrey", "black", "pink", "coral", "maroon1")
#color_choices <- rainbow(10)

#color_choices <- rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
#color_choices <- heat.colors(n, alpha = 1)
#color_choices <- terrain.colors(n, alpha = 1)
#color_choices <- topo.colors(n, alpha = 1)
#color_choices <- cm.colors(n, alpha = 1)
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
data_input_type <- as.character(design_info[3,2])
data_output_type <- as.character(design_info[4,2])
phos_peptide_only <- as.logical(design_info[5,2])

protein_peptide_input <- FALSE
peptide_to_protein <- FALSE
psm_peptide_fdr <- FALSE
tmt_spqc_sets <- FALSE



if ((data_input_type == "Protein Peptide") & (data_output_type == "Protein")) {
    protein_peptide_input <- TRUE
    peptide_to_protein <- TRUE
  }else if(data_input_type == "Protein" & data_output_type == "Protein") {
    protein_peptide_input <- FALSE
    peptide_to_protein <- FALSE
  }else if(data_input_type == "Peptide" & data_output_type == "Peptide") {
    protein_peptide_input <- FALSE
    peptide_to_protein <- FALSE
  }else if(data_input_type == "Peptide" & data_output_type == "Protein") {
    protein_peptide_input <- FALSE
    peptide_to_protein <- TRUE
  }else if (data_input_type == "PSM Peptide Decoy" & data_output_type == "Peptide") {
    psm_peptide_fdr <- TRUE
   protein_peptide_input <- FALSE
    peptide_to_protein <- FALSE
  }else if(data_input_type == "PSM Peptide Decoy" & data_output_type == "Protein") {
    psm_peptide_fdr <- TRUE
    protein_peptide_input <-FALSE
    peptide_to_protein <- TRUE
  }else if (data_input_type == "TMT SPQC Protein Peptide" & data_output_type == "Protein") {
    tmt_spqc_sets <- TRUE
    protein_peptide_input <- TRUE
    peptide_to_protein  <- FALSE
  }else if(data_input_type == "TMT SPQC Protein" & data_output_type == "Protein") {
    tmt_spqc_sets <- TRUE
    protein_peptide_input <- FALSE
    peptide_to_protein  <-FALSE
  }else{
    msgBox <- tkmessageBox(title = "Whoops",
                           message = "Invalid input output in design file", icon = "info", type = "ok")
  }



normalize_to_protein <- as.logical(design_info[7,2])
pair_comp <- as.logical(design_info[9,2])
adh_spike <- as.logical(design_info[10,2])
impute_method <- as.character(design_info[11,2])
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
protein1_list <- as.character(na.omit(t(design_info[22,2:11])))
protein2_list <- as.character(na.omit(t(design_info[23,2:11])))
protein3_list <- as.character(na.omit(t(design_info[24,2:11])))
protein4_list <- as.character(na.omit(t(design_info[25,2:11])))
protein_norm_list <- as.character(na.omit(t(design_info[8,2:11])))

misaligned_filter <- as.logical(design_info[26,2])
#impute_method <- as.character(design_info[27,2])

use_ti <-  as.logical(design_info[28,2])
use_ai <-  as.logical(design_info[29,2])
use_mi <-  as.logical(design_info[30,2])
use_lr <-  as.logical(design_info[31,2])
use_vsn <-  as.logical(design_info[32,2])
use_quantile <-  as.logical(design_info[33,2])
use_loess <-  as.logical(design_info[34,2])
use_localsl <-  as.logical(design_info[35,2])



if (psm_peptide_fdr == TRUE){
  msgBox <- tkmessageBox(title = "Data",
                         message = "Load forward psm data", icon = "info", type = "ok")
  forward_psm_name <- file.choose()
  msgBox <- tkmessageBox(title = "Data",
                         message = "Load decoy psm data", icon = "info", type = "ok")
  decoy_psm_name <- file.choose()  
  msgBox <- tkmessageBox(title = "Data",
                         message = "Load forward peptide data", icon = "info", type = "ok")
  forward_peptide_name <- file.choose()  
  msgBox <- tkmessageBox(title = "Data",
                         message = "Load decoy peptide data", icon = "info", type = "ok")
  decoy_peptide_name <- file.choose() 
  forward_psm <- read_excel(forward_psm_name, 1)
  decoy_psm <- read_excel(decoy_psm_name, 1)
  forward_peptide <- read_excel(forward_peptide_name, 1)
  decoy_peptide <- read_excel(decoy_peptide_name, 1)
}else{
  msgBox <- tkmessageBox(title = "Data",
                         message = "Load raw data file", icon = "info", type = "ok")
  input_data <- file.choose()
  forward_data <- read_excel(input_data, 1)
}



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
  comp_groups$exactTest[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_ExactTest")  
  comp_groups$Ncount <- sample_groups$end[comp_N] - sample_groups$start[comp_N] +1
  comp_groups$Dcount <- sample_groups$end[comp_D] - sample_groups$start[comp_D] +1
  }

testme <- rep(1,3)

excel_order <- sample_info$PD_Order
#sample_list <- sample_info$ID
group_list <- sample_groups$Group
treatment_groups <- sample_info$Group

color_choices <- rainbow(group_number, s = 1, v = 1, start = 0, end = max(1, group_number - 1)/group_number, alpha = 1)
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

#create factor vector
group_factor <- rep(1, sample_groups$Count[1])
for(i in 2:nrow(sample_groups)) {
  group_factor <- c(group_factor, rep(i, sample_groups$Count[i]))
}
group_factor <-factor(group_factor)


#Create DEP summarized experiment
dep_se <- data.frame(cbind(sample_info$Label, sample_info$Group, sample_info$Replicate))
colnames(dep_se) <- c("label", "condition", "replicate")
dep_se$label <- as.character(dep_se$label)
dep_se$condition <- as.character(dep_se$condition)
dep_se$replicate <- as.numeric(dep_se$replicate)


#create dataframe to hold cv summaries for normalization strategies
summary_cv <- data.frame(sample_groups$Group)

#create shiny list of norms used and proteins of interest
shiny_norm <- shiny_norm_list()
shiny_norm_final <- shiny_norm_list_final()
shiny_protein <- c(adh_list, bait_list, avidin_list, casein_list, carbox_list, bira_list, 
                   protein1_list, protein2_list, protein3_list, protein4_list, protein_norm_list)
shiny_protein_name <- c("ADH_", "Bait_", "Avidin_", "Casein_", "Carbox_", "BirA_", "User_", "User_", "User_", "User_", "Norm_")
names(shiny_protein) <- str_c(shiny_protein_name, shiny_protein)
shiny_comp <- comp_groups$comp_name
names(shiny_comp) <- comp_groups$comp_name


# create subdirectory to store csv and plots
data_dir <- file_prefix
output_dir <- create_dir(data_dir)
file_prefix1 <- str_c(output_dir, file_prefix)

output_dir2 <- create_dir(str_c(data_dir,"//Backup"))

output_dir3 <- create_dir(str_c(data_dir,"//Extra"))
file_prefix3 <- str_c(output_dir3, file_prefix)

if (psm_peptide_fdr == FALSE){file.copy(input_data, str_c(output_dir2, basename(input_data)))}
file.copy(study_design, str_c(output_dir2, basename(study_design)))
file.copy("Quant Setup_v6.R", str_c(output_dir2, "Quant Setup_v6.R"))
file.copy("Quant Main_v6.R", str_c(output_dir2, "Quant Main_v6.R"))
file.copy("Quant Functions Misc v6.R", str_c(output_dir2, "Quant Functions Misc v6.R"))
file.copy("Quant Functions Impute v6.R", str_c(output_dir2, "Quant Functions Impute v6.R"))
file.copy("Quant Functions Stats v6.R", str_c(output_dir2, "Quant Functions Stats v6.R"))
file.copy("Quant Functions Plots v6.R", str_c(output_dir2, "Quant Functions Plots v6.R"))
file.copy("Quant Functions Norm v6.R", str_c(output_dir2, "Quant Functions Norm v6.R"))
app_dir <- create_dir(str_c(data_dir,"//Backup//Quant"))
file.copy("Quant//app.R", str_c(app_dir, "app.R"))
#rmarkdown::render("Quant Main_v5.R")

