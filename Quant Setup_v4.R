library(stringr)
library(readxl)
library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(gridExtra)


#-------------------------------------------
# User input
#-------------------------------------------

# excel export from PD
forward_data <- read_excel("5116 protein peptide 091118.xlsx", 1)
# excel file contains 4 columns, Output order from PD (Protein/peptide order on export, PSM numerical order), ID, Label, Group
sample_info <- read_excel("5116_SampleList_091118.xlsx", 1)
# excel export from PD for decoy data
# decoy_data <- read_excel("4983_MS2_decoyPSM_041918.xlsx", 1)
# file prefix for excel outputs
file_prefix <- "5116_091118"


psm_input <- FALSE
psm_to_peptide <- FALSE

protein_peptide_input <- TRUE #PD export with nested protein/peptide, need specific columns
peptide_to_protein <- TRUE  #collapse from

normalize_to_protein <- TRUE  #set accession numbers below for list of target proteins
pair_comp <- FALSE  # pairwise comparison

log_normalize <- TRUE #log2 intensitites prior to normalization, will unlog after to report unlog intensities


holes <- "Impute"  # "Impute", Average", "Minimum", "Floor"
intensity_cutoff <- 5000000  # if a replicate group has >50% missing values and a measured value is above this cutoff then the measured value is a misalignment, value will be removed
area_floor <- 1000


pvalue_cutoff <- 0.05  # extra worksheet is created in excel for each comparison, subset with these cutoff values
fc_cutoff <- 2

color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black")


# Protein accession numbers for plots and normalization
adh_list <- c("P00330")
bait_list <- c("P02662","P02663")
avidin_list <- c("P02701")
#avidin_list <- c("P22629")  #steptavidin
casein_list <- c("P02662", "P02663")
#carbox_list <- c("Q05920", "Q91ZA3", "Q5SWU9") #mouse
carbox_list <- c("F1QPL7", "A0A0R4IFJ4", "F1QM37") #zebrafish
#bira_list <- c("P06709")
bira_list <- c("O66837") #biotin ligase

protein_norm_list <- c("F1QPL7", "A0A0R4IFJ4", "F1QM37") #use these if normalize_to_protein <- TRUE



source("Quant Functions_v4.R")

source("Quant Main_v4.R")
