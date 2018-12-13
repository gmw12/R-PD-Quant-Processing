library(stringr)
library(readxl)
library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(gridExtra)
library(MASS)


# create final excel documents
excel_gw <- function(df, filename) {
  require(openxlsx)
  filename <- str_c(file_prefix, filename)
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet=1, df)  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# input PD protein/peptide export and output peptide (keeps master protein assignment)
protein_to_peptide <- function(data_in){
  data_in[,1][is.na(data_in[,1])] <- ""  
  
  for(i in 1:nrow(data_in)) {
    if(data_in[i,1] == "High") {
      set_accession <- data_in[i,3]
      set_description <- data_in[i,4]
    }else{
      data_in[i,3] <- set_accession
      data_in[i,4] <- set_description
    }
  }
  peptide_header <- data_in[2,]
  data_in <- subset(data_in, Master=="High")
  colnames(data_in) <- peptide_header
  colnames(data_in)[3] <- "Accession"
  colnames(data_in)[4] <- "Description"
  data_in <- data_in[,-1]
  colnames(data_in)[colnames(data_in) == 'Quan Usage'] <- 'Used'
  data_in <- subset(data_in, Used=="Used")
  while (sapply("Abundance", grepl, colnames(data_in[ ,ncol(data_in)]))==FALSE){
    data_in <- data_in[, -ncol(data_in)]
  }
  data_in[, (ncol(data_in)-sample_number+1):ncol(data_in)] <- 
    sapply(data_in[, (ncol(data_in)-sample_number+1):ncol(data_in)], as.numeric)
  return(data_in)
}


# create column to flag and delete rows with no data
delete_empty_row <- function(data_in){
  data_in$na_count <- apply(data_in, 1, function(x) sum(is.na(x)))
  data_in <- subset(data_in, data_in$na_count != tmt_kit)
  data_in <- data_in[1:(info_count + tmt_kit)]
  return(data_in)
}


#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide <- function(peptide_data){
  colnames(peptide_data)[5] <- "Peptides"
  peptide_data$Peptides <- 1
  peptide_data$Peptides <- as.numeric(peptide_data$Peptides)
  test1 <- peptide_data[ , -1]
  test1 <- test1[ , -3]
  colnames(test1)[colnames(test1) == 'Master Protein Accessions'] <- 'Accession'
  colnames(test1)[colnames(test1) == 'Accessions'] <- 'Accession'
  test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(funs(sum))
  return(test2)
}


#Bar plot-------------------------------------------------
barplot_gw <- function(x,y) {
  x[is.na(x)] <- 0.0
  x <- colSums(x)
  png(filename=str_c(output_dir, y, "_barplot.png"), width = 888, height = 571)  
  barplot(x, 
          col = color_list,
          names.arg = treatment_groups,
          cex.names = .8,
          las=2,
          main = y)
  dev.off()
}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
input_data <- "4227 Cortex Combined ProteinPeptideMouse_110918.xlsx"
input_sample_info <- "4227 TMT Cortex Combined 110718 Sample Info.xlsx"

all_data <- read_excel(input_data, 1)
sample_info <- read_excel(input_sample_info, 1)
sets_tmt <- 4
tmt_kit <- 11
sd_factor <- 2  
sample_number <- nrow(sample_info)
file_prefix <- "4227_Cortex"
treatment_groups <- sample_info$Group
color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black", "dodgerblue", "gold", "lightblue", "palegreen")

# create unqiue group dataframe with sample number
sample_groups <- sample_info[5:ncol(sample_info)]
sample_groups <- sample_groups %>% distinct(Group, .keep_all = TRUE)
group_number <- nrow(sample_groups)
group_color <- color_choices[1:group_number]
sample_groups$colorlist <- color_choices[1:group_number]
sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info$Group, Group)])
sample_groups$title <- str_c(sample_groups$Group,"(",sample_groups$colorlist,")")
color_list<- sample_info$colorlist

# create subdirectory to store csv and plots
if(dir.exists(file.path(".", "outputTMT_files"))) { unlink(file.path(".", "outputTMT_files"), recursive = TRUE, force=TRUE)}
ifelse(!dir.exists(file.path(".", "outputTMT_files")), dir.create(file.path(".", "outputTMT_files")), FALSE)
output_dir <- ".//outputTMT_files//"
file_prefix <- str_c(output_dir, file_prefix, "_")


file.copy(input_data, str_c(output_dir, input_data))
file.copy(input_sample_info, str_c(output_dir, input_sample_info))
file.copy("TMT IRS_v3.R", str_c(output_dir, "TMT IRS_v3.R"))

#break down PD output to peptide only with protein assignments
all_peptide <- protein_to_peptide(all_data)
excel_order <- sample_info$PD_Order
info_count <- ncol(all_peptide) - sample_number
info_columns <- all_peptide[1:info_count]
all_peptide <-  all_peptide[(info_count+1):ncol(all_peptide)][, (excel_order)]
barplot_gw(all_peptide, "Raw_Peptide")
raw_peptide <- cbind(info_columns, all_peptide)
excel_gw(raw_peptide, "TMT_raw_peptide.xlsx")


# split combined data from PD into individual sets
start_col<-1
for(i in 1:sets_tmt){
  irs_set <- subset(sample_info, sample_info$Set==i)
  irs_set$setorder <- seq(start_col, (nrow(irs_set)+start_col-1))
  irs_set$sequence <- seq(1, nrow(irs_set))
  temp_data <- cbind(info_columns, all_peptide[ ,irs_set$setorder])
  temp_data <- delete_empty_row(temp_data)
  temp_data[is.na(temp_data)] <- 0.0
  assign(str_c("Raw",i,"_peptide"), temp_data)
  assign(str_c("IRS",i,"_set"), irs_set)
  treatment_groups<-irs_set$Group
  color_list<-irs_set$colorlist
  barplot_gw(temp_data[(info_count+1):ncol(temp_data)], str_c("Raw_Peptide_",i))
  start_col <- start_col + tmt_kit
}


#identify common proteins #1
common_protein <- unlist(data.frame(get("Raw1_peptide")[ , "Accession"]))
for(i in 2:sets_tmt){
  common_protein2 <- unlist(data.frame(get(str_c("Raw",i,"_peptide"))[ , "Accession"]))
  common_protein <- intersect(common_protein, common_protein2)
}

for(i in 1:sets_tmt){
  temp_peptide <- get(str_c("Raw",i,"_peptide"))
  temp_peptide <- subset(temp_peptide, temp_peptide$Accession %in% common_protein)
  temp_peptide <- temp_peptide[order(temp_peptide$Accession),]
  barplot_gw(temp_peptide[(ncol(temp_peptide)-tmt_kit+1):ncol(temp_peptide)], str_c("CommonProtein_Filtered_",i))
  assign(str_c("Raw",i,"_common"), subset(temp_peptide, temp_peptide$Accession %in% common_protein))
}


# SL normalize sets, using target created from all TMT sets
temp_data <- get(str_c("Raw",1,"_common"))
target <- colSums(temp_data[(info_count+1):ncol(temp_data)])
for(i in 2:sets_tmt){
  temp_data <- get(str_c("Raw",i,"_common"))
  target <- c(target, colSums(temp_data[(info_count+1):ncol(temp_data)]))
}
target <- mean(target)
for(i in 1:sets_tmt){
  temp_data <- get(str_c("Raw",i,"_common"))
  irs_info <- temp_data[1:info_count]
  temp_data <- temp_data[(info_count+1):ncol(temp_data)]
  norm_facs <- target / colSums(temp_data)
  temp_sl <- sweep(temp_data, 2, norm_facs, FUN = "*")
  color_list<-irs_set$colorlist
  barplot_gw(temp_sl, str_c("SL_Peptide_",i))
  temp_sl <- cbind(irs_info, temp_sl)
  assign(str_c("Norm",i,"_SL_peptide"), temp_sl)
  excel_gw(temp_sl, str_c("Norm",i, "_SL_peptide.xlsx"))
}

#Filter peptides based on CV of SPQC
for(i in 1:sets_tmt){
  temp_data <- get(str_c("Norm",i,"_SL_peptide"))
  irs_set <- get(str_c("IRS",i,"_set"))
  spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
  spqc_count <- nrow(spqc_set)
  temp_data$average <- apply(temp_data[(info_count+1):(info_count+spqc_count)], 1, FUN = mean)
  temp_data$stddev <- apply(temp_data[(info_count+1):(info_count+spqc_count)], 1, FUN = sd)
  temp_data$cv <- temp_data$stddev / temp_data$average * 100
  temp_data$bin <- ntile(temp_data$average, 5)  
  sd_info <- subset(temp_data) %>% group_by(bin) %>% 
    summarize(min = min(average), max = max(average), min_cv= min(cv), max_cv = max(cv), av_cv = mean(cv), sd_cv = sd(cv))
  temp_data$maxcv <- NA  
  for (j in 1:nrow(temp_data)){
    temp_data$maxcv[j] <- sd_info$av_cv[temp_data$bin[j]] + (sd_factor * sd_info$sd_cv[temp_data$bin[j]])
  }
  temp_data <- subset(temp_data, temp_data$cv < temp_data$maxcv)
  temp_data <- temp_data[1:(info_count+tmt_kit)]
  color_list<-irs_set$colorlist
  barplot_gw(temp_data[(info_count+1):ncol(temp_data)], str_c("Norm_Peptide_Filtered_",i))
  assign(str_c("Norm",i,"_SL_filtered_peptide"), temp_data)
  excel_gw(temp_data, str_c('Norm_filtered',i,'_peptide.xlsx'))
  #collapse to protein level
  temp_protein <- collapse_peptide(temp_data)
  #barplot_gw(temp_protein[], str_c("SL_Protein_Filtered_",i))
  assign(str_c("IRS",i,"_protein"), temp_protein)
  barplot_gw(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_Protein_Filtered_",i))
  excel_gw(temp_protein, str_c("SL_Filtered",i, "_protein.xlsx"))
  }


#identify common proteins
common_protein <- unlist(data.frame(get("IRS1_protein")[ , "Accession"]))
for(i in 2:sets_tmt){
  common_protein2 <- unlist(data.frame(get(str_c("IRS",i,"_protein"))[ , "Accession"]))
  common_protein <- intersect(common_protein, common_protein2)
}

for(i in 1:sets_tmt){
  temp_protein <- get(str_c("IRS",i,"_protein"))
  temp_protein <- subset(temp_protein, temp_protein$Accession %in% common_protein)
  temp_protein <- temp_protein[order(temp_protein$Accession),]
  barplot_gw(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_CommonProtein_Filtered_",i))
  assign(str_c("IRS",i,"_common"), subset(temp_protein, temp_protein$Accession %in% common_protein))
}

#reset info columns for protein level data
irs_info <- data.frame(IRS1_common[1:(ncol(IRS1_common)-tmt_kit)])
info_count <- ncol(irs_info)

#create SPQC data frames for each set, with average
for(i in 1:sets_tmt){
  spqc_data <- get(str_c("IRS",i,"_common"))
  irs_set <- get(str_c("IRS",i,"_set"))
  spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
  spqc_set$sequence <- spqc_set$sequence + info_count
  spqc_data <- spqc_data[, (spqc_set$sequence)]
  spqc_data$average <- rowMeans(spqc_data)
  assign(str_c("SPQC",i,"_data"), spqc_data)
  }

#create dataframe with SPQC averages
spqc_all <- data.frame(SPQC1_data$average)
for(i in 2:sets_tmt){
  spqc_data <- get(str_c("SPQC",i,"_data"))
  spqc_all <- cbind(spqc_all, spqc_data$average)
}
spqc_all$average <- rowMeans(spqc_all)

for(i in 1:sets_tmt){
  spqc_data <- get(str_c("SPQC",i,"_data"))
  spqc_data$fact <- spqc_data$average/spqc_all$average
  assign(str_c("SPQC",i,"_data"), spqc_data)
}

#apply IRS norm to each set
for(i in 1:sets_tmt){
  temp_data <- get(str_c('IRS',i,'_common'))[(info_count+1):(info_count+tmt_kit)]
  spqc_data <- get(str_c('SPQC',i,'_data'))
  temp_data <- temp_data/spqc_data$fact
  barplot_gw(temp_data[(ncol(temp_data)-tmt_kit+1):ncol(temp_data)], str_c("SL_IRS_CommonProtein_Filtered_",i))
  assign(str_c("IRS",i,"_final"), temp_data)
  excel_gw(temp_data, str_c("IRS_SL_Filtered",i, "_protein.xlsx"))
}

#recombine sets for output
final_SL_IRS <- IRS1_final
for(i in 2:sets_tmt){
  final_SL_IRS <- cbind(final_SL_IRS, get(str_c("IRS",i,"_final")))
}
#reorder combined set in original order
set_original <- IRS1_set
for(i in 2:sets_tmt){
  irs_set <- get(str_c('IRS',i,'_set'))
  set_original <- rbind(set_original, irs_set) 
}
set_original$reorder <- seq(1,(tmt_kit*sets_tmt))
set_original <- set_original[order(set_original$PD_Order),]
final_SL_IRS <- final_SL_IRS[ ,set_original$reorder]
final_SL_IRS <- cbind(irs_info, final_SL_IRS)

excel_gw(final_SL_IRS, "TMT_SL_IRS_Norm.xlsx")

treatment_groups<-sample_info$Group
color_list<-sample_info$colorlist
barplot_gw(final_SL_IRS[(info_count+1):ncol(final_SL_IRS)], "Final_SL_IRS_CommonProtein_Filtered_")






