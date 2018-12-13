#global options, numbers, sig digits
options(scipen = 999)
options(digits=3)

# !diagnostics off

# number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file

sample_columns <- ncol(sample_info)
sample_number <- nrow(sample_info)
comp_number <- sample_columns - 4 #subtract 3 non comparison columns

#count number of samples in each group - add column
sample_info$Count <- sapply(sample_info$Group, function(string) sum(string==sample_info$Group))

# create unqiue group dataframe with sample number
sample_groups <- sample_info[4:(sample_columns+1)]
sample_groups <- sample_groups %>% distinct(Group, .keep_all = TRUE)

# comparison groups in pairs, c(3,2) 3/2, c(3,2,4,2) 3/2, 4/2
group_comp <- vector(mode="numeric", length=0)
for (i in 2:(comp_number+1)){
  group_comp <- c(group_comp, grep("N", sample_groups[[i]], ignore.case=TRUE))
  group_comp <- c(group_comp, grep("D", sample_groups[[i]], ignore.case=TRUE))
}

group_number <- nrow(sample_groups)
excel_order <- sample_info$PD_Order
sample_list <- sample_info$ID
group_list <- sample_groups$Group
group_rep <- sample_groups$Count
treatment_groups <- sample_info$Group
group_color <- color_choices[1:group_number]
sample_groups$colorlist <- color_choices[1:group_number]
sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info$Group, Group)])
sample_groups$title <- str_c(sample_groups$Group,"(",sample_groups$colorlist,")")

# create a list of groups to be compared, fold change (inverse for neg), foldchange, pval
comp_fc_groups <- NULL
comp_fc2_groups <- NULL
comp_pval_groups <- NULL
comp_groups <- NULL
comp_header <- NULL
z <- 1
for(i in 1:comp_number) {
  fc1 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC", collapse = " ")
  fc2 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC2", collapse = " ")
  pval1 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_pval", collapse = " ")
  comp <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]], collapse = " ")
  comp_fc_groups <- c(comp_fc_groups, fc1 )
  comp_fc2_groups <- c(comp_fc2_groups, fc2)
  comp_pval_groups <- c(comp_pval_groups, pval1)
  comp_groups <- c(comp_groups, comp)
  comp_header <- c(comp_header, fc1, fc2, pval1)
  z <- z+2
}


sample_groups$cv <- str_c(sample_groups$Group, "_CV")  
sample_groups$log <- str_c(sample_groups$Group, "_log")
#need to replace group_cv with sample_groups$cv and log
group_cv <- sample_groups$cv
group_log <- sample_groups$log

color_list<- sample_info$colorlist

group_title <- paste(sample_groups$title, collapse=" ")
group_title <- gsub("\\)", "\\),",  group_title)

#assign start and end columns for each group in pd output
group_startcol <- 1
for(i in 1:(group_number-1)) group_startcol <- c(group_startcol, group_startcol[i]+group_rep[i])
group_endcol <- group_startcol + group_rep
group_endcol <- group_endcol -1


#organize column headers for final output
sample_header <- str_c(sample_list[1]," ", treatment_groups[1])
for(i in 2:sample_number){
  sample_header <- c(sample_header, str_c(sample_list[i], " ", treatment_groups[i]))
}

for(i in 1:sample_number){
  sample_header <- c(sample_header, str_c(sample_list[i], " ", treatment_groups[i], " Normalized"))
}

sample_header <- c(sample_header, group_cv, comp_header)


# create subdirectory to store csv and plots
if(dir.exists(file.path(".", "output_files"))) { unlink(file.path(".", "output_files"), recursive = TRUE, force=TRUE)}

ifelse(!dir.exists(file.path(".", "output_files")), dir.create(file.path(".", "output_files")), FALSE)
output_dir <- ".//output_files//"
file_prefix <- str_c(output_dir, file_prefix)


#format data from nested protein/peptide output to peptide only
if (protein_peptide_input){
  forward_data$`Protein FDR Confidence: Combined`[is.na(forward_data$`Protein FDR Confidence: Combined`)] <- ""  
  for(i in 1:nrow(forward_data)) {
    if(forward_data[i,1] == "High") {set_accession <- forward_data[i,3]
    }else{forward_data[i,3] <- set_accession}
  }
  
  protein_names <- forward_data
  colnames(protein_names)[1] <- "filter"
  protein_names <- subset(protein_names, filter=="High")
  protein_names <- protein_names[ , c(3:4)]
  peptide_header <- forward_data[2,]
  forward_data <- subset(forward_data, Master=="High")
  colnames(forward_data) <- peptide_header
  forward_data <- forward_data[,-1]
  colnames(forward_data)[colnames(forward_data) == 'Quan Usage'] <- 'Used'
  forward_data <- subset(forward_data, Used=="Used")
  forward_data <- forward_data[-ncol(forward_data)]
  forward_data[(ncol(forward_data)-sample_number+1):ncol(forward_data)] <- as.numeric(unlist(forward_data[(ncol(forward_data)-sample_number+1):ncol(forward_data)]))
  colnames(forward_data)[2] <- "Master Protein Accessions"
}


#----- edit column headers
col_headers <- colnames(forward_data) 
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Combined", "Confidence")
col_headers <- str_replace(col_headers, "Annotated Sequence", "Annotated_Sequence")
col_headers <- str_replace(col_headers, "Master Protein Accessions", "Accessions")
col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
col_headers <- str_replace(col_headers,"\\[", "")
col_headers <- str_replace(col_headers,"\\]", "")
colnames(forward_data) <- col_headers
total_columns <- ncol(forward_data)
info_columns <- total_columns - sample_number
info_headers <- colnames(forward_data[1:info_columns])
final_sample_header <- c(info_headers, sample_header)





if (psm_input) {
  data_ready <- psm_decoy(forward_data, decoy_data)
  }else{
    data_ready <- forward_data 
  }  

# create column to flag and delete rows with no data
data_ready$na_count <- apply(data_ready, 1, function(x) sum(is.na(x)))
data_ready <- subset(data_ready, data_ready$na_count < (sample_number))
data_ready <- data_ready[1:total_columns]

# save the annotation columns for later and remove from data frame
annotate_df <- data_ready[1:info_columns]
data_ready <- data_ready[(info_columns+1):total_columns]
#row.names(data_ready) <- annotate_df$`Accessions` delete later if no errors 


#reorder data if needed, will use PD_order from sample list
data_ready <- order_columns(data_ready)


# create histogram of all measured values in data set, compare against hole fill confidence value
histogram_gw(data_ready,"Total_Set_Intensity","Log2 Intensity Distribution")

#save copy of raw data ready, if protein/peptide input will be peptide data
write.csv(cbind(annotate_df, data_ready), file= str_c(file_prefix, "_Raw.csv", collapse = " "))


#replaces NA, no data, with assigned value
if (holes == "Floor") {
  data_ready[is.na(data_ready)] <- area_floor
  write.csv(cbind(annotate_df, data_ready), file= str_c(file_prefix, "_Floor.csv", collapse = " "))
} else if (holes == "Average") {
  data_ready <- hole_average(data_ready)
  write.csv(cbind(annotate_df, data_ready), file= str_c(file_prefix, "_Average.csv", collapse = " "))
} else if (holes == "Minimium") {
  data_ready <- hole_minimum(data_ready)
  write.csv(cbind(annotate_df, data_ready), file= str_c(file_prefix, "_Minimum.csv", collapse = " "))
} else {
  data_ready[is.na(data_ready)] <- 0.0}


#change column names to from PD default to user defined
colnames(data_ready) <- sample_header[1:sample_number]


#save dataframe with ADH peptides for QC meteric
ADH_data_raw <-subset(cbind(annotate_df, data_ready), Accession %in% adh_list)
ADH_data_raw$holes <- rowSums(ADH_data_raw[(info_columns+1):ncol(ADH_data_raw)] ==0)
ADH_data_raw <- subset(ADH_data_raw, holes==0)
ADH_data_raw <- ADH_data_raw[ , -ncol(ADH_data_raw)]
ADH_data_raw$cv <- percentCV_gw(ADH_data_raw[(info_columns+1):ncol(ADH_data_raw)])
ADH_data_raw$sd <- apply(ADH_data_raw[(info_columns+1):(ncol(ADH_data_raw)-1)], 1, FUN = function(x) {sd(x)})
ADH_data_raw$av <- apply(ADH_data_raw[(info_columns+1):(ncol(ADH_data_raw)-1)], 1, FUN = function(x) {mean(x)})
ADH_data_raw$id <- seq(1, nrow(ADH_data_raw), by=1)
  
# log2 data for normalization
if (log_normalize){
  data_ready <- log(data_ready,2)
  data_ready[data_ready==-Inf] = 0}  # fix log2 of 0


#normalize on data with no missing values, create new data frame
data_ready$holes <- rowSums(data_ready == 0.0)
data_normalize <- subset(data_ready, holes==0)
data_normalize <- data_normalize[1:sample_number]
data_ready <- data_ready[1:sample_number]

#---------------------------------------------
# LR Normalized 
#--------------------------------------------
data_ready_lr <- lr_normalize(data_normalize,data_ready)
write.csv(cbind(annotate_df, data.frame(2^data_ready_lr)), file= str_c(file_prefix, "_LR_Normalized.csv", collapse = " "))

#---------------------------------------------
# SL Normalized 
#--------------------------------------------
# global scaling value, sample loading normalization
target <- mean(colSums(data_normalize))
norm_facs <- target / colSums(data_normalize)
data_ready_sl <- sweep(data_ready, 2, norm_facs, FUN = "*")
write.csv(cbind(annotate_df, data.frame(2^data_ready_sl)), file= str_c(file_prefix, "_SL_Normalized.csv", collapse = " "))

#---------------------------------------------
# TMM Normalized 
#--------------------------------------------
raw_tmm <- calcNormFactors(data_normalize, method = "TMM", sumTrim = 0.1)
data_ready_tmm <- sweep(data_ready, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale
write.csv(cbind(annotate_df, data.frame(2^data_ready_tmm)), file= str_c(file_prefix, "_TMM_Normalized.csv", collapse = " "))

#---------------------------------------------
# TMM & SL Normalized, step moved to after imputation of TMM
#--------------------------------------------
if (holes != "Impute"){
  sl_tmm <- calcNormFactors(data_ready_sl, method = "TMM", sumTrim = 0.1)
  data_ready_sl_tmm <- sweep(data_ready_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
  write.csv(cbind(annotate_df, data.frame(2^data_ready_sl_tmm)), file= str_c(file_prefix, "_SLTMM_Normalized.csv", collapse = " "))
}

#---------------------------------------------
# Specific Protein Normalized Data, ie carboxylases
#--------------------------------------------
# global scaling value, sample loading normalization
if (normalize_to_protein) {
  protein_norm_raw <- cbind(annotate_df, data_ready)
  protein_norm_raw <-subset(protein_norm_raw, Accessions %in% protein_norm_list)
  protein_norm_raw <- protein_norm_raw[(info_columns+1):ncol(protein_norm_raw)]
  protein_norm_raw$holes <- rowSums(protein_norm_raw == 0.0)
  protein_norm_raw <- subset(protein_norm_raw, holes==0)
  protein_norm_raw <- protein_norm_raw[1:sample_number]
  target <- mean(colSums((protein_norm_raw)))
  norm_facs <- target / colSums(protein_norm_raw)
  data_ready_protein_norm <- sweep(data_ready, 2, norm_facs, FUN = "*")
  write.csv(cbind(annotate_df, data.frame(2^data_ready_protein_norm)), file= str_c(file_prefix, "_Protein_Normalized.csv", collapse = " "))
  
}


#----------------------------------
# fill missing data, "unlog"
#----------------------------------------
if (holes == "Impute"){
  data_ready_fill <- hole_fill(data_ready)
  data_ready_sl <- hole_fill(data_ready_sl)
  data_ready_tmm <- hole_fill(data_ready_tmm)
  data_ready_lr <- hole_fill(data_ready_lr)
  if(normalize_to_protein) {
    data_ready_protein_norm_missing <- data_ready_protein_norm
    data_ready_protein_norm <- hole_fill(data_ready_protein_norm)
    }

  sl_tmm <- calcNormFactors(data_ready_sl, method = "TMM", sumTrim = 0.1)
  data_ready_sl_tmm <- sweep(data_ready_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
}

if (log_normalize){
  data_ready <- data.frame(2^(data_ready))
  data_ready_fill <- data.frame(2^(data_ready_fill))
  data_ready_sl <- data.frame(2^(data_ready_sl))
  data_ready_tmm <- data.frame(2^(data_ready_tmm))
  data_ready_sl_tmm <- data.frame(2^(data_ready_sl_tmm))
  data_ready_lr <- data.frame(2^(data_ready_lr))
  if(normalize_to_protein) {
    data_ready_protein_norm <- data.frame(2^(data_ready_protein_norm))
    }
}

# fix unlog of 0
data_ready[data_ready ==1 ] <- 0


write.csv(cbind(annotate_df, data_ready_fill), file= str_c(file_prefix, "_Raw_Impute.csv", collapse = " "))
write.csv(cbind(annotate_df, data_ready_sl), file= str_c(file_prefix, "_SL_Impute.csv", collapse = " "))
write.csv(cbind(annotate_df, data_ready_tmm), file= str_c(file_prefix, "_TMM_Impute.csv", collapse = " "))
write.csv(cbind(annotate_df, data_ready_sl_tmm), file= str_c(file_prefix, "_SLTMM_Impute.csv", collapse = " "))
write.csv(cbind(annotate_df, data_ready_lr), file= str_c(file_prefix, "_LR_Impute.csv", collapse = " "))


#-------------------------------------------
# Plots, boxplot, MDS, density, barplot
#------------------------------------------

Plot_All_gw(data_ready, "Raw Data")
Plot_All_gw(data_ready_lr, "Lin Reg")
Plot_All_gw(data_ready_fill, "Raw Impute")
Plot_All_gw(data_ready_sl, "SL Normalized")
Plot_All_gw(data_ready_tmm, "TMM Normalized")
Plot_All_gw(data_ready_sl_tmm, "SL TMM Normalized")
if (normalize_to_protein) {Plot_All_gw(data_ready, "Protein Normalized")}


#--recombine annotation and data
data_ready <- data.frame(annotate_df, data_ready)
data_ready_lr <- data.frame(annotate_df, data_ready_lr)
data_ready_fill <- data.frame(annotate_df, data_ready_fill)
data_ready_sl <- data.frame(annotate_df, data_ready_sl)  
data_ready_sl_tmm <- data.frame(annotate_df, data_ready_sl_tmm)  
data_ready_tmm <- data.frame(annotate_df, data_ready_tmm) 
if (normalize_to_protein) {data_ready_protein_norm <- data.frame(annotate_df, data_ready_protein_norm)}




#collapse psm to peptide-----------------------------------------------
if (psm_input){
  data_ready <- collapse_psm(data_ready)
  data_ready_fill <- collapse_psm(data_ready_fill)
  data_ready_sl <- collapse_psm(data_ready_sl)
  data_ready_sl_tmm <- collapse_psm(data_ready_sl_tmm)
  data_ready_tmm <- collapse_psm(data_ready_tmm)
  #data_list <- collapse_psm(data_ready)
  #ata_ready <- data_list[[1]]
  #annotate_df <- data_list[[2]]
}


# collapse peptide to protein
if (peptide_to_protein){
  data_ready <- collapse_peptide(data_ready)
  data_ready_lr <- collapse_peptide(data_ready_lr)
  data_ready_fill <- collapse_peptide(data_ready_fill)
  data_ready_sl <- collapse_peptide(data_ready_sl)
  data_ready_sl_tmm <- collapse_peptide(data_ready_sl_tmm)
  data_ready_tmm <- collapse_peptide(data_ready_tmm)
  if (normalize_to_protein) {data_ready_protein_norm <- collapse_peptide(data_ready_protein_norm)}
  # reassign column names, format change from peptide to protein
  info_columns <- ncol(data_ready) - sample_number
  info_headers <- colnames(data_ready[1:info_columns])
  final_sample_header <- c(info_headers, sample_header)
}


#-----------------------------------------------------------------------------------------
# stats
#-----------------------------------------------------------------------------------------

#PCA_gw(data_ready[(info_columns+1):ncol(data_ready)], "Raw Data")
PCA_gw(data_ready_fill[(info_columns+1):ncol(data_ready_fill)], "Raw Filled")
PCA_gw(data_ready_sl[(info_columns+1):ncol(data_ready_sl)], "SL Normalized")
PCA_gw(data_ready_tmm[(info_columns+1):ncol(data_ready_tmm)], "TMM Normalized")
PCA_gw(data_ready_sl_tmm[(info_columns+1):ncol(data_ready_sl_tmm)], "TMM SL Normalized")
if (normalize_to_protein) {PCA_gw(data_ready_protein_norm[(info_columns+1):ncol(data_ready_protein_norm)], 
                                  "Specific Protein Normalized")}



data_ready_fill_final <- stat_test_gw(data_ready_fill[1:info_columns], 
                                 data_ready_fill[(info_columns+1):ncol(data_ready_fill)],
                                 "Impute Data")
data_ready_sl_final <- stat_test_gw(data_ready_sl[1:info_columns], 
                                    data_ready_sl[(info_columns+1):ncol(data_ready_sl)],
                                    "SL Normalized")
data_ready_tmm_final <- stat_test_gw(data_ready_tmm[1:info_columns], 
                                     data_ready_tmm[(info_columns+1):ncol(data_ready_tmm)],
                                     "TMM Normalized")
data_ready_sl_tmm_final <- stat_test_gw(data_ready_sl_tmm[1:info_columns], 
                                        data_ready_sl_tmm[(info_columns+1):ncol(data_ready_sl_tmm)],
                                        "TMM SL Normalized")
if (normalize_to_protein) {data_ready_protein_norm_final <- stat_test_gw(data_ready_protein_norm[1:info_columns], 
                                       data_ready_protein_norm[(info_columns+1):ncol(data_ready_protein_norm)],
                                     "Specific Protein Normalized")}


#fix headers
#colnames(data_ready_final) <- final_sample_header
colnames(data_ready_fill_final) <- final_sample_header
colnames(data_ready_sl_final) <- final_sample_header
colnames(data_ready_tmm_final) <- final_sample_header
colnames(data_ready_sl_tmm_final) <- final_sample_header
if (normalize_to_protein) {colnames(data_ready_protein_norm_final) <- final_sample_header}

  #--csv output
#write.csv(data.frame(data_ready_final), file= str_c(file_prefix, "_final.csv", collapse = " "))
#write.csv(data.frame(data_ready_fill_final), file= str_c(file_prefix, "_final.csv", collapse = " "))
#write.csv(data.frame(data_ready_sl_final), file= str_c(file_prefix, "_sl_final.csv", collapse = " "))
#write.csv(data.frame(data_ready_tmm_final), file= str_c(file_prefix, "_tmm_final.csv", collapse = " "))
#write.csv(data.frame(data_ready_sl_tmm_final), file= str_c(file_prefix, "_sl_tmm_final.csv", collapse = " "))

# create plots for BirA, Carbox, Avidin, Bait, ADH
#BioID_normalize_gw(data_ready, "Raw")
BioID_normalize_gw(data_ready, "Raw")
BioID_normalize_gw(data_ready_fill, "Fill")
BioID_normalize_gw(data_ready_sl, "SL")
BioID_normalize_gw(data_ready_tmm, "TMM")
BioID_normalize_gw(data_ready_sl_tmm, "SL TMM")
if (normalize_to_protein) {BioID_normalize_gw(data_ready_protein_norm, "Protein Norm")}


cv_stats()


# final excel tables with tabs for each comparison, subset by pval, fc
#Final_Excel_gw(data_ready_final, "_raw_final.xlsx")
Final_Excel_gw(data_ready_fill_final, "_fill_final.xlsx")
Final_Excel_gw(data_ready_sl_final, "_sl_final.xlsx")
Final_Excel_gw(data_ready_tmm_final, "_tmm_final.xlsx")
Final_Excel_gw(data_ready_sl_tmm_final, "_sl_tmm_final.xlsx")
if (normalize_to_protein) {Final_Excel_gw(data_ready_protein_norm_final, "_specific_protein_final.xlsx")}












