#global options, numbers, sig digits
options(scipen = 999)
options(digits=3)

#import data from csv file
if (psm_input) {
  data_raw <- psm_decoy(forward_data, decoy_data)
  }else{
    data_raw <- forward_data 
  }  

#replaces NA, no data, with assigned value
data_raw[is.na(data_raw)] <- area_floor

#----- edit column headers
col_headers <- colnames(data_raw) 
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Combined", "Confidence")
col_headers <- str_replace(col_headers, "Annotated Sequence", "Annotated_Sequence")
col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
col_headers <- str_replace(col_headers,"\\[", "")
col_headers <- str_replace(col_headers,"\\]", "")
colnames(data_raw) <- col_headers
total_columns <- ncol(data_raw)
info_columns <- total_columns - sample_number
info_headers <- colnames(data_raw[1:info_columns])
final_sample_header <- c(info_headers, sample_header)


# filter list, only "High", delete proteins with no ID in any sample, only if PSM is false
if (!psm_input){
  data_ready <-data_raw[grepl("High", data_raw$Confidence, ignore.case=TRUE),] 
}else{
  data_ready <- data_raw
}

# create column to sum abundances, flag and delete rows with no data
total_row <- rowSums(data_ready[(info_columns+1):total_columns])
total_row <- data_frame(total_row)
data_ready <- cbind(data_ready, total_row)
data_ready <- subset(data_ready, total_row > (sample_number*area_floor)) 
data_ready <- data_ready[1:total_columns]

# save the annotation columns for later and remove from data frame
annotate_df <- data_ready[1:info_columns]
data_ready <- data_ready[(info_columns+1):total_columns]
row.names(data_ready) <- annotate_df$`Accession`

#remove unused data 
try(rm(data_raw, forward_data, decoy_data, total_row), silent = TRUE)

#arrange columns if psm data
if (psm_input){
  data_ready <- order_columns(data_ready)
}

colnames(data_ready) <- sample_header[1:sample_number]

#---------------------------------------------
# SL Normalized 
#--------------------------------------------
# global scaling value, sample loading normalization
target <- mean(colSums(data_ready))
norm_facs <- target / colSums(data_ready)
data_ready_sl <- sweep(data_ready, 2, norm_facs, FUN = "*")

#---------------------------------------------
# TMM Normalized 
#--------------------------------------------
raw_tmm <- calcNormFactors(data_ready, method = "TMM", sumTrim = 0.1)
data_ready_tmm <- sweep(data_ready, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale

#---------------------------------------------
# TMM & SL Normalized
#--------------------------------------------
# see exactly what TMM does with SL data
sl_tmm <- calcNormFactors(data_ready_sl, method = "TMM", sumTrim = 0.1)
data_ready_sl_tmm <- sweep(data_ready_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale

#---------------------------------------------
# Specific Protein Normalized Data, ie carboxylases
#--------------------------------------------
# global scaling value, sample loading normalization
if (normalize_to_protein == TRUE) {
  protein_norm_raw <-subset(data_ready, rownames(data_ready) %in% protein_norm_list)
  target <- mean(colSums((protein_norm_raw)))
  norm_facs <- target / colSums(protein_norm_raw)
  data_ready_protein_norm <- sweep(data_ready, 2, norm_facs, FUN = "*")
}

#-------------------------------------------
# Plots, boxplot, MDS, density, barplot
#------------------------------------------

Plot_All_gw(data_ready, "Raw Data")
Plot_All_gw(data_ready_sl, "SL Normalized")
Plot_All_gw(data_ready_tmm, "TMM Normalized")
Plot_All_gw(data_ready_sl_tmm, "SL TMM Normalized")
if (normalize_to_protein == TRUE) {Plot_All_gw(data_ready, "Protein Normalized")}


#--recombine annotation and data
data_ready <- data.frame(annotate_df, data_ready)
data_ready_sl <- data.frame(annotate_df, data_ready_sl)  
data_ready_sl_tmm <- data.frame(annotate_df, data_ready_sl_tmm)  
data_ready_tmm <- data.frame(annotate_df, data_ready_tmm) 
colnames(data_ready) <- col_headers
colnames(data_ready_sl) <- col_headers
colnames(data_ready_sl_tmm) <- col_headers
colnames(data_ready_tmm) <- col_headers
if (normalize_to_protein == TRUE) {
  data_ready_protein_norm <- data.frame(annotate_df, data_ready_protein_norm)   
  colnames(data_ready_protein_norm) <- col_headers}

#collapse psm to peptide-----------------------------------------------
if (psm_input){
  data_ready <- collapse_psm(data_ready)
  data_ready_sl <- collapse_psm(data_ready_sl)
  data_ready_sl_tmm <- collapse_psm(data_ready_sl_tmm)
  data_ready_tmm <- collapse_psm(data_ready_tmm)
  
  #data_list <- collapse_psm(data_ready)
  #ata_ready <- data_list[[1]]
  #annotate_df <- data_list[[2]]
}


if (phos_peptide_only){
  data_peptide <- data_ready
  data_peptide_sl <- data_ready_sl
  data_peptide_sl_tmm <- data_ready_sl_tmm
  data_peptide_tmm <- data_ready_tmm
  
  data_ready <- data_peptide[grepl("Phospho", data_peptide$Modifications, ignore.case=TRUE),]
  data_ready_sl <- data_peptide_sl[grepl("Phospho", data_peptide_sl$Modifications, ignore.case=TRUE),]
  data_ready_sl_tmm <- data_peptide_sl_tmm[grepl("Phospho", data_peptide_sl_tmm$Modifications, ignore.case=TRUE),]
  data_ready_tmm <- data_peptide_tmm[grepl("Phospho", data_peptide_tmm$Modifications, ignore.case=TRUE),]
}


#-----------------------------------------------------------------------------------------
# stats
#-----------------------------------------------------------------------------------------

PCA_gw(data_ready[(info_columns+1):ncol(data_ready)], "Raw Data")
PCA_gw(data_ready_sl[(info_columns+1):ncol(data_ready_sl)], "SL Normalized")
PCA_gw(data_ready_tmm[(info_columns+1):ncol(data_ready_tmm)], "TMM Normalized")
PCA_gw(data_ready_sl_tmm[(info_columns+1):ncol(data_ready_sl_tmm)], "TMM SL Normalized")



data_ready_final <- stat_test_gw(data_ready[1:info_columns], 
                                 data_ready[(info_columns+1):ncol(data_ready)],
                                 "Raw Data")
data_ready_sl_final <- stat_test_gw(data_ready_sl[1:info_columns], 
                                    data_ready_sl[(info_columns+1):ncol(data_ready_sl)],
                                    "SL Normalized")
data_ready_tmm_final <- stat_test_gw(data_ready_tmm[1:info_columns], 
                                     data_ready_tmm[(info_columns+1):ncol(data_ready_tmm)],
                                     "TMM Normalized")
data_ready_sl_tmm_final <- stat_test_gw(data_ready_sl_tmm[1:info_columns], 
                                        data_ready_sl_tmm[(info_columns+1):ncol(data_ready_sl_tmm)],
                                        "TMM SL Normalized")
if (normalize_to_protein == TRUE) {
  PCA_gw(data_ready_protein_norm[(info_columns+1):ncol(data_ready_protein_norm)], "Specific Protein Normalized")
  data_ready_protein_norm_final <- stat_test_gw(data_ready_protein_norm[1:info_columns], 
                                       data_ready_protein_norm[(info_columns+1):ncol(data_ready_protein_norm)],
                                     "Specific Protein Normalized")}


#fix headers
colnames(data_ready_final) <- final_sample_header
colnames(data_ready_sl_final) <- final_sample_header
colnames(data_ready_tmm_final) <- final_sample_header
colnames(data_ready_sl_tmm_final) <- final_sample_header

#--csv output
write.csv(data.frame(data_ready_final), file= str_c(file_prefix, "_final.csv", collapse = " "))
write.csv(data.frame(data_ready_sl_final), file= str_c(file_prefix, "_sl_final.csv", collapse = " "))
write.csv(data.frame(data_ready_tmm_final), file= str_c(file_prefix, "_tmm_final.csv", collapse = " "))
write.csv(data.frame(data_ready_sl_tmm_final), file= str_c(file_prefix, "_sl_tmm_final.csv", collapse = " "))

# create plots for BirA, Carbox, Avidin, Bait, ADH
BioID_normalize_gw(data_ready, "Raw")
BioID_normalize_gw(data_ready_sl, "SL")
BioID_normalize_gw(data_ready_tmm, "TMM")
BioID_normalize_gw(data_ready_sl_tmm, "SL TMM")


if (normalize_to_protein == TRUE) {
  colnames(data_ready_protein_norm_final) <- final_sample_header
  write.csv(data.frame(data_ready_protein_norm_final), file= str_c(file_prefix, "_protein_norm_final.csv", collapse = " "))
  BioID_normalize_gw(data_ready_protein_norm, "Protein Norm")}



# create summary table for CV's for groups under different normalize conditions
cv_start <- info_columns+sample_number+sample_number+1
raw_avg_cv <- colMeans(data_ready_final[cv_start:(cv_start+group_number-1)])
sl_avg_cv <- colMeans(data_ready_sl_final[cv_start:(cv_start+group_number-1)])
tmm_avg_cv <- colMeans(data_ready_tmm_final[cv_start:(cv_start+group_number-1)])
sl_tmm_avg_cv <- colMeans(data_ready_sl_tmm_final[cv_start:(cv_start+group_number-1)])
if (normalize_to_protein == TRUE) {
  protein_avg_cv <- colMeans(data_ready_protein_norm_final[cv_start:(cv_start+group_number-1)])
  final_cv_avg <- rbind(raw_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv, protein_avg_cv )
}else{
  final_cv_avg <- rbind(raw_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv)
}

png(filename=str_c(output_dir, "Final_CV_Avg.png"), width = 888, height = 571)
p<-tableGrob(final_cv_avg)
grid.arrange(p)
dev.off()






# final excel tables with tabs for each comparison, subset by pval, fc
Final_Excel_gw(data_ready_final, "_raw_final.xlsx")
Final_Excel_gw(data_ready_sl_final, "_sl_final.xlsx")
Final_Excel_gw(data_ready_tmm_final, "_tmm_final.xlsx")
Final_Excel_gw(data_ready_sl_tmm_final, "_sl_tmm_final.xlsx")
if (normalize_to_protein == TRUE) {Final_Excel_gw(data_ready_protein_norm_final, "_specific_protein_final.xlsx")}





