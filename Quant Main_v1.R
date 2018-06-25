#global options, numbers, sig digits
options(scipen = 999)
options(digits=3)

if (psm_input) {
  data_raw <- psm_decoy(forward_data, decoy_data)
  }else{
    data_raw <- forward_data 
  }  

#set NA to 1 for protein, 0.01 for TMT
data_raw[is.na(data_raw)] <- area_floor

#----- edit column headers
col_headers <- colnames(data_raw) 
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Combined", "Confidence")
col_headers <- str_replace(col_headers, "Annotated Sequence", "Annotated_Sequence")
col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
col_headers <- str_replace(col_headers,"\\[", "")
col_headers <- str_replace(col_headers,"\\]", "")
#col_headers <- str_replace(col_headers, "Abundance: ", "")
#col_headers <- str_replace(col_headers, "Sample, ", "")
#col_headers <- str_replace(col_headers, "F13: ", "")
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
data_ready <- subset(data_ready, total_row > sample_number) 
data_ready <- data_ready[1:total_columns]

# save the annotation columns for later and remove from data frame
annotate_df <- data_ready[1:info_columns]
data_ready <- data_ready[(info_columns+1):total_columns]
row.names(data_ready) <- annotate_df$`Accession`

#remove unused data 
rm(data_raw, forward_data, decoy_data, total_row)

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
# Carbox Normalized Data
#--------------------------------------------
# global scaling value, sample loading normalization
#carbox_list <- c("Q05920", "Q91ZA3")
carbox_list <- c("P11498", "P05166", "Q96RQ3", "Q13085")
carbox_raw <-subset(data_ready, rownames(data_ready) %in% carbox_list)
target <- mean(colSums((carbox_raw)))
norm_facs <- target / colSums(carbox_raw)
data_ready_carbox <- sweep(data_ready, 2, norm_facs, FUN = "*")

#-------------------------------------------
# Plots
#------------------------------------------

# Raw, NOT normalized ----------------------
boxplot_gw(data_ready, "Raw Data")
plotMDS_gw(data_ready,"Raw Data_Multidimension Scaling")
data_ready_bar <- colSums(data_ready)
barplot_gw(data_ready_bar, "Raw Data")
plotDensities_gw(data_ready, "Raw data")

#--Total, Normalize Data---------------------------------
boxplot_gw(data_ready_sl, "SL Normalized")
plotMDS_gw(data_ready_sl,"SL Normalized Data_Multidimension Scaling")
data_ready_bar <- colSums(data_ready_sl)
barplot_gw(data_ready_bar, "SL Normalized")
plotDensities_gw(data_ready_sl, "SL Normalized")

#----TMM from Raw--------------------------
boxplot_gw(data_ready_tmm, "TMM Normalized")
plotMDS_gw(data_ready_tmm,"TMM Normalized_Multidimension Scaling")
data_ready_bar <- colSums(data_ready_tmm)
barplot_gw(data_ready_bar, "TMM Data")
plotDensities_gw(data_ready_tmm, "TMM Normalized")

#---TMM, SL Normalized Data-------------------------------------------
boxplot_gw(data_ready_sl_tmm, "TMM SL Normalized")
plotMDS_gw(data_ready_sl_tmm,"TMM SL Normalized Data_Multidimension Scaling")
data_ready_bar <- colSums(data_ready_sl_tmm)
barplot_gw(data_ready_bar, "TMM SL Normalized")
plotDensities_gw(data_ready_sl_tmm, "TMM SL Normalized")

#--Carbox, Normalized Data---------------------------------
boxplot_gw(data_ready_carbox, "Carbox, Normalized")
plotMDS_gw(data_ready_carbox,"Carbox, Normalized Data_Multidimension Scaling")
data_ready_bar <- colSums(data_ready_carbox)
barplot_gw(data_ready_bar, "Carbox, Normalized")
plotDensities_gw(data_ready_carbox, "Carbox, Normalized")


#--recombine annotation and data
data_ready <- data.frame(annotate_df, data_ready)
data_ready_sl <- data.frame(annotate_df, data_ready_sl)  
data_ready_sl_tmm <- data.frame(annotate_df, data_ready_sl_tmm)  
data_ready_tmm <- data.frame(annotate_df, data_ready_tmm) 
data_ready_carbox <- data.frame(annotate_df, data_ready_carbox) 
colnames(data_ready) <- col_headers
colnames(data_ready_sl) <- col_headers
colnames(data_ready_sl_tmm) <- col_headers
colnames(data_ready_tmm) <- col_headers
colnames(data_ready_carbox) <- col_headers

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
PCA_gw(data_ready_carbox[(info_columns+1):ncol(data_ready_carbox)], "Carbox Normalized")


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
data_ready_carbox_final <- stat_test_gw(data_ready_carbox[1:info_columns], 
                                     data_ready_carbox[(info_columns+1):ncol(data_ready_carbox)],
                                     "Carbox Normalized")


#fix headers
colnames(data_ready_final) <- final_sample_header
colnames(data_ready_sl_final) <- final_sample_header
colnames(data_ready_tmm_final) <- final_sample_header
colnames(data_ready_sl_tmm_final) <- final_sample_header
colnames(data_ready_carbox_final) <- final_sample_header

#--csv for large peptide output
write.csv(data.frame(data_ready_final), file= str_c(file_prefix, "_final.csv", collapse = " "))
write.csv(data.frame(data_ready_sl_final), file= str_c(file_prefix, "_sl_final.csv", collapse = " "))
write.csv(data.frame(data_ready_tmm_final), file= str_c(file_prefix, "_tmm_final.csv", collapse = " "))
write.csv(data.frame(data_ready_sl_tmm_final), file= str_c(file_prefix, "_sl_tmm_final.csv", collapse = " "))
write.csv(data.frame(data_ready_carbox_final), file= str_c(file_prefix, "_carbox_final.csv", collapse = " "))

# csv for all peptides if filtering for phos for report
if (phos_peptide_only){
  write.csv(data.frame(data_peptide), file= str_c(file_prefix, "_peptide.csv", collapse = " "))
  write.csv(data.frame(data_peptide_sl), file= str_c(file_prefix, "_sl_peptide.csv", collapse = " "))
  write.csv(data.frame(data_peptide_tmm), file= str_c(file_prefix, "_tmm_peptide.csv", collapse = " "))
  write.csv(data.frame(data_peptide_sl_tmm), file= str_c(file_prefix, "_sl_tmm_peptide.csv", collapse = " ")) 
}

BioID_normalize_gw(data_ready, "Raw")
BioID_normalize_gw(data_ready_sl, "SL")
BioID_normalize_gw(data_ready_tmm, "TMM")
BioID_normalize_gw(data_ready_sl_tmm, "SL TMM")
BioID_normalize_gw(data_ready_carbox, "Carbox")





















#--- directly to excel for protein projects
wb = createWorkbook()
sheet = createSheet(wb, "Sheet 1")
addDataFrame(data.frame(data_ready_sl_tmm_list[1]), sheet=sheet, startColumn=1, row.names=FALSE)
sheet = createSheet(wb, "Sheet 2")
addDataFrame(data.frame(data_ready_sl_tmm_list[2]), sheet=sheet, startColumn=1, row.names=FALSE)
saveWorkbook(wb, "4227_TMT_MS3_MRM_033018_sl_tmm.xlsx")

wb = createWorkbook()
sheet = createSheet(wb, "Sheet 1")
addDataFrame(data.frame(data_ready_sl_list[1]), sheet=sheet, startColumn=1, row.names=FALSE)
sheet = createSheet(wb, "Sheet 2")
addDataFrame(data.frame(data_ready_sl_list[2]), sheet=sheet, startColumn=1, row.names=FALSE)
saveWorkbook(wb, "4227_TMT_MS3_MRM_032318_sl.xlsx")
