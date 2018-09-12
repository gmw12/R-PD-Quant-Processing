#global options, numbers, sig digits
options(scipen = 999)
options(digits=3)


if (psm_input) {
  data_ready <- psm_decoy(forward_data, decoy_data)
  }else{
    data_ready <- forward_data 
  }  

#replaces NA, no data, with assigned value
if (holes == "Floor") {
  data_ready[is.na(data_ready)] <- area_floor
} else if (holes == "Average") {
  data_ready <- hole_average(data_ready)
} else if (holes == "Minimium") {
  data_ready <- hole_minimum(data_ready)
} else {
  data_ready[is.na(data_ready)] <- 0.0}


# create column to sum abundances, flag and delete rows with no data
data_ready$test <- rowSums(data_ready[(info_columns+1):total_columns])
data_ready <- subset(data_ready, data_ready$test > (sample_number * area_floor))
data_ready <- data_ready[1:total_columns]


# save the annotation columns for later and remove from data frame
annotate_df <- data_ready[1:info_columns]
data_ready <- data_ready[(info_columns+1):total_columns]
#row.names(data_ready) <- annotate_df$`Accessions` delete later if no errors 

#arrange columns if psm data
if (psm_input){
  data_ready <- order_columns(data_ready)
}

#change column names to from PD default to user defined
colnames(data_ready) <- sample_header[1:sample_number]

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
# SL Normalized 
#--------------------------------------------
# global scaling value, sample loading normalization
target <- mean(colSums(data_normalize))
norm_facs <- target / colSums(data_normalize)
data_ready_sl <- sweep(data_ready, 2, norm_facs, FUN = "*")

#---------------------------------------------
# TMM Normalized 
#--------------------------------------------
raw_tmm <- calcNormFactors(data_normalize, method = "TMM", sumTrim = 0.1)
data_ready_tmm <- sweep(data_ready, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale

#---------------------------------------------
# TMM & SL Normalized, step moved to after imputation of TMM
#--------------------------------------------
if (holes != "Impute"){
  sl_tmm <- calcNormFactors(data_ready_sl, method = "TMM", sumTrim = 0.1)
  data_ready_sl_tmm <- sweep(data_ready_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
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
}


#----------------------------------
# fill missing data, "unlog"
#----------------------------------------
if (holes == "Impute"){
  data_ready_fill <- hole_fill(data_ready)
  data_ready_sl <- hole_fill(data_ready_sl)
  data_ready_tmm <- hole_fill(data_ready_tmm)
  if(normalize_to_protein) {data_ready_protein_norm <- hole_fill(data_ready_protein_norm)}

  sl_tmm <- calcNormFactors(data_ready_sl, method = "TMM", sumTrim = 0.1)
  data_ready_sl_tmm <- sweep(data_ready_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
}

if (log_normalize){
  data_ready <- data.frame(2^(data_ready))
  data_ready_fill <- data.frame(2^(data_ready_fill))
  data_ready_sl <- data.frame(2^(data_ready_sl))
  data_ready_tmm <- data.frame(2^(data_ready_tmm))
  data_ready_sl_tmm <- data.frame(2^(data_ready_sl_tmm))
  if(normalize_to_protein) {data_ready_protein_norm <- data.frame(2^(data_ready_protein_norm))}
}

# fix unlog of 0
data_ready[data_ready ==1 ] <- 0

#-------------------------------------------
# Plots, boxplot, MDS, density, barplot
#------------------------------------------

Plot_All_gw(data_ready, "Raw Data")
Plot_All_gw(data_ready_fill, "Raw Impute")
Plot_All_gw(data_ready_sl, "SL Normalized")
Plot_All_gw(data_ready_tmm, "TMM Normalized")
Plot_All_gw(data_ready_sl_tmm, "SL TMM Normalized")
if (normalize_to_protein) {Plot_All_gw(data_ready, "Protein Normalized")}


#--recombine annotation and data
data_ready <- data.frame(annotate_df, data_ready)
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
BioID_normalize_gw(data_ready_fill, "Fill")
BioID_normalize_gw(data_ready_sl, "SL")
BioID_normalize_gw(data_ready_tmm, "TMM")
BioID_normalize_gw(data_ready_sl_tmm, "SL TMM")
if (normalize_to_protein) {BioID_normalize_gw(data_ready_protein_norm, "Protein Norm")}



# create summary table for CV's for groups under different normalize conditions
cv_start <- info_columns+sample_number+sample_number+1
#raw_avg_cv <- colMeans(data_ready_final[cv_start:(cv_start+group_number-1)])
fill_avg_cv <- colMeans(data_ready_fill_final[cv_start:(cv_start+group_number-1)])
sl_avg_cv <- colMeans(data_ready_sl_final[cv_start:(cv_start+group_number-1)])
tmm_avg_cv <- colMeans(data_ready_tmm_final[cv_start:(cv_start+group_number-1)])
sl_tmm_avg_cv <- colMeans(data_ready_sl_tmm_final[cv_start:(cv_start+group_number-1)])
if (normalize_to_protein) {
  protein_avg_cv <- colMeans(data_ready_protein_norm_final[cv_start:(cv_start+group_number-1)])
  final_cv_avg <- rbind(fill_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv, protein_avg_cv )
}else{
  final_cv_avg <- rbind(fill_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv)
}

png(filename=str_c(output_dir, "Final_CV_Avg.png"), width = 888, height = 571)
p<-tableGrob(final_cv_avg)
grid.arrange(p)
dev.off()


# final excel tables with tabs for each comparison, subset by pval, fc
#Final_Excel_gw(data_ready_final, "_raw_final.xlsx")
Final_Excel_gw(data_ready_fill_final, "_fill_final.xlsx")
Final_Excel_gw(data_ready_sl_final, "_sl_final.xlsx")
Final_Excel_gw(data_ready_tmm_final, "_tmm_final.xlsx")
Final_Excel_gw(data_ready_sl_tmm_final, "_sl_tmm_final.xlsx")
if (normalize_to_protein) {Final_Excel_gw(data_ready_protein_norm_final, "_specific_protein_final.xlsx")}












