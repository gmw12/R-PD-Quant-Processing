# !diagnostics off
options(digits=3)

#format data from nested protein/peptide output to peptide only
if (protein_peptide_input){forward_data <- protein_to_peptide(forward_data)}


# clean column headers and save 
total_columns <- ncol(forward_data)
info_columns <- total_columns - sample_number
final_sample_headers <- clean_headers(colnames(forward_data[1:info_columns]))
colnames(forward_data) <- final_sample_headers

data_raw <- forward_data

# create column to flag and delete rows with no data or only one data value
data_raw$na_count <- apply(data_raw[(info_columns+1):ncol(data_raw)], 1, function(x) sum(!is.na(x)))
data_raw <- subset(data_raw, data_raw$na_count > 1)
data_raw <- data_raw[1:total_columns]
Simple_Excel(data_raw, str_c(file_prefix3, "_raw_input.xlsx", collapse = " "))

# save the annotation columns for later and remove from data frame
annotate_df <- data_raw[1:info_columns]
data_raw <- data_raw[(info_columns+1):total_columns]

#reorder data if needed, will use PD_order from sample list
data_raw <- order_columns(data_raw)

# create histogram of all measured values in data set, compare against hole fill confidence value
histogram_gw(data_raw,"Total_Set_Intensity","Log2 Intensity Distribution")

#save dataframe with hole locations
#hole_df <- data_raw
#hole_df[hole_df > 0] <- 1
#hole_df[is.na(hole_df)] <- 0


#data_raw = no holes filled, data_fil = fill holes according to user choice
data_fill <- data_raw

# save raw and process plots
data_raw[is.na(data_raw)] <- 0.0
Simple_Excel(cbind(annotate_df, data_raw), str_c(file_prefix3, "_Peptide_Raw.xlsx", collapse = " "))
Plot_All_gw(data_raw, "Raw")
data_raw <- cbind(annotate_df, data_raw)
if (peptide_to_protein){data_raw <- collapse_peptide(data_raw)}
# if peptide to protein need to reassign info columns and headers
if (peptide_to_protein){
  info_columns_protein <- ncol(data_raw) - sample_number
  info_headers_protein <- colnames(data_raw[1:info_columns_protein])
  final_sample_header_protein <- c(info_headers_protein, sample_info$Header1)
  final_sample_header_protein_norm <- c(info_headers_protein, sample_info$Header2)
}else{
  info_columns_protein <- info_columns
  info_headers_protein <- colnames(data_raw[1:info_columns_protein])
  final_sample_header_protein <- c(info_headers_protein, sample_info$Header1)
  final_sample_header_protein_norm <- c(info_headers_protein, sample_info$Header2)
}


#create dataframe for CV compairson box plots,  needs protein accession list to start
total_cv <- data.frame(data_raw$Accession)

# Impute strategy
if (holes == "Floor") {
  data_fill[is.na(data_fill)] <- area_floor
  Simple_Excel(cbind(annotate_df, data_fill), str_c(file_prefix3, "_Floor.csv", collapse = " "))
} else if (holes == "Average") {
  data_fill <- hole_average(data_fill)
  Simple_Excel(cbind(annotate_df, data_fill), str_c(file_prefix3, "_Average.csv", collapse = " "))
} else if (holes == "Minimium") {
  data_fill <- hole_minimum(data_fill)
  Simple_Excel(cbind(annotate_df, data_fill), str_c(file_prefix3, "_Minimum.csv", collapse = " "))
} else {
  data_fill[is.na(data_fill)] <- 0.0}


#save dataframe with ADH peptides for QC meteric
if (adh_spike){
  ADH_data_fill <-subset(cbind(annotate_df, data_fill), Accession %in% adh_list)
  ADH_data_fill$holes <- rowSums(ADH_data_fill[(info_columns+1):ncol(ADH_data_fill)] ==0)
  ADH_data_fill <- subset(ADH_data_fill, holes==0)
  ADH_data_fill <- ADH_data_fill[ , -ncol(ADH_data_fill)]
  ADH_data_fill$cv <- percentCV_gw(ADH_data_fill[(info_columns+1):ncol(ADH_data_fill)])
  ADH_data_fill$sd <- apply(ADH_data_fill[(info_columns+1):(ncol(ADH_data_fill)-1)], 1, FUN = function(x) {sd(x)})
  ADH_data_fill$av <- apply(ADH_data_fill[(info_columns+1):(ncol(ADH_data_fill)-1)], 1, FUN = function(x) {mean(x)})
  ADH_data_fill$id <- seq(1, nrow(ADH_data_fill), by=1)
}

#place holder that would allow a subset of data to be used for normalization ...values with no missing data...
data_normalize <- data_fill


#---------------------------------------------------------------------------------------------
#Normalize Section, normalize -> plots -> peptide to protein -> stats/PCA
#---------------------------------------------------------------------------------------------

# SL - global scaling value, sample loading normalization I * (avgSumAvgI/AvgI)
data_sl <- sl_normalize(data_normalize, data_fill, "SL_Norm")
data_sl_final <- dostats(data_sl, "SL_Norm")

# TI - global scaling value, sample loading normalization I * (medianSumI/SumI)
data_ti <- ti_normalize(data_normalize, data_fill, "TI_Norm")
data_ti_final <- dostats(data_ti, "TI_Norm")

# MI - global scaling value, sample loading normalization  I * (meanMedianI/MedianI)
data_mi <- mi_normalize(data_normalize, data_fill, "MI_Norm")
data_mi_final <- dostats(data_mi, "MI_Norm")

# AI - global scaling value, sample loading normalization  I * (meanMeanI/MeanI)
data_ai <- ai_normalize(data_normalize, data_fill, "AI_Norm")
data_ai_final <- dostats(data_ai, "AI_Norm")

# TMM Normalized - trimmed mean normalization, coded for middle 80 (-10 top and bottom)
data_tmm <- tmm_normalize(data_normalize, data_fill, "TMM_Norm")
data_tmm_final <- dostats(data_tmm, "TMM_Norm")

# SLTMM Normalized - SL first then trimmed mean normalization, coded for middle 80 (-10 top and bottom)
data_sltmm <- sltmm_normalize(data_normalize, data_fill, "SLTMM_Norm")
data_sltmm_final <- dostats(data_sltmm, "SLTMM_Norm")

# LR Normalized  - simple linear regression normalization
data_lr <- lr_normalize(data_fill, "LinReg_Norm")
data_lr_final <- dostats(data_lr, "LinReg_Norm")

# VSN Normalized
data_vsn <- vsn_normalize(data_fill, "VSN_Norm")
data_vsn_final <- dostats(data_vsn, "VSN_Norm")

# Quantile Normalized
data_quantile <- quantile_normalize(data_fill, "Quantile_Norm")
data_quantile_final <- dostats(data_quantile, "Quantile_Norm")

# LOESS Normalized
data_loess <- loess_normalize(data_fill, "LOESS_Norm")
data_loess_final <- dostats(data_loess, "LOESS_Norm")

# Local SL - local scaling value, sample loading normalization
data_localsl <- sl_local_normalize(data_normalize, data_fill, "LocalSL_Norm")
data_localsl_final <- dostats(data_localsl, "LocalSL_Norm")

# ProteinNorm - Specific Protein Normalized Data, ie carboxylases
if (normalize_to_protein) {data_protein_norm <- protein_normalize(data_fill, "Protein_Norm")
  data_protein_norm_final <- dostats(data_protein_norm, "Protein_Norm")}

# Raw - impute raw or fill raw plots
if (holes=="Impute") {
  data_raw_impute <- hole_fill(data_fill)
  Simple_Excel(cbind(annotate_df, data_raw_impute), str_c(file_prefix3, "_Peptide_Raw_Impute.xlsx", collapse = " "))
  Plot_All_gw(data_raw_impute, "Raw Impute")
  data_raw_impute <- cbind(annotate_df, data_raw_impute)
  if (peptide_to_protein){data_raw_impute <- collapse_peptide(data_raw_impute)}
}else{
  Simple_Excel(cbind(annotate_df, data_fill), str_c(file_prefix3, "_Peptide_FilledRaw.xlsx", collapse = " "))
  Plot_All_gw(data_fill, "FilledRaw")
  data_fill <- cbind(annotate_df, data_fill)
  if (peptide_to_protein){data_fill <- collapse_peptide(data_fill)}
}


cv_stats(summary_cv, total_cv)


data_raw_final <- dostats(data_raw, "Raw")









