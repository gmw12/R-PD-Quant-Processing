# !diagnostics off
start_time <- Sys.time()
options(digits=3)

#format data from nested protein/peptide output to peptide only
if (protein_peptide_input && !tmt_spqc_sets){data_raw <- protein_to_peptide(forward_data)}

if (psm_peptide_fdr == TRUE){data_raw <- peptide_psm_set_fdr()}

if (tmt_spqc_sets == TRUE){data_raw <- tmt_spqc_normalize()}

if (!protein_peptide_input && !peptide_to_protein){data_raw <- forward_data}

# data stat info
total_columns <- ncol(data_raw)
info_columns <- total_columns - sample_number

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

#clean column headers
sample_headers <- clean_headers(colnames(annotate_df))
colnames(annotate_df) <- sample_headers[1:info_columns]
colnames(data_raw) <- sample_info$Header1


#save dataframe with ADH peptides for QC meteric
if (adh_spike){
  ADH_data_fill <-subset(cbind(annotate_df, data_raw), Accession %in% adh_list)
  ADH_data_fill$missings <- rowSums(ADH_data_fill[(info_columns+1):ncol(ADH_data_fill)] ==0)
  ADH_data_fill <- subset(ADH_data_fill, missings==0)
  ADH_data_fill <- ADH_data_fill[ , -ncol(ADH_data_fill)]
  ADH_data_fill$cv <- percentCV_gw(ADH_data_fill[(info_columns+1):ncol(ADH_data_fill)])
  ADH_data_fill$sd <- apply(ADH_data_fill[(info_columns+1):(ncol(ADH_data_fill)-1)], 1, FUN = function(x) {sd(x)})
  ADH_data_fill$av <- apply(ADH_data_fill[(info_columns+1):(ncol(ADH_data_fill)-1)], 1, FUN = function(x) {mean(x)})
  ADH_data_fill$id <- seq(1, nrow(ADH_data_fill), by=1)
  ADH_data_fill$color <- color_list[1:nrow(ADH_data_fill)]
  barplot_adh(ADH_data_fill$cv[1:12], ADH_data_fill$Sequence[1:12],"ADH Peptide CV's", file_prefix1)
}

# create histogram of all measured values in data set, compare against missing fill confidence value
intensity_cutoff <- histogram_gw(data_raw,"Total_Set_Intensity","Log2 Intensity Distribution")

#save dataframe with missing locations
#missing_df <- data_raw
#missing_df[missing_df > 0] <- 1
#missing_df[is.na(missing_df)] <- 0

#data_raw = no missings filled, data_fil = fill missings according to user choice
data_to_normalize <- data_raw
data_to_normalize[is.na(data_to_normalize)] <- 0.0

# save raw and process plots
# if collapse will set new column strucdture
data_raw[is.na(data_raw)] <- 0.0
Simple_Excel(cbind(annotate_df, data_raw), str_c(file_prefix3, "_Peptide_Raw.xlsx", collapse = " "))
Plot_All_gw(data_raw, "Raw")
data_raw <- cbind(annotate_df, data_raw)
if (peptide_to_protein){data_raw <- collapse_peptide(data_raw)}


# if peptide to protein need to reassign info columns and headers
if (peptide_to_protein){
  info_columns_final <- ncol(data_raw) - sample_number
  info_headers_final <- colnames(data_raw[1:info_columns_final])
  sample_header_final <- c(info_headers_final, sample_info$Header1)
  sample_header_final_norm <- c(info_headers_final, sample_info$Header2)
}else{
  info_columns_final <- info_columns
  info_headers_final <- colnames(data_raw[1:info_columns_final])
  sample_header_final <- c(info_headers_final, sample_info$Header1)
  sample_header_final_norm <- c(info_headers_final, sample_info$Header2)
}

#create dataframe for CV compairson box plots,  needs protein accession list to start
total_cv <- data.frame(data_raw$Accession)

#---------------------------------------------------------------------------------------------
#Normalize Section, normalize -> plots -> peptide to protein -> stats/PCA
#---------------------------------------------------------------------------------------------
# Raw Impute strategy only
data_fill_raw <- fill_only(data_to_normalize, data_to_normalize, "Fill_Raw")
data_fill_raw_final <- dostats(data_fill_raw, "Fill_Raw")

# SL - global scaling value, sample loading normalization I * (avgSumAvgI/AvgI)
data_sl <- sl_normalize(data_to_normalize, data_to_normalize, "SL_Norm")
data_sl_final <- dostats(data_sl, "SL_Norm")

# TMM Normalized - trimmed mean normalization, coded for middle 80 (-10 top and bottom)
data_tmm <- tmm_normalize(data_to_normalize, data_to_normalize, "TMM_Norm")
data_tmm_final <- dostats(data_tmm, "TMM_Norm")

# SLTMM Normalized - SL first then trimmed mean normalization, coded for middle 80 (-10 top and bottom)
data_sltmm <- sltmm_normalize(data_to_normalize, data_to_normalize, "SLTMM_Norm")
data_sltmm_final <- dostats(data_sltmm, "SLTMM_Norm")

# TI - global scaling value, sample loading normalization I * (medianSumI/SumI)
if (use_ti){
  data_ti <- ti_normalize(data_to_normalize, data_to_normalize, "TI_Norm")
  data_ti_final <- dostats(data_ti, "TI_Norm")
}

# MI - global scaling value, sample loading normalization  I * (meanMedianI/MedianI)
if (use_mi){
  data_mi <- mi_normalize(data_to_normalize, data_to_normalize, "MI_Norm")
  data_mi_final <- dostats(data_mi, "MI_Norm")
}
# AI - global scaling value, sample loading normalization  I * (meanMeanI/MeanI)
if (use_ai){
  data_ai <- ai_normalize(data_to_normalize, data_to_normalize, "AI_Norm")
  data_ai_final <- dostats(data_ai, "AI_Norm")
}

# LR Normalized  - simple linear regression normalization
if (use_lr){
  data_lr <- lr_normalize(data_to_normalize, "LinReg_Norm")
  data_lr_final <- dostats(data_lr, "LinReg_Norm")
}

# VSN Normalized
if (use_vsn){
  data_vsn <- vsn_normalize(data_to_normalize, "VSN_Norm")
  data_vsn_final <- dostats(data_vsn, "VSN_Norm")
}

# Quantile Normalized
if (use_quantile){
  data_quantile <- quantile_normalize(data_to_normalize, "Quantile_Norm")
  data_quantile_final <- dostats(data_quantile, "Quantile_Norm")
}
# LOESS Normalized
if (use_loess){
  data_loess <- loess_normalize(data_to_normalize, "LOESS_Norm")
  data_loess_final <- dostats(data_loess, "LOESS_Norm")
}

# Local SL - local scaling value, sample loading normalization
if (use_localsl){
  data_localsl <- sl_local_normalize(data_to_normalize, data_to_normalize, "LocalSL_Norm")
  data_localsl_final <- dostats(data_localsl, "LocalSL_Norm")
}

# ProteinNorm - Specific Protein Normalized Data, ie carboxylases
if (normalize_to_protein) {data_protein_norm <- protein_normalize(data_to_normalize, "Protein_Norm")
  data_protein_norm_final <- dostats(data_protein_norm, "Protein_Norm")}


# Raw only - fill with area floor
data_raw_floor <- data_raw
data_raw_floor[data_raw_floor==0] <- area_floor
data_raw_final <- dostats(data_raw_floor, "Raw")


cv_stats(summary_cv, total_cv)

print(Sys.time() - start_time)





