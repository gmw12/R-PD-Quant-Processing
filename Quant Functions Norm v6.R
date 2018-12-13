#Function List for Quant
#normfactor_data<-data_normalize 
#data_to_norm<-data_fill 
#data_title <- "eraseme"

# global scaling value, sample loading normalization
sl_normalize <- function(normfactor_data, data_to_norm, data_title){
  target <- mean(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SL_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SL_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


# average global scaling value, sample loading normalization
ai_normalize <- function(normfactor_data, data_to_norm, data_title){
  target <- mean(colMeans(normfactor_data))
  norm_facs <- target / colMeans(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_AI_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_AI_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
ti_normalize <- function(normfactor_data, data_to_norm, data_title){
  target <- median(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TI_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TI_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
mi_normalize <- function(normfactor_data, data_to_norm, data_title){
  intensity_medians <- colMedians(data.matrix(normfactor_data))
  target <- mean(intensity_medians)
  norm_facs <- target / intensity_medians
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_MI_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_MI_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


# global scaling value, sample loading normalization
vsn_normalize <- function(data_to_norm, data_title){
  data_to_norm <- data_to_norm * 10
  data_to_norm <- log2(data_to_norm)
  data_out <- justvsn(data.matrix(data_to_norm))
  data_out < data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[is.na(data_out)] <- 0.0
  data_out <- data_out / 10
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_VSN_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_VSN_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


# global scaling value, sample loading normalization
quantile_normalize <- function(data_to_norm, data_title){
  #data_to_norm <- log2(data_to_norm)
  data_out <- normalize.quantiles(data.matrix(data_to_norm))
  data_out <- data.frame(data_out)
  #data_out <- data.frame(2^data_out)
  #data_out[data_out==-Inf] = 0  # fix log2 of 0  
  #data_out[is.na(data_out)] <- 0.0
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_Quantile_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_Quantile_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


# global scaling value, sample loading normalization
loess_normalize <- function(data_to_norm, data_title){
  #data_to_norm <- data_to_norm * 10
  data_to_norm <- log2(data_to_norm)
  data_out <- normalizeCyclicLoess(data_to_norm, weights = NULL, span=0.7, iterations = 3, method = "fast")
  data_out <- data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[is.na(data_out)] <- 0.0
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LOESS_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LOESS_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}



# global scaling value, sample loading normalization
sl_local_normalize <- function(normfactor_data, data_to_norm, data_title){
  data_out <- data_to_norm[ ,1]
  for(i in 1:group_number) {
    target <- mean(colSums(normfactor_data[group_startcol[i]:group_endcol[i]]))
    norm_facs <- target / colSums(normfactor_data[group_startcol[i]:group_endcol[i]])
    sl <- sweep(data_to_norm[group_startcol[i]:group_endcol[i]], 2, norm_facs, FUN = "*")   
    data_out <- cbind(data_out, sl)
  }
  data_out <- data_out[,-1]
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LocalSL_Normalized.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LocalSL_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)  
}
 

# TMM Normalized 
tmm_normalize <- function(normfactor_data, data_to_norm, data_title){
  raw_tmm <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
  data_out <- sweep(data_to_norm, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


# TMM Normalized 
sltmm_normalize <- function(normfactor_data, data_to_norm, data_title){
  target <- mean(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  sl <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM1_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){
    sl <- hole_fill(sl)
    Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM_Norm_Impute.xlsx", collapse = " ")) 
    sl_tmm <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
    data_out <- sweep(sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SLTMM2_Norm.xlsx", collapse = " "))
    }else{
      sl_tmm <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
      data_out <- sweep(sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
      Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SLTMM2_Norm.xlsx", collapse = " "))
    }
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
  }


 

# global scaling value, sample loading normalization
protein_normalize <- function(data_to_norm, data_title){
  protein_norm_raw <- cbind(annotate_df, data_to_norm)
  protein_norm_raw <-subset(protein_norm_raw, Accession %in% protein_norm_list)
  protein_norm_raw <- protein_norm_raw[(info_columns+1):ncol(protein_norm_raw)]
  protein_norm_raw$holes <- rowSums(protein_norm_raw == 0.0)
  protein_norm_raw <- subset(protein_norm_raw, holes==0)
  protein_norm_raw <- protein_norm_raw[1:sample_number]
  target <- mean(colSums((protein_norm_raw)))
  norm_facs <- target / colSums(protein_norm_raw)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_ProteinNorm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_Protein_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}


#linear regression normalization
lr_normalize <- function(data_in, data_title) {
  #normalize lr on data with no missing values, create new data frame
  data_nomissing <- data_in
  data_nomissing$missingvalues <- rowSums(data_nomissing ==0)
  data_nomissing <- subset(data_nomissing, missingvalues==0)
  data_nomissing <- data_nomissing[1:sample_number]
  #log2 data
  data_nomissing <- log(data_nomissing,2)
  data_nomissing[data_nomissing==-Inf] = 0  # fix log2 of 0}
  data_in <- log(data_in,2)
  data_in[data_in==-Inf] = 0  # fix log2 of 0}  
  
  data_in[data_in==0] <- NA
  data_out <- data_in
  
  #reorders data by intensity independent of indentification
  for(i in 1:sample_number){
    temp <- data_nomissing[,i]
    colnames(temp) <- "test"
    temp <- arrange(temp, test)
    data_nomissing[,i] <- temp
  }
  
  colnames(data_nomissing) <- seq(from=1, to=sample_number)
  data_nomissing$avg <- apply(data_nomissing, 1, FUN = function(x) {median(x[x > 0])})
  
  for(i in 1:sample_number){
    i=1
    data_test <- cbind(data_nomissing[,i], data_nomissing$avg)
    colnames(data_test) <- c("x", "y")
    LMfit <- rlm(x~y, data_test, na.action=na.exclude)
    Coeffs <- LMfit$coefficients
    m <- Coeffs[2] # y = mX + b
    b <- Coeffs[1] 
    normdata <- (data_in[,i] - b) / m
    data_out[,i] <- normdata
  }
  
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[is.na(data_out)] <- 0.0
  
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LR_Norm.xlsx", collapse = " "))
  if (holes == "Impute"){data_out <- hole_fill(data_out)
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_LR_Norm_Impute.xlsx", collapse = " "))}
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- final_sample_header_protein_norm
  return(data_out)
}