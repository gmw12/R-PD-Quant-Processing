#Function List for Quant
#normfactor_data<-data_normalize 
#data_to_norm<-data_fill 
#data_title <- "eraseme"

finish_norm <-  function(data_out, excel_name, data_title){
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, excel_name, collapse = " "))
  data_out <- impute_only(data_out)
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, impute_method, " ", excel_name, collapse = " "))
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- sample_header_final_norm  
  return(data_out)
}


impute_only <-  function(data_out){
  if (impute_method == "Duke" || impute_method == "KNN" || impute_method =="LocalLeastSquares"){
    data_out <- missing_fill(data_out)
  } else if (impute_method== "Floor") {
    data_out[data_out==0] <- area_floor
  } else if (impute_method == "Average") {
    data_out <- missing_average(data_out)
  } else if (impute_method == "Minimium") {
    data_out <- missing_minimum(data_out)
  } else if (impute_method == "MLE") {
    data_out <- missing_mle(data_out)  
  } else {
    data_out[is.na(data_out)] <- 0.0}
  return(data_out)
}

# no normalization, fill missing only
fill_only <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_Fill_Raw.xlsx"
  data_out <- finish_norm(data_to_norm, excel_name, data_title)
  return(data_out)
}

# global scaling value, sample loading normalization
sl_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_SL_Norm.xlsx"
  target <- mean(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}

# average global scaling value, sample loading normalization
ai_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_AI_Norm.xlsx"
  target <- mean(colMeans(normfactor_data))
  norm_facs <- target / colMeans(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
ti_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_TI_Norm.xlsx"
  target <- median(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
mi_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_MI_Norm.xlsx"
  intensity_medians <- colMedians(data.matrix(normfactor_data))
  target <- mean(intensity_medians)
  norm_facs <- target / intensity_medians
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


# global scaling value, sample loading normalization
vsn_normalize <- function(data_to_norm, data_title){
  excel_name <- "_Peptide_VSN_Norm.xlsx"
  data_to_norm[data_to_norm==0] <- NA
  #data_to_norm <- data_to_norm * 10
  #data_to_norm <- log2(data_to_norm)
  data_out <- normalizeVSN(data.matrix(data_to_norm))
  data_out < data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[is.na(data_out)] <- 0.0
  data_out <- data_out / 10
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


# global scaling value, sample loading normalization
quantile_normalize <- function(data_to_norm, data_title){
  excel_name <- "_Peptide_Quantile_Norm.xlsx"
  #data_to_norm <- log2(data_to_norm)
  data_out <- normalize.quantiles(data.matrix(data_to_norm))
  data_out <- data.frame(data_out)
  #data_out <- data.frame(2^data_out)
  #data_out[data_out==-Inf] = 0  # fix log2 of 0  
  #data_out[is.na(data_out)] <- 0.0
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


# global scaling value, sample loading normalization
loess_normalize <- function(data_to_norm, data_title){
  excel_name <- "_Peptide_LOESS_Norm.xlsx"
  #data_to_norm <- data_to_norm * 10
  data_to_norm <- log2(data_to_norm)
  data_out <- normalizeCyclicLoess(data_to_norm, weights = NULL, span=0.7, iterations = 3, method = "fast")
  data_out <- data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[is.na(data_out)] <- 0.0
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


# global scaling value, sample loading normalization
sl_local_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_LocalSL_Normalized.xlsx"
  data_out <- data_to_norm[ ,1]
  for(i in 1:group_number) {
    target <- mean(colSums(normfactor_data[sample_groups$start[i]:sample_groups$end[i]]))
    norm_facs <- target / colSums(normfactor_data[sample_groups$start[i]:sample_groups$end[i]])
    sl <- sweep(data_to_norm[sample_groups$start[i]:sample_groups$end[i]], 2, norm_facs, FUN = "*")   
    data_out <- cbind(data_out, sl)
  }
  data_out <- data_out[,-1]
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)  
}
 

# TMM Normalized 
tmm_normalize_old <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_TMM_Norm.xlsx"
  raw_tmm <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
  data_out <- sweep(data_to_norm, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}

#TMM Normalized 
tmm_normalize_old2 <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_TMM1_Norm.xlsx"
  if (impute_method == "Duke" || impute_method == "KNN" || impute_method =="LocalLeastSquares"){
    tmm <- missing_fill(data_to_norm)
    Simple_Excel(cbind(annotate_df, tmm), str_c(file_prefix3, "_Peptide_TMM_Norm_Impute.xlsx", collapse = " ")) 
    tmm_factor <- calcNormFactors(tmm, method = "TMM", sumTrim = 0.1)
    data_out <- sweep(tmm, 2, tmm_factor, FUN = "/") # this is data after SL and TMM on original scale
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM2_Norm.xlsx", collapse = " "))
  }else{
    tmm_factor <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
    data_out <- sweep(data_to_norm, 2, tmm_factor, FUN = "/") # this is data after SL and TMM on original scale
    Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM2_Norm.xlsx", collapse = " "))
  }
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- sample_header_final_norm
  return(data_out)
}


#TMM Normalized 
tmm_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_TMM1_Norm.xlsx"
  tmm <- impute_only(data_to_norm)
  Simple_Excel(cbind(annotate_df, tmm), str_c(file_prefix3, "_Peptide_TMM_Norm_Impute.xlsx", collapse = " ")) 
  tmm_factor <- calcNormFactors(tmm, method = "TMM", sumTrim = 0.1)
  data_out <- sweep(tmm, 2, tmm_factor, FUN = "/") # this is data after SL and TMM on original scale
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM2_Norm.xlsx", collapse = " "))
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- sample_header_final_norm
  return(data_out)
}

# SL TMM Normalized 
sltmm_normalize <- function(normfactor_data, data_to_norm, data_title){
  excel_name <- "_Peptide_SLTMM1_Norm.xlsx"
  target <- mean(colSums(normfactor_data))
  norm_facs <- target / colSums(normfactor_data)
  sl <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM1_Norm.xlsx", collapse = " "))
  sl <- impute_only(sl)
  Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM_Norm_Impute.xlsx", collapse = " ")) 
  sl_tmm <- calcNormFactors(normfactor_data, method = "TMM", sumTrim = 0.1)
  data_out <- sweep(sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SLTMM2_Norm.xlsx", collapse = " "))
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- sample_header_final_norm
  return(data_out)
  }
 

# global scaling value, sample loading normalization
protein_normalize <- function(data_to_norm, data_title){
  excel_name <- "_Peptide_ProteinNorm.xlsx"
  protein_norm_raw <- cbind(annotate_df, data_to_norm)
  protein_norm_raw <-subset(protein_norm_raw, Accession %in% protein_norm_list)
  protein_norm_raw <- protein_norm_raw[(info_columns+1):ncol(protein_norm_raw)]
  protein_norm_raw$missings <- rowSums(protein_norm_raw == 0.0)
  protein_norm_raw <- subset(protein_norm_raw, missings==0)
  protein_norm_raw <- protein_norm_raw[1:sample_number]
  target <- mean(colSums((protein_norm_raw)))
  norm_facs <- target / colSums(protein_norm_raw)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


#linear regression normalization
lr_normalize <- function(data_in, data_title) {
  excel_name <- "_Peptide_LR_Norm.xlsx"
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
  data_out <- finish_norm(data_out, excel_name, data_title)
  return(data_out)
}


#-------------------------------------------------------------------------------------
# TMT SPQC - normalize TMT, SL, IRS
#-------------------------------------------------------------------------------------
tmt_spqc_normalize <- function(){
  
  sets_tmt <- as.numeric(design_info[36,2])
  tmt_kit <- as.numeric(design_info[37,2])
  filter_peptides <- as.logical(design_info[38,2])
  sd_factor <- as.numeric(design_info[39,2])
  
  all_data <- forward_data
  
  # create subdirectory to store csv and plots
  TMT_dir <- create_dir(str_c(data_dir,"//TMT_IRS"))
  TMT_prefix <- TMT_dir #str_c(TMT_dir, "_")
  testthis <- str_c(TMT_dir, basename(input_data))
  #file.copy(input_data, str_c(TMT_dir, basename(input_data)))
  #file.copy(study_design, str_c(TMT_dir, basename(study_design)))
  
  #function Bar plot-------------------------------------------------
  barplot_TMT <- function(x,y) {
    x[is.na(x)] <- 0.0
    x <- colSums(x)
    png(filename=str_c(TMT_dir, y, "_barplot.png"), width = 888, height = 571)  
    barplot(x, 
            col = color_list,
            names.arg = treatment_groups,
            cex.names = .8,
            las=2,
            main = y)
    dev.off()
  }
  
  # function create column to flag and delete rows with no data
  delete_empty_row <- function(data_in){
    data_in$na_count <- apply(data_in, 1, function(x) sum(is.na(x)))
    data_in <- subset(data_in, data_in$na_count != tmt_kit)
    data_in <- data_in[1:(info_count + tmt_kit)]
    return(data_in)
  }
  
  #break down PD output to peptide only with protein assignments
  all_peptide <- protein_to_peptide(all_data)
  excel_order <- sample_info$PD_Order
  info_count <- ncol(all_peptide) - sample_number
  info_columns <- all_peptide[1:info_count]
  #all_peptide <-  all_peptide[(info_count+1):ncol(all_peptide)][, (excel_order)]
  all_peptide <-  all_peptide[(info_count+1):ncol(all_peptide)]
  barplot_TMT(all_peptide, "Raw_Peptide")
  raw_peptide <- cbind(info_columns, all_peptide)
  Simple_Excel(raw_peptide, str_c(TMT_prefix, "TMT_raw_peptide.xlsx"))
  
  i=1
  # split combined data from PD into individual sets
  start_col<-1
  for(i in 1:sets_tmt){
    irs_set <- subset(sample_info, sample_info$Set==i)
    irs_set$setorder <- seq(start_col, (nrow(irs_set)+start_col-1))
    irs_set$sequence <- seq(1, nrow(irs_set))
    temp_data <- cbind(info_columns, all_peptide[ ,irs_set$PD_Order])
    temp_data <- delete_empty_row(temp_data)
    temp_data[is.na(temp_data)] <- 0.0
    assign(str_c("Raw",i,"_peptide"), temp_data)
    assign(str_c("IRS",i,"_set"), irs_set)
    treatment_groups<-irs_set$Group
    color_list<-irs_set$colorlist
    barplot_TMT(temp_data[(info_count+1):ncol(temp_data)], str_c("Raw_Peptide_",i))
    start_col <- start_col + tmt_kit
  }
  
  
  # SL normalize sets, using target created from all TMT sets
  temp_data <- get(str_c("Raw",1,"_peptide"))
  target <- colSums(temp_data[(info_count+1):ncol(temp_data)])
  for(i in 2:sets_tmt){
    temp_data <- get(str_c("Raw",i,"_peptide"))
    target <- c(target, colSums(temp_data[(info_count+1):ncol(temp_data)]))
  }
  target <- mean(target)
  i=1
  for(i in 1:sets_tmt){
    temp_data <- get(str_c("Raw",i,"_peptide"))
    irs_info <- temp_data[1:info_count]
    temp_data <- temp_data[(info_count+1):ncol(temp_data)]
    norm_facs <- target / colSums(temp_data)
    temp_sl <- sweep(temp_data, 2, norm_facs, FUN = "*")
    color_list<-irs_set$colorlist
    barplot_TMT(temp_sl, str_c("SL_Peptide_",i))
    temp_sl <- cbind(irs_info, temp_sl)
    assign(str_c("Norm",i,"_SL_peptide"), temp_sl)
    Simple_Excel(temp_sl, str_c(TMT_prefix, "Norm",i, "_SL_peptide.xlsx"))
  }
  
  #Filter peptides based on CV of SPQC if requested, collapse peptide to protein
  for(i in 1:sets_tmt){
    temp_data <- get(str_c("Norm",i,"_SL_peptide"))
    if (filter_peptides){
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
      barplot_TMT(temp_data[(info_count+1):ncol(temp_data)], str_c("Norm_Peptide_Filtered_",i))
      assign(str_c("Norm",i,"_SL_filtered_peptide"), temp_data)
      Simple_Excel(temp_data, str_c(TMT_prefix, 'Norm_filtered',i,'_peptide.xlsx'))
    }
    #collapse to protein level
    temp_protein <- collapse_peptide(temp_data)
    #barplot_TMT(temp_protein[], str_c("SL_Protein_Filtered_",i))
    assign(str_c("IRS",i,"_protein"), temp_protein)
    if (filter_peptides){
      barplot_TMT(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_Protein_Filtered_",i))
      Simple_Excel(temp_protein, str_c(TMT_prefix, "SL_Filtered",i, "_protein.xlsx"))
    }else {
      barplot_TMT(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_Protein_",i))
      Simple_Excel(temp_protein, str_c(TMT_prefix, "SL_",i, "_protein.xlsx"))
    }
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
    barplot_TMT(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_CommonProtein_Filtered_",i))
    assign(str_c("IRS",i,"_common"), subset(temp_protein, temp_protein$Accession %in% common_protein))
  }
  
  #reset info columns for protein level data
  irs_info <- data.frame(IRS1_common[1:(ncol(IRS1_common)-tmt_kit)])
  info_count <- ncol(irs_info)
  
  #create SPQC data frames for each set, with sums
  for(i in 1:sets_tmt){
    spqc_data <- get(str_c("IRS",i,"_common"))
    irs_set <- get(str_c("IRS",i,"_set"))
    spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
    spqc_set$sequence <- spqc_set$sequence + info_count
    spqc_data <- spqc_data[, (spqc_set$sequence)]
    spqc_data$sums <- rowSums(spqc_data)
    assign(str_c("SPQC",i,"_data"), spqc_data)
  }
  
  #create dataframe with SPQC sums, add column with geometric mean (of sums of each protein per set)
  spqc_all <- data.frame(SPQC1_data$sums)
  for(i in 2:sets_tmt){
    spqc_data <- get(str_c("SPQC",i,"_data"))
    spqc_all <- cbind(spqc_all, spqc_data$sums)
  }
  spqc_all$geoavg <- apply(spqc_all, 1, function(x) exp(mean(log(x))))
  
  
  for(i in 1:sets_tmt){
    spqc_data <- get(str_c("SPQC",i,"_data"))
    spqc_data$fact <- spqc_data$sums/spqc_all$geoavg
    assign(str_c("SPQC",i,"_data"), spqc_data)
  }
  
  #apply IRS norm to each set
  for(i in 1:sets_tmt){
    temp_data <- get(str_c('IRS',i,'_common'))[(info_count+1):(info_count+tmt_kit)]
    spqc_data <- get(str_c('SPQC',i,'_data'))
    temp_data <- temp_data/spqc_data$fact
    barplot_TMT(temp_data[(ncol(temp_data)-tmt_kit+1):ncol(temp_data)], str_c("SL_IRS_CommonProtein_Filtered_",i))
    assign(str_c("IRS",i,"_final"), temp_data)
    Simple_Excel(temp_data, str_c(TMT_prefix, "IRS_SL_Filtered",i, "_protein.xlsx"))
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
  
  Simple_Excel(final_SL_IRS, str_c(TMT_prefix, "TMT_SL_IRS_Norm.xlsx"))
  
  treatment_groups<-sample_info$Group
  color_list<-sample_info$colorlist
  barplot_TMT(final_SL_IRS[(info_count+1):ncol(final_SL_IRS)], "Final_SL_IRS_CommonProtein_Filtered_")
  
  return(final_SL_IRS)
}