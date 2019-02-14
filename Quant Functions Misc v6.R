create_dir <- function(name){
  if(dir.exists(file.path(".", name))) { unlink(file.path(".", name), recursive = TRUE, force=TRUE)}
  ifelse(!dir.exists(file.path(".", name)), dir.create(file.path(".", name)), FALSE)
  return(str_c(".//", name, "//"))
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

clean_headers <- function(col_headers){
  #----- edit column headers
  col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
  col_headers <- str_replace(col_headers, "Protein FDR Confidence: Combined", "Confidence")
  col_headers <- str_replace(col_headers, "Annotated Sequence", "Annotated_Sequence")
  col_headers <- str_replace(col_headers, "Master Protein Accessions", "Accessions")
  col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
  col_headers <- str_replace(col_headers,"\\[", "")
  col_headers <- str_replace(col_headers,"\\]", "")
  col_headers <- c(col_headers, sample_info$Header1)
  return(col_headers)
}

# Rearrange columns if raw data is psm, PD does not organize
order_columns <- function(data_in){
  #use setup data to reorder columns
  #order_list <- seq(1:sample_number)
  #order_frame <- data.frame(order_list, excel_order)
  # will rearragne based on setup sheet excel order
  #data_out <- data_in[, c(order_frame$excel_order)]
  data_out <- data_in[, (excel_order)]
  return(data_out)
}



# combines forward and reverse PSM data and returns forward only correct <=1% FDR
psm_decoy <- function(forward_data, decoy_data){
  forward_data$fdr <- rep("forward", nrow(forward_data))
  decoy_data$fdr <- rep("decoy", nrow(decoy_data))
  
  isv1 <- forward_data$`Ions Score`
  isv2 <- decoy_data$`Ions Score`
  isv3 <- c(isv1, isv2)
  
  fdr1 <- forward_data$fdr
  fdr2 <- decoy_data$fdr
  fdr3 <- c(fdr1, fdr2)
  
  forward_data[nrow(forward_data)+nrow(decoy_data),] <- NA
  forward_data$Ions_Score <- isv3
  forward_data$fdr <- fdr3
  
  rm(isv1, isv2, isv3, fdr1, fdr2, fdr3, decoy_data)
  
  forward_data <- forward_data[order(forward_data$Ions_Score),]
  
  rcount <- nrow(forward_data)
  
  for (i in 1:rcount){
    testthis <- data.frame(table(forward_data$fdr[i:rcount]))
    test_fdr <- testthis$Freq[1]/testthis$Freq[2]*100
    if (test_fdr<=1.0000) {break}
  }
  fdr_data <- forward_data[i:rcount,]
  fdr_data <- fdr_data[order(-fdr_data$Ions_Score),]
  fdr_data <- fdr_data[fdr_data$fdr=="forward",]
  fdr_data$fdr <- NULL
  fdr_data$Ions_Score <- NULL
  
  return(fdr_data)
}  





#--- collapse psm to peptide-------------------------------------------------------------
collapse_psm <- function(psm_data){
  
  psm_data$longname <- paste(psm_data$`Annotated_Sequence`, psm_data$Modifications)
  psm_data$duplicate <- duplicated(psm_data$longname)
  annotate_df$longname <- paste(psm_data$`Annotated_Sequence`, psm_data$Modifications)
  annotate_df$duplicate <- duplicated(psm_data$longname)
  
  #create data frame name for each sample
  name<- rep(NA, sample_number) 
  for (i in 1:sample_number){
    name[i]<-str_c("temp",i)
  }
  #create data frames for each sample
  for (i in 1:sample_number){
    tempframe <- cbind(psm_data$longname, psm_data[info_columns+i])
    colnames(tempframe) <- c("peptide", "abundance")
    assign(name[i], tempframe)
  }
  
  #collapse psm to peptide for each sample dataframe
  for (i in 1:sample_number){
    tempframe <- get(name[i])
    tempframe <- tempframe %>%
      group_by(peptide) %>%
      summarise(abundance = sum(abundance))
    assign(name[i],tempframe)
  }
  
  #merge sample dataframes
  merge_data <- get(name[1])
  colnames(merge_data) <- c('peptide', '1')
  for (i in 2:sample_number){
    merge_data <- merge(merge_data, get(name[i]), by.x="peptide", by.y = "peptide") 
    colnames(merge_data) <- c('peptide', seq(1:i))
  }
  
  #remove duplicates from annotation data frames, merge with data to insure order
  annotate_df <- annotate_df[!annotate_df$duplicate,]
  annotate_df$duplicate <- NULL
  final_data <- merge(annotate_df, merge_data, by.x = "longname", by.y = "peptide") 
  final_data$longname <- NULL
  annotate_df<-final_data[1:info_columns] # reassign due to sorting
  final_data2 <- final_data[(info_columns+1):ncol(final_data)]
  colnames(final_data2) <- sample_header[1:sample_number]
  
  final_data2 <- cbind(annotate_df, final_data2)
  return(final_data2)
}




#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide_old <- function(peptide_data){
  test1 <- peptide_data[ , c(2:3, (info_columns+1):ncol(peptide_data))]
  colnames(test1)[colnames(test1) == 'Master Protein Accessions'] <- 'Accessions'
  colnames(test1)[colnames(test1) == 'Master Protein Descriptions'] <- 'Descriptions'
  test1$TotalPeptide <- 1.0
  
  test2 <- test1[ ,c(1, ncol(test1), 3:(sample_number+2)    )]
  test2$Accessions <- gsub(" ", "", test2$Accessions)
  
  test3 <- test2 %>% group_by(Accessions) %>% summarise_all(funs(sum))
  
  test4 <- merge( protein_names, test3,  by.x = "Accession", by.y = "Accessions")
  
  return(test4)
}

#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide <- function(peptide_data){
  colnames(peptide_data)[5] <- "Peptides"
  peptide_data$Peptides <- 1
  peptide_data$Peptides <- as.numeric(peptide_data$Peptides)
  test1 <- peptide_data[ , -1]
  test1 <- test1[ , -3]
  test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(funs(sum))
  test2 <- ungroup(test2)
  test2 %>% mutate_if(is.numeric, ~round(., 1))
  return(test2)
} 

# create final excel documents
Simple_Excel <- function(df, filename) {
  require(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet =1, df)  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# create final excel documents
Final_Excel_gw <- function(df, filename) {
  require(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet =1, df)  
  z=1
  for(i in 1:comp_number)  {
    comp_string <- comp_groups$comp_name[i]
    #assign(comp_string, subset(df, get(comp_pval_groups[i])<=pvalue_cutoff & (get(comp_fc_groups[i])>=fc_cutoff | get(comp_fc_groups[i])<= -fc_cutoff)) )
    filtered_df <- subset(df, df[ , comp_groups$pval[i]] <= pvalue_cutoff &  
                            (df[ , comp_groups$fc[i]] >= fc_cutoff | df[ , comp_groups$fc[i]] <= -fc_cutoff)) 
    addWorksheet(wb, comp_string)
    writeData(wb, sheet = i+1, filtered_df)
    z <- z+2
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}


shiny_norm_list <- function(){
  shiny_norm <- c("data_raw", "data_fill_raw", "data_sl", "data_sltmm", "data_tmm")
  if (normalize_to_protein) {shiny_norm <-c(shiny_norm, "data_protein_norm")} 
  if (use_ti) {shiny_norm <- c(shiny_norm, "data_ti")} 
  if (use_ai) {shiny_norm <- c(shiny_norm, "data_ai")} 
  if (use_mi) {shiny_norm <- c(shiny_norm, "data_mi")} 
  if (use_lr) {shiny_norm <- c(shiny_norm, "data_lr")} 
  if (use_vsn) {shiny_norm <- c(shiny_norm, "data_vsn")} 
  if (use_quantile) {shiny_norm <- c(shiny_norm, "data_quantile")} 
  if (use_loess) {shiny_norm <- c(shiny_norm, "data_loess")}   
  if (use_localsl) {shiny_norm <- c(shiny_norm, "data_localsl")}   
  names(shiny_norm) <- shiny_norm
return(shiny_norm)
}

shiny_norm_list_final <- function(){
  shiny_norm <- c("data_raw_final", "data_fill_raw_final", "data_sl_final", "data_sltmm_final", "data_tmm_final")
  if (normalize_to_protein) {shiny_norm <-c(shiny_norm, "data_protein_norm_final")} 
  if (use_ti) {shiny_norm <- c(shiny_norm, "data_ti_final")} 
  if (use_ai) {shiny_norm <- c(shiny_norm, "data_ai_final")} 
  if (use_mi) {shiny_norm <- c(shiny_norm, "data_mi_final")} 
  if (use_lr) {shiny_norm <- c(shiny_norm, "data_lr_final")} 
  if (use_vsn) {shiny_norm <- c(shiny_norm, "data_vsn_final")} 
  if (use_quantile) {shiny_norm <- c(shiny_norm, "data_quantile_final")} 
  if (use_loess) {shiny_norm <- c(shiny_norm, "data_loess_final")}   
  if (use_localsl) {shiny_norm <- c(shiny_norm, "data_localsl_final")}   
  names(shiny_norm) <- shiny_norm
return(shiny_norm)
}

#-------------------------------------------------------------------------------------
# Use psm and peptide forward and decoy data to export peptide list with 1% FDR based on psm
#-------------------------------------------------------------------------------------
peptide_psm_set_fdr <- function(){
  psmfdr_dir <- create_dir(str_c(data_dir,"//PSM_FDR"))
  psm_prefix <- str_c(psmfdr_dir, file_prefix)
  
  file.copy(forward_psm_name, str_c(psmfdr_dir, basename(forward_psm_name)))
  file.copy(decoy_psm_name, str_c(psmfdr_dir, basename(decoy_psm_name)))
  file.copy(forward_peptide_name, str_c(psmfdr_dir, basename(forward_peptide_name)))
  file.copy(decoy_peptide_name, str_c(psmfdr_dir, basename(decoy_peptide_name)))
  
  forward_psm <- subset(forward_psm, forward_psm$`Ions Score`>10)
  decoy_psm <- subset(decoy_psm, decoy_psm$`Ions Score`>10)
  
  psm_record <- c(nrow(forward_psm), nrow(decoy_psm), nrow(forward_peptide), nrow(decoy_peptide))
  psm_record_names <- c("psm_forward", "psm_decoy", "peptide_forward", "peptide_decoy")
  
  forward_psm$fdr <- rep("forward", nrow(forward_psm))
  decoy_psm$fdr <- rep("decoy", nrow(decoy_psm))
  
  isv1 <- forward_psm$`Ions Score`
  isv2 <- decoy_psm$`Ions Score`
  isv3 <- c(isv1, isv2)
  
  fdr1 <- forward_psm$fdr
  fdr2 <- decoy_psm$fdr
  fdr3 <- c(fdr1, fdr2)
  
  combo_psm <-forward_psm  
  combo_psm[nrow(forward_psm)+nrow(decoy_psm),] <- NA
  combo_psm$Ions_Score <- isv3
  combo_psm$fdr <- fdr3
  
  combo_psm <- combo_psm[order(combo_psm$Ions_Score),]
  rcount <- nrow(combo_psm)
  
  for (i in 1:rcount) {
    testthis <- data.frame(table(combo_psm$fdr[i:rcount]))
    test_fdr <- testthis$Freq[1] / testthis$Freq[2] * 100
    if (test_fdr <= 1.0000) {
      break
    }
  }
  
  fdr_psm <- combo_psm[i:rcount, ]
  fdr_psm <- fdr_psm[order(-fdr_psm$Ions_Score), ]
  psm_record_names <- c(psm_record_names, "ForwardDecoy", "IonScore")
  psm_record <- c(psm_record, nrow(fdr_psm), min(fdr_psm$Ions_Score))
  fdr_psm <- fdr_psm[fdr_psm$fdr == "forward", ]
  fdr_psm$fdr <- NULL
  fdr_psm$Ions_Score <- NULL
  psm_record_names <- c(psm_record_names, "Forward")
  psm_record <- c(psm_record, nrow(fdr_psm))
  
  Simple_Excel(fdr_psm, str_c(psmfdr_dir, "_PSM.xlsx"))
  
  fdr_psm$Match <- str_c(fdr_psm$`Sequence`,fdr_psm$`m/z [Da]`)
  peptide_psmfdr <- forward_peptide
  peptide_psmfdr$Match <- str_c(peptide_psmfdr$`Sequence`,peptide_psmfdr$`m/z [Da] (by Search Engine): Mascot`)
  
  unique_fdr_psm <- fdr_psm$Match
  unique_fdr_psm <- unique(unique_fdr_psm)
  
  peptide_psmfdr <- subset(peptide_psmfdr, Match %in% unique_fdr_psm)
  Simple_Excel(peptide_psmfdr, str_c(psm_prefix, "_Peptide.xlsx"))
  if (phos_peptide_only){
    peptide_psmfdr <- peptide_psmfdr[grep("Phospho", peptide_psmfdr$Modifications),]
    Simple_Excel(peptide_psmfdr, str_c(psm_prefix, "_PhosPeptide.xlsx"))
  }

  # find fdr of peptides directly with ion score from psm
  psm_cuttoff <- min(fdr_psm$`Ions Score`)
  
  forward_peptide$fdr <- "forward"
  decoy_peptide$fdr <- "decoy)"
  all_peptide1 <- c(forward_peptide$fdr, decoy_peptide$fdr)
  all_peptide2 <- c(forward_peptide$`Ions Score (by Search Engine): Mascot`, decoy_peptide$`Ions Score (by Search Engine): Mascot`)
  all_peptide <- cbind(all_peptide1, all_peptide2)
  all_peptide <- subset(all_peptide, all_peptide2 > psm_cuttoff)
  
  psm_record_names <- c(psm_record_names, "PeptideForwardDecoy")
  psm_record <- c(psm_record, nrow(all_peptide))
  
  all_peptide <- data.frame(all_peptide)
  all_count <- nrow(all_peptide)
  testforward <- subset(all_peptide, all_peptide1 == "forward")
  peptideFDR <- (1 - (all_count/nrow(testforward)))*100
  
  psm_record_names <- c(psm_record_names, "PeptideForward", "PeptideFDR")
  psm_record <- c(psm_record, nrow(all_peptide), peptideFDR )
  psm_stat <- cbind(psm_record_names, psm_record)
  Simple_Excel(psm_stat, str_c(psm_prefix, "_fdr_stats.xlsx"))
  
  peptide_psmfdr$Match <- NULL
  
  return(peptide_psmfdr)
}



