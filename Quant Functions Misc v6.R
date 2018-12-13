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


