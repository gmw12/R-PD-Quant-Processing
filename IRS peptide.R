IRS1_data <- read_excel("4227 Cortex Syngap1 110518.xlsx", 1)
IRS2_data <- read_excel("4227 Cortex Ube3A 110518.xlsx", 1)
IRS3_data <- read_excel("4227 Cortex Shank2 110518.xlsx", 1)
IRS4_data <- read_excel("4227 Cortex Shank3 110518.xlsx", 1)

IRS1_prefix <- "4227_Cortex_Syngap1"
IRS2_prefix <- "4227_Cortex_Ube3A"
IRS3_prefix <- "4227_Cortex_Shank2"
IRS4_prefix <- "4227_Cortex_Shank3"

sets_tmt <- 4
tmt_kit <- 11

IRS1_peptide <- protein_to_peptide(IRS1_data)
IRS2_peptide <- protein_to_peptide(IRS2_data)
IRS3_peptide <- protein_to_peptide(IRS3_data)
IRS4_peptide <- protein_to_peptide(IRS4_data)


IRS1_sl <- tmt_sl_normalize(IRS1_peptide)
IRS2_sl <- tmt_sl_normalize(IRS2_peptide)
IRS3_sl <- tmt_sl_normalize(IRS3_peptide)
IRS4_sl <- tmt_sl_normalize(IRS4_peptide)


IRS1_sl$unique <- str_c(IRS1_sl$`Master Protein Accessions`, IRS1_sl$`Annotated Sequence`, IRS1_sl$Modifications, sep = ":")
IRS2_sl$unique <- str_c(IRS2_sl$`Master Protein Accessions`, IRS2_sl$`Annotated Sequence`, IRS2_sl$Modifications, sep = ":")
IRS3_sl$unique <- str_c(IRS3_sl$`Master Protein Accessions`, IRS3_sl$`Annotated Sequence`, IRS3_sl$Modifications, sep = ":")
IRS4_sl$unique <- str_c(IRS4_sl$`Master Protein Accessions`, IRS4_sl$`Annotated Sequence`, IRS4_sl$Modifications, sep = ":")

IRS1_sl$unique <- IRS1_sl$`Annotated Sequence`
IRS2_sl$unique <- IRS2_sl$`Annotated Sequence`
IRS3_sl$unique <- IRS3_sl$`Annotated Sequence`
IRS4_sl$unique <- IRS4_sl$`Annotated Sequence`

common_peptide <- Reduce(intersect, list(IRS1_sl$unique, IRS2_sl$unique, IRS3_sl$unique, IRS4_sl$unique))
common_peptide2 <- Reduce(intersect, list(IRS1_sl$unique, IRS2_sl$unique))



IRS1_irs <- subset(IRS1_sl, unique %in% common_peptide)
IRS2_irs <- subset(IRS2_sl, unique %in% common_peptide)
IRS3_irs <- subset(IRS3_sl, unique %in% common_peptide)
IRS4_irs <- subset(IRS4_sl, unique %in% common_peptide)


testthis <- unique(IRS1_irs[ , 2])


IRS1_irs <- IRS1_irs[order(IRS1_irs$Accession),]
IRS2_irs <- IRS2_irs[order(IRS2_irs$Accession),]
IRS3_irs <- IRS3_irs[order(IRS3_irs$Accession),]


IRS_info <- IRS1_sl[1:(ncol(IRS1_sl)-tmt_kit)]


IRS1_spqc <- IRS1_sl[(ncol(IRS1_sl)-tmt_kit+1):ncol(IRS1_sl)]
IRS2_spqc <- IRS2_sl[(ncol(IRS2_sl)-tmt_kit+1):ncol(IRS2_sl)]
IRS3_spqc <- IRS3_sl[(ncol(IRS3_sl)-tmt_kit+1):ncol(IRS3_sl)]



spqc <- data.frame(rowSums(IRS1_spqc), rowSums(IRS2_spqc), rowSums(IRS3_spqc))
spqc$average <- apply(spqc, 1, FUN = function(x) {mean(x)})


spqc$fact1 <- spqc$average/spqc$rowSums.IRS1_spqc.
spqc$fact2 <- spqc$average/spqc$rowSums.IRS2_spqc.
spqc$fact3 <- spqc$average/spqc$rowSums.IRS3_spqc.

IRS1_sl_irs <- IRS1_spqc * spqc$fact1
IRS2_sl_irs <- IRS2_spqc * spqc$fact2
IRS3_sl_irs <- IRS3_spqc * spqc$fact3

IRS_all_sl_irs <- cbind(IRS1_sl_irs, IRS2_sl_irs, IRS3_sl_irs)

#---------------------------------------------
# TMM Normalized 
#--------------------------------------------
raw_tmm <- calcNormFactors(IRS_all_sl_irs, method = "TMM", sumTrim = 0.1)
IRS_all_sl_irs_tmm <- sweep(IRS_all_sl_irs, 2, raw_tmm, FUN = "/") 
IRS1_sl_irs_tmm <- sweep(IRS1_sl_irs, 2, raw_tmm, FUN = "/")
IRS2_sl_irs_tmm <- sweep(IRS2_sl_irs, 2, raw_tmm, FUN = "/")
IRS3_sl_irs_tmm <- sweep(IRS3_sl_irs, 2, raw_tmm, FUN = "/")

# create final excel documents
Final_Excel_gw <- function(df, filename, description) {
  require(openxlsx)
  tempfile <- str_c(filename, description, ".xlsx" )
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet =1, df)  
  saveWorkbook(wb, tempfile, overwrite = TRUE)
  
}

Final_Excel_gw(IRS1_sl_irs, IRS1_prefix, "_sl_irs")
Final_Excel_gw(IRS2_sl_irs, IRS2_prefix, "_sl_irs")
Final_Excel_gw(IRS3_sl_irs, IRS3_prefix, "_sl_irs")




# global scaling value, sample loading normalization
tmt_sl_normalize <- function(data_in){
  data_info <- data_in[1:(ncol(data_in)-tmt_kit)]
  data_in <- data_in[(ncol(data_in)-tmt_kit+1):ncol(data_in)]
  data_in[is.na(data_in)] <- 1
  target <- mean(colSums(data_in))
  norm_facs <- target / colSums(data_in)
  data_out_sl <- sweep(data_in, 2, norm_facs, FUN = "*")
  data_out_sl <- cbind(data_info, data_out_sl)
  return(data_out_sl)
}


# input PD protein/peptide export and output peptide (keeps master protein assignment)
protein_to_peptide <- function(data_in){
  
  data_in$`Protein FDR Confidence: Mascot`[is.na(data_in$`Protein FDR Confidence: Mascot`)] <- ""  
  
  for(i in 1:nrow(data_in)) {
    if(data_in[i,1] == "High") {set_accession <- data_in[i,3]
    }else{data_in[i,3] <- set_accession}
  }
  
  protein_names <- data_in
  colnames(protein_names)[1] <- "filter"
  protein_names <- subset(protein_names, filter=="High")
  protein_names <- protein_names[ , c(3:4)]
  peptide_header <- data_in[2,]
  data_in <- subset(data_in, Master=="High")
  colnames(data_in) <- peptide_header
  data_in <- data_in[,-1]
  colnames(data_in)[colnames(data_in) == 'Quan Usage'] <- 'Used'
  data_in <- subset(data_in, Used=="Used")
  data_in <- data_in[-ncol(data_in)]
  data_in <- data_in[-ncol(data_in)]
  data_in <- data_in[-ncol(data_in)]
  data_in <- data_in[-ncol(data_in)]
  data_in[(ncol(data_in)-tmt_kit+1):ncol(data_in)] <- as.numeric(unlist(data_in[(ncol(data_in)-tmt_kit+1):ncol(data_in)]))
  colnames(data_in)[2] <- "Master Protein Accessions"
  return(data_in)
}