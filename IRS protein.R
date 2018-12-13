IRS1_data <- read_excel("4227 Cortex Syngap1 110518.xlsx", 1)
IRS2_data <- read_excel("4227 Cortex Ube3A 110518.xlsx", 1)
IRS3_data <- read_excel("4227 Cortex Shank2 110518.xlsx", 1)
IRS4_data <- read_excel("4227 Cortex Shank3 110518.xlsx", 1)

IRS1_prefix <- "4227_Cortex_Syngap1"
IRS2_prefix <- "4227_Cortex_Ube3A"
IRS3_prefix <- "4227_Cortex_Shank2"
IRS4_prefix <- "4227_Cortex_Shank3"

sets_tmt <- 3
tmt_kit <- 11


# global scaling value, sample loading normalization
tmt_sl_normalize <- function(data_in){
  data_info <- data_in[1:(ncol(data_in)-tmt_kit)]
  data_in <- data_in[(ncol(data_in)-tmt_kit+1):ncol(data_in)]
  target <- mean(colSums(data_in))
  norm_facs <- target / colSums(data_in)
  data_out_sl <- sweep(data_in, 2, norm_facs, FUN = "*")
  data_out_sl <- cbind(data_info, data_out_sl)
  return(data_out_sl)
}

IRS1_sl <- tmt_sl_normalize(IRS1_data)
IRS2_sl <- tmt_sl_normalize(IRS2_data)
IRS3_sl <- tmt_sl_normalize(IRS3_data)

 
common_protein <- Reduce(intersect, list(IRS1_sl$Accession, IRS2_sl$Accession, IRS3_sl$Accession))


IRS1_sl <- subset(IRS1_data, Accession %in% common_protein)
IRS2_sl <- subset(IRS2_data, Accession %in% common_protein)
IRS3_sl <- subset(IRS3_data, Accession %in% common_protein)

IRS1_sl <- IRS1_sl[order(IRS1_sl$Accession),]
IRS2_sl <- IRS2_sl[order(IRS2_sl$Accession),]
IRS3_sl <- IRS3_sl[order(IRS3_sl$Accession),]

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



