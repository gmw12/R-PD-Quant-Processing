#global options, numbers, sig digits
options(scipen = 999)
options(digits=3)


data_raw <- my_data
#set NA to 1
data_raw[is.na(data_raw)] <- 1

# filter list, only "High", delete proteins with no ID in any sample
total_columns <- ncol(data_raw)
info_columns <- total_columns - sample_number
#data_high <-data_raw[grepl("High", data_raw$Confidence, ignore.case=TRUE),] 
data_high <- data_raw

total_row <- rowSums(data_high[(info_columns+1):total_columns])
total_row <- data_frame(total_row)
data_high <- cbind(data_high, total_row)
data_high <- subset(data_high, total_row > sample_number)  # remove all lines without data !=0)
data_high <- data_high[1:total_columns]


#----- edit column headers
col_headers <- colnames(data_raw) 
col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
col_headers <- str_replace(col_headers,"\\[", "")
col_headers <- str_replace(col_headers,"\\]", "")
#col_headers <- str_replace(col_headers, "Abundance: ", "")
#col_headers <- str_replace(col_headers, "Sample, ", "")
#col_headers <- str_replace(col_headers, "F13: ", "")
colnames(data_raw) <- col_headers
info_headers <- colnames(data_raw[1:info_columns])
final_sample_header <- c(info_headers, sample_header)

# save the annotation columns (gene symbol and protein accession) for later and remove from data frame
annotate_df <- data_high[1:info_columns]
data_high <- data_high[(info_columns+1):total_columns]
row.names(data_high) <- annotate_df$`Accession`


#---------------------------------------------
# SL Normalize Data
#--------------------------------------------
# global scaling value, sample loading normalization
target <- mean(colSums((data_high)))
norm_facs <- target / colSums(data_high)
data_high_sl <- sweep(data_high, 2, norm_facs, FUN = "*")

#---------------------------------------------
# TMM from Raw
#--------------------------------------------
raw_tmm <- calcNormFactors(data_high, method = "TMM", sumTrim = 0.1)
data_high_tmm <- sweep(data_high, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale

#---------------------------------------------
# TMM & SL Normalized Data
#--------------------------------------------
# see exactly what TMM does with SL data
sl_tmm <- calcNormFactors(data_high_sl, method = "TMM", sumTrim = 0.1)
data_high_sl_tmm <- sweep(data_high_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale


#-------------------------------------------
# Plots
#------------------------------------------

# Raw, NOT normalized ----------------------
boxplot_gw(data_high, "Raw Data")
plotMDS_gw(data_high,"Raw Data_Multidimension Scaling")
data_high_bar <- colSums(data_high)
barplot_gw(data_high_bar, "Raw Data")
plotDensities_gw(data_high, "Raw data")
PCA_gw(data_high, "Raw Data")


#--Total, Normalize Data---------------------------------
boxplot_gw(data_high_sl, "SL Normalized")
plotMDS_gw(data_high_sl,"SL Normalized Data_Multidimension Scaling")
data_high_bar <- colSums(data_high_sl)
barplot_gw(data_high_bar, "SL Normalized")
plotDensities_gw(data_high_sl, "SL Normalized")
PCA_gw(data_high_sl, "SL Normalized")


#----TMM from Raw--------------------------
boxplot_gw(data_high_tmm, "TMM Normalized")
plotMDS_gw(data_high_tmm,"TMM Normalized_Multidimension Scaling")
data_high_bar <- colSums(data_high_tmm)
barplot_gw(data_high_bar, "TMM Data")
plotDensities_gw(data_high_tmm, "TMM Normalized")
PCA_gw(data_high_tmm, "TMM Normalized")


#---TMM, SL Normalized Data-------------------------------------------
boxplot_gw(data_high_sl_tmm, "TMM SL Normalized")
plotMDS_gw(data_high_sl_tmm,"TMM SL Normalized Data_Multidimension Scaling")
data_high_bar <- colSums(data_high_sl_tmm)
barplot_gw(data_high_bar, "TMM SL Normalized")
plotDensities_gw(data_high_sl_tmm, "TMM SL Normalized")
PCA_gw(data_high_sl_tmm, "TMM SL Normalized")



#-----------------------------------------------------------------------------------------
# stats
#-----------------------------------------------------------------------------------------

data_high_final <- stat_test_gw(data_high, "Raw Data")
data_high_sl_final <- stat_test_gw(data_high_sl, "SL Normalized")
data_high_tmm_final <- stat_test_gw(data_high_tmm, "TMM Normalized")
data_high_sl_tmm_final <- stat_test_gw(data_high_sl_tmm, "TMM SL Normalized")

#fix headers
colnames(data_high_final) <- final_sample_header
colnames(data_high_sl_final) <- final_sample_header
colnames(data_high_tmm_final) <- final_sample_header
colnames(data_high_sl_tmm_final) <- final_sample_header

#--csv for large peptide output
write.csv(data.frame(data_high_final), file= str_c(file_prefix, "_final.csv", collapse = " "))
write.csv(data.frame(data_high_sl_final), file= str_c(file_prefix, "_sl_final.csv", collapse = " "))
write.csv(data.frame(data_high_tmm_final), file= str_c(file_prefix, "_tmm_final.csv", collapse = " "))
write.csv(data.frame(data_high_sl_tmm_final), file= str_c(file_prefix, "_sl_tmm_final.csv", collapse = " "))



#--- directly to excel for protein projects
wb = createWorkbook()
sheet = createSheet(wb, "Sheet 1")
addDataFrame(data.frame(data_high_sl_tmm_list[1]), sheet=sheet, startColumn=1, row.names=FALSE)
sheet = createSheet(wb, "Sheet 2")
addDataFrame(data.frame(data_high_sl_tmm_list[2]), sheet=sheet, startColumn=1, row.names=FALSE)
saveWorkbook(wb, "4227_TMT_MS3_MRM_033018_sl_tmm.xlsx")

wb = createWorkbook()
sheet = createSheet(wb, "Sheet 1")
addDataFrame(data.frame(data_high_sl_list[1]), sheet=sheet, startColumn=1, row.names=FALSE)
sheet = createSheet(wb, "Sheet 2")
addDataFrame(data.frame(data_high_sl_list[2]), sheet=sheet, startColumn=1, row.names=FALSE)
saveWorkbook(wb, "4227_TMT_MS3_MRM_032318_sl.xlsx")


test_data <- data_high
test_data <- data_high_sl
test_data <- data_high_tmm
test_data <- data_high_sl_tmm
#--------------------------------BirA-----O66837-------------------------------
bira_list <- c("O66837")
bira_test_raw <-subset(test_data, rownames(test_data) %in% bira_list)
bira_raw_bar <- colSums(bira_test_raw)
barplot(bira_raw_bar, 
        col = color_list,
        main = "Bira")
#--------------------------------Carboxylase-----Q05920, Q91ZA3-------------------------------
carbox_list <- c("Q05920", "Q91ZA3")
carbox_test_raw <-subset(test_data, rownames(test_data) %in% carbox_list)
carbox_raw_bar <- colSums(carbox_test_raw)
barplot(carbox_raw_bar, 
        col = color_list,
        main = "Carbox")

#--------------------------------Avidin------------------------------------
avidin_list <- c("P02701")
avidin_test_raw <-subset(test_data, rownames(test_data) %in% avidin_list)
avidin_raw_bar <- colSums(avidin_test_raw)
barplot(avidin_raw_bar, 
        col = color_list,
        main = "Avidin")

#--------------------------------Bait-------------------------------
bait_list <- c("Q8JZP2")
bait_test_raw <-subset(test_data, rownames(test_data) %in% bait_list)
bait_raw_bar <- colSums(bait_test_raw)
barplot(bait_raw_bar, 
        col = color_list,
        main = "Bait")

#--------------------------------ADH-------------------------------
adh_list <- c("P00330")
adh_test_raw <-subset(test_data, rownames(test_data) %in% adh_list)
adh_raw_bar <- colSums(adh_test_raw)
barplot(adh_raw_bar, 
        col = color_list,
        main = "ADH")

