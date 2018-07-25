#-------------------------------------------
# User input
#-------------------------------------------

psm_input <- FALSE
psm_to_peptide <- FALSE
phos_peptide_only <- FALSE
peptide_to_protein <- FALSE
normalize_to_protein <- TRUE
area_floor <- 1000
pvalue_cutoff <- 0.05
fc_cutoff <- 2
color_choices <- c("red", "green", "blue", "yellow", "grey", "orange", "purple", "black")
file_prefix <- "5026_072518"
# excel export from PD
forward_data <- read_excel("5026_Protein_072518.xlsx", 1)
# excel export from PD for decoy data
# decoy_data <- read_excel("4983_MS2_decoyPSM_041918.xlsx", 1)
# excel file contains 3 columns, Output order from PD (Protein/peptide order on export, PSM numerical order), ID, Group
sample_info <- read_excel("5026_SampleList_072518.xlsx", 1)

# Protein accession numbers for plots and can normalize by 
adh_list <- c("P00330")
bait_list <- c("Q8JZL2","P30935")
avidin_list <- c("P02701")
carbox_list <- c("Q05920", "Q91ZA3", "Q5SWU9")
protein_norm_list <- c("Q05920", "Q91ZA3", "Q5SWU9") #use these if normalize_to_protein <- TRUE
# carbox_list <- c("Q05920", "Q91ZA3")
bira_list <- c("O66837")






# number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file

sample_columns <- ncol(sample_info)
sample_number <- nrow(sample_info)
comp_number <- sample_columns - 3 #subtract 3 non comparison columns

#count number of samples in each group - add column
sample_info$Count <- sapply(sample_info$Group, function(string) sum(string==sample_info$Group))

# create unqiue group dataframe with sample number
sample_groups <- sample_info[3:(sample_columns+1)]
sample_groups <- sample_groups %>% distinct(Group, .keep_all = TRUE)

# comparison groups in pairs, c(3,2) 3/2, c(3,2,4,2) 3/2, 4/2
group_comp <- vector(mode="numeric", length=0)
for (i in 2:(comp_number+1)){
group_comp <- c(group_comp, grep("N", sample_groups[[i]], ignore.case=TRUE))
group_comp <- c(group_comp, grep("D", sample_groups[[i]], ignore.case=TRUE))
}

group_number <- nrow(sample_groups)
excel_order <- sample_info$PD_Order
sample_list <- sample_info$ID
group_list <- sample_groups$Group
group_rep <- sample_groups$Count
treatment_groups <- sample_info$Group
group_color <- color_choices[1:group_number]
sample_groups$colorlist <- color_choices[1:group_number]
sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info$Group, Group)])
sample_groups$title <- str_c(sample_groups$Group,"(",sample_groups$colorlist,")")

# create a list of groups to be compared, fold change (inverse for neg), foldchange, pval
comp_fc_groups <- NULL
comp_fc2_groups <- NULL
comp_pval_groups <- NULL
comp_header <- NULL
z <- 1
for(i in 1:comp_number) {
  fc1 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC", collapse = " ")
  fc2 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC2", collapse = " ")
  pval1 <- str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_pval", collapse = " ")
  comp_fc_groups <- c(comp_fc_groups, fc1 )
  comp_fc2_groups <- c(comp_fc2_groups, fc2)
  comp_pval_groups <- c(comp_pval_groups, pval1)
  comp_header <- c(comp_header, fc1, fc2, pval1)
  z <- z+2
}


sample_groups$cv <- str_c(sample_groups$Group, "_CV")  
sample_groups$log <- str_c(sample_groups$Group, "_log")
#need to replace group_cv with sample_groups$cv and log
group_cv <- sample_groups$cv
group_log <- sample_groups$log

color_list<- sample_info$colorlist

group_title <- paste(sample_groups$title, collapse=" ")
group_title <- gsub("\\)", "\\),",  group_title)

#assign start and end columns for each group in pd output
group_startcol <- 1
for(i in 1:(group_number-1)) group_startcol <- c(group_startcol, group_startcol[i]+group_rep[i])
group_endcol <- group_startcol + group_rep
group_endcol <- group_endcol -1


#organize column headers for final output
sample_header <- str_c(sample_list[1]," ", treatment_groups[1])
for(i in 2:sample_number){
  sample_header <- c(sample_header, str_c(sample_list[i], " ", treatment_groups[i]))
}

for(i in 1:sample_number){
  sample_header <- c(sample_header, str_c(sample_list[i], " ", treatment_groups[i], " Normalized"))
}

sample_header <- c(sample_header, group_cv, comp_header)


# create subdirectory to store csv and plots
ifelse(!dir.exists(file.path(".", "output_files")), dir.create(file.path(".", "output_files")), FALSE)
output_dir <- ".//output_files//"
file_prefix <- str_c(output_dir, file_prefix)


