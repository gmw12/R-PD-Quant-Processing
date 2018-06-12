#-------------------------------------------
# User input
#-------------------------------------------

forward_data <- read_excel("4903 Protein 061118.xlsx", 1)
#decoy_data <- read_excel("4983_MS2_decoyPSM_041918.xlsx", 1)
file_prefix <- "4903_061118"
psm_input <- FALSE
psm_to_peptide <- FALSE
phos_peptide_only <- FALSE
peptide_to_protein <- FALSE
#sample list in PD order if protein output, numerical order if PSM
sample_list <- c("38592_01","38592_02","38593_03","38332","38333","38334", "38335", "38336", "38337")
                 
#sample numbers as outputed by PD
excel_order <- c(1,2,3,4,5,6,7,8,9)
# order sorted by PD
group_list <- c("QCPool", "TTBK2", "GFP" )
group_rep <- c(3, 3, 3)
group_color <- c("red", "green", "blue")
group_comp <- c(1,2)  # comparison groups in pairs, c(3,2) 3/2, c(3,2,4,2) 3/2, 4/2
area_floor <- 1000


#-------------------------------------------------------------------------------------------------
#Number of variables are not hardcoded, below creates lists of variables needed for processing
#-------------------------------------------------------------------------------------------------
sample_number <- sum(group_rep)
group_number <- length(group_rep)
comp_number <- length(group_comp)/2

comp_fc_groups <- ""
comp_fc2_groups <- ""
comp_pval_groups <- ""
z <- 1
for(i in 1:comp_number) {
  comp_fc_groups <- c(comp_fc_groups, str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC", collapse = " "))
  comp_fc2_groups <- c(comp_fc2_groups, str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_FC2", collapse = " "))
  comp_pval_groups <- c(comp_pval_groups, str_c(group_list[group_comp[z]],"_",group_list[group_comp[z+1]],"_pval", collapse = " "))
  z <- z+2
  }
  comp_fc_groups <- comp_fc_groups[-1]
  comp_fc2_groups <- comp_fc2_groups[-1]
  comp_pval_groups <- comp_pval_groups[-1]
  
treatment_groups<-c(rep(group_list[1], each=group_rep[1]))
for(i in 2:length(group_rep)) treatment_groups <- c(treatment_groups, rep(group_list[i], each=group_rep[i]))

group_cv <- str_c(group_list[1], "_CV", collapse = " ")
group_log <- str_c(group_list[1], "_log", collapse = " ")
for(i in 2:length(group_rep)) {
  group_cv <- c(group_cv, str_c(group_list[i], "_CV", collapse = " "))
  group_log <- c(group_log, str_c(group_list[i], "_log", collapse = " "))  
}

color_list<- c(rep(group_color[1], each=group_rep[1]))
for(i in 2:length(group_rep)) color_list <- c(color_list, rep(group_color[i], each=group_rep[i]))

group_title<-c(group_list[1],"(", group_color[1],"),")
for(i in 2:group_number) group_title <- c(group_title, group_list[i],"(", group_color[i],"),")
group_title <- str_c(group_title, collapse = " ")
group_title <- gsub(" \\( ", "\\(",  group_title)
group_title <- gsub("\\ )", "\\)",  group_title)
group_title <- substr(group_title, 1, nchar(group_title)-1)

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

sample_header <- c(sample_header, group_cv)

for(i in 1:comp_number) {
  sample_header <- c(sample_header, comp_fc_groups[i], comp_fc2_groups[i], comp_pval_groups[i])
}



# create subdirectory to store csv and plots
ifelse(!dir.exists(file.path(".", "output_files")), dir.create(file.path(".", "output_files")), FALSE)
output_dir <- ".//output_files//"
file_prefix <- str_c(output_dir, file_prefix)


