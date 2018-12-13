#Function List for Quant


# Rearrange columns if raw data is psm, PD does not organize
order_columns <- function(data_in){
  #use setup data to reorder columns
  #order_list <- seq(1:sample_number)
  #order_frame <- data.frame(order_list, excel_order)
  # will rearragne based on setup sheet excel order
  #data_out <- data_in[, c(order_frame$excel_order)]
  data_out <- data_in[, (excel_order)]
  colnames(data_out) <- sample_header[1:sample_number]
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




# fix holes - by replicate group if more than half of the values exist, replace holes with average, 
#if less than half values exist all values go to area_floor
hole_average <- function(data_in){
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(data_in[c(group_startcol[i]:group_endcol[i])]))
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- group_rep[i]
    df$min <- group_rep[i]/2
    df$holes <- rowSums(df[1:group_rep[i]] == "0")
    df$average <- df$sum / (df$rep - df$holes)
    
    for (j in 1:nrow(data_in)){
      for (k in 1:group_rep[i]){
        if (df[j,k] == "0" && df$holes[j] <= df$min[j]) { df[j,k] = df$average[j]}
        else if (df$holes[j] > df$min[j]) { df[j,k] = area_floor}
      }
    }
    assign(group_list[i], df[1:group_rep[i]])
  }
  
  
  
  for(i in 1:group_number)  {data_ready <- cbind(data_ready, get(group_list[i]))}
  data_ready2 <- data_ready[(sample_number+1):(sample_number*2)]

  return(data_ready2)
}




# fix holes - by replicate group if more than half of the values exist, replace holes with average, if less than half values exist all values go to area_floor
hole_minimum <- function(df){
    
    df$minimum <- apply(df, 1, FUN = function(x) {min(x[x > 0])})
    
    for (j in 1:nrow(df)){
      for (k in 1:sample_number){
        if (df[j,k] == "0") {df[j,k] = df$minimum[j]}
      }
    }
  
  return(df[1:sample_number])
}

# imputation of missing data

hole_fill <- function(data_in){
  
  count_align <-0
  count_onehole <- 0
  count_censor <- 0
  
  for(i in 1:group_number){
    assign(group_list[i], data.frame(data_in[c(group_startcol[i]:group_endcol[i])]))
    
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- group_rep[i]
    df$min <- group_rep[i]/2
    df$holes <- rowSums(df[1:group_rep[i]] == 0.0)
    df$average <- apply(df[1:group_rep[i]], 1, FUN = function(x) {mean(x[x > 0])})
    df$sd <- apply(df[1:group_rep[i]], 1, FUN = function(x) {sd(x[x > 0])})
    df$bin <- ntile(df$average, 20)  
    
    sd_info <- subset(df, holes ==0) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
    for (x in 1:19){sd_info$max[x] <- sd_info$min[x+1]}
    sd_info$max[nrow(sd_info)] <- 100
    sd_info$min2 <- sd_info$min
    sd_info$min2[1] <- 0
    sd_info <- sd_info[-21,]
    
  
    
    for (j in 1:nrow(data_in)){
      for (k in 1:group_rep[i]){
        if (df[j,k] == 0.0 && df$holes[j] <= df$min[j]) {
          findsd <- sd_info %>% filter(df$average[j] >= min2, df$average[j]<= max)
          nf <-  rnorm(1, 0, 1)
          testsd <- findsd$sd 
          df[j,k] = df$average[j] + (nf*findsd$sd)
          count_onehole <- count_onehole+1
        }
      }
    }
    
    
    for (j in 1:nrow(data_in)){
      for (k in 1:group_rep[i]){
        if (df[j,k] > 1 && df$holes[j] > df$min[j]) {
          if (df[j,k] > log(intensity_cutoff,2)) {
            df[j,k] = 0.0
            count_align <- count_align+1}
        }
      }
    }    
    
    
    for (j in 1:nrow(data_in)){
      for (k in 1:group_rep[i]){
        if (df[j,k] == 0.0) {
          df[j,k] = runif(1, sd_info$min[1], sd_info$max[1])
          count_censor <- count_censor +1}
      }
    }    
    
    
    assign(group_list[i], df[1:group_rep[i]])
  }
  
  for(i in 1:group_number)  {data_ready <- cbind(data_ready, get(group_list[i]))}
  data_ready2 <- data_ready[(sample_number+1):(sample_number*2)]
  
  return(data_ready2)
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
collapse_peptide <- function(peptide_data){
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

 

#t.test ---------------------------------
ttest_gw <- function(x,y) {
  if (!pair_comp){
  ttest_pvalue = try(t.test(x,y, 
                            alternative="two.sided",
                            var.equal = FALSE,
                            paired = FALSE), silent=TRUE)
  }else{
    ttest_pvalue = t.test(x,y, paired = TRUE)
  }
  if (is(ttest_pvalue, "try-error")) return(NA) else return(signif((ttest_pvalue$p.value), digits = 3))
}

#cohensD ---------------------------------
cohend_gw <- function(x,y) {
  cohend_est = try(cohen.d(as.numeric(x),as.numeric(y), 
                           na.rm = TRUE, pooled=FALSE, paired=FALSE))
  if (is(cohend_est, "try-error")) return(NA) else return(signif((cohend_est$estimate), digits = 3))
}


#fold change ---------------------------------
foldchange_gw <- function(x,y) {
  if(!pair_comp){
    ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x/ave_y
}else{
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn){
    indiv_fc[i] <- (x[i]/y[i])
  }
  test <- rowMeans(indiv_fc)
}
   fc <- ifelse ((test >= 1), test, -1/test)
  return(signif(fc, digits = 3))
}

#fold change pair---------------------------------
foldchange_pair_gw <- function(x,y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn){
    indiv_fc[i] <- (x[i]/y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- ifelse ((test >= 1), test, -1/test)
  return(signif(fc, digits = 3))
}



#fold change decimal ---------------------------------
foldchange_decimal_gw <- function(x,y) {
  if(!pair_comp){
  ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x/ave_y
  }else{
    sn <- ncol(x)
    indiv_fc <- x
    for (i in 1:sn){
      indiv_fc[i] <- (x[i]/y[i])
    }
    test <- rowMeans(indiv_fc) 
  }
  fc <- test
  return(signif(fc, digits = 3))
}

#fold change pair decimal---------------------------------
foldchange_pair_decimal_gw <- function(x,y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn){
    indiv_fc[i] <- (x[i]/y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- test
  return(signif(fc, digits = 3))
}


#Percent CV ---------------------------------
percentCV_gw <- function(x) {
  ave <- rowMeans(x)
  n <- ncol(x)
  sd <- apply(x[1:n], 1, sd)
  cv <- (100 * sd / ave)
  return(signif(cv, digits = 3))
}


# plot all at once
Plot_All_gw <- function(df, y) {
  boxplot_gw(df, y)
  plotMDS_gw(df, y) #str_c(y, "Multidimension Scaling"))
  df_bar <- colSums(df)
  barplot_gw(df_bar, y)
  plotDensities_gw(df, y)
}


#Box plot-------------------------------------------------
boxplot_gw <- function(x,y) {
  png(filename=str_c(output_dir, y, "_boxplot.png"), width = 888, height = 571)
  data_box <- log2(x)
  data_box[data_box ==-Inf ] <- NA
  boxplot(data_box, 
          col = color_list, 
          notch = TRUE, 
          boxwex = 0.8,
          main = c(y, group_title),
          axes=TRUE,
          horizontal = TRUE,
          las=1,
          par(mar=c(8,10,4,2)))
  dev.off()
}

#MDS Plot-------------------------------------------------
plotMDS_gw <- function(x,y) {
  png(filename=str_c(output_dir, y, "_MDS.png"), width = 888, height = 571)  
  plotMDS(log2(x), 
        col = color_list,
        main = y)
 dev.off()
}

#Bar plot-------------------------------------------------
barplot_gw <- function(x,y) {
  png(filename=str_c(output_dir, y, "_barplot.png"), width = 888, height = 571)  
  barplot(x, 
          col = color_list,
          names.arg = treatment_groups,
          cex.names = .8,
          las=2,
          main = y)
  dev.off()
}

plotDensities_gw <- function(x,y) {
  png(filename=str_c(output_dir, y, "_Density.png"), width = 888, height = 571)  
  plotDensities(log2(x), 
                group = treatment_groups,   
                col = group_color, 
                main = y)
  dev.off()
}

?snapshotPCA3d

#PCA 2D 3D-------------------------------------------------
PCA_gw <- function(x,y) {
    require(pca3d)
    require(rgl)
    x_transpose <- t(x)
    x_transpose <-data.frame(x_transpose)
    row.names(x_transpose) <- NULL
    x_transpose <-cbind(treatment_groups, x_transpose)
    x_pca <- prcomp(x_transpose[,-1], scale=TRUE)
    test_this <-x_transpose[,1]
    x_gr <- factor(unlist(test_this))
    summary(x_gr)
    pca3d(x_pca, 
        group=x_gr,
        legend = "right",
        radius = 2,
        title = y)
    snapshotPCA3d(file=str_c(output_dir, y, "_PCA3d", ".png"))
    
    df_out <- as.data.frame(x_pca$x)
    file_name <- str_c(output_dir, y, "_PCA2D.png")
    ggplot(df_out,aes(x=PC1,y=PC2,color=x_gr )) +
      geom_point(size =2) +
      theme(legend.title=element_blank()) +
      ggtitle(y) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file_name, width=5, height=4)
    
    return("done")
}


#Volcano-------------------------------------------------
volcano_gw <- function(x,y,z,title,plottitle)
{
  file_name <- str_c(output_dir, title, "_volcano.png")
  ggplot(x, aes(log2(y), -log10(z))  ) +
    geom_point(alpha=0.4, size=2, aes(color = z)) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    scale_colour_gradient(low = "blue", high = "black") +
    ggtitle(plottitle)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")  
  #legend.text=element_text(size=10, vjust=-1), 
  #legend.title = element_text(size=10, vjust = 0.2),
  #legend.title = element_blank(),
  #legend.position=c(.9,0.8))
  ggsave(file_name, width=5, height=4)
}



#Volcano-------------------------------------------------
histogram_gw <- function(x,title,plottitle)
{
  file_name <- str_c(output_dir, title, "_histogram.png")
  testgather <- gather(x)
  testgather <- subset(testgather, testgather$value>0)
  testgather$value <- log2(testgather$value)
  ggplot(testgather, aes(value))+
    geom_histogram(bins=100, fill="blue") +
    geom_vline(aes(xintercept = log2(intensity_cutoff))) +
    ggtitle(plottitle)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") 
  ggsave(file_name, width=5, height=4)
}


#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_test_gw <- function(annotate_in, data_in, title) {
  title <- title
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(data_in[c(group_startcol[i]:group_endcol[i])]))
    assign(group_cv[i], percentCV_gw(get(group_list[i])))
  } 
  
  # Fold Change-------------------------------------
  z <- 1
  for(i in 1:comp_number){
    assign(comp_fc_groups[i], foldchange_gw(get(group_list[group_comp[z]]),get(group_list[group_comp[z+1]]) ))
    assign(comp_fc2_groups[i], foldchange_decimal_gw(get(group_list[group_comp[z]]),get(group_list[group_comp[z+1]]) ))
    z <-z+2
  }
  
  # PValue-----------------------------------------------------
  for(i in 1:group_number) {
    assign(group_log[i], log(get(group_list[i]),2))
  }
  
  z <- 1
  pval_rows <- nrow(get(group_log[group_comp[z]]))
  for(i in 1:comp_number)
  {
    temp_pval <- rep(NA, pval_rows)
    temp_comp1 <- get(group_log[group_comp[z]])
    temp_comp2 <- get(group_log[group_comp[z+1]])
    for(y in 1:pval_rows) 
    {
      temp_pval[y] <- ttest_gw(as.numeric(temp_comp1[y,]), as.numeric(temp_comp2[y,]) )   
    }
    assign(comp_pval_groups[i], temp_pval)
    z <- z +2
  }
  
  # volcano
  for(i in 1:comp_number)
  {
    newtitle <- str_c(title,"_", i)
    plottitle <- str_c(title,"_", comp_groups[i])
    volcano_gw(data_in, get(comp_fc2_groups[i]), get(comp_pval_groups[i]), newtitle, plottitle)
  }
  
  # Create tables for excel--------------------------------------------------
  data_table <- cbind(annotate_in, data_ready[(info_columns+1):ncol(data_ready)], data_in)
  for(i in 1:group_number)  {data_table <- cbind(data_table, get(group_cv[i]))}
  for(i in 1:comp_number) {
    data_table <- cbind(data_table, get(comp_fc_groups[i]), get(comp_fc2_groups[i]), get(comp_pval_groups[i]))}
  return(data_table)
}




BioID_normalize_gw <- function(test_data, title) {

  #--------------------------------BirA-----O66837-------------------------------
  bira_test_raw <-subset(test_data, test_data$Accession %in% bira_list)
  bira_test_raw <- bira_test_raw[(info_columns+1):ncol(data_ready)]
  bira_raw_bar <- colSums(bira_test_raw)
  barplot_gw(bira_raw_bar, str_c("BirA_",title))

  
  #--------------------------------Carboxylase-----Q05920, Q91ZA3-------------------------------
  carbox_test_raw <-subset(test_data, test_data$Accession %in% carbox_list)
  carbox_test_raw <- carbox_test_raw[(info_columns+1):ncol(data_ready)]
  carbox_raw_bar <- colSums(carbox_test_raw)
  barplot_gw(carbox_raw_bar, str_c("Carbox ", title))
  
  #--------------------------------Avidin------------------------------------
  avidin_test_raw <-subset(test_data, test_data$Accession %in% avidin_list)
  avidin_test_raw <- avidin_test_raw[(info_columns+1):ncol(data_ready)]
  avidin_raw_bar <- colSums(avidin_test_raw)
  barplot_gw(avidin_raw_bar, str_c("Avidin ", title))

  #--------------------------------Bait-------------------------------
  bait_test_raw <-subset(test_data, test_data$Accession %in% bait_list)
  bait_test_raw <- bait_test_raw[(info_columns+1):ncol(data_ready)]
  bait_raw_bar <- colSums(bait_test_raw)
  barplot_gw(bait_raw_bar, str_c("Bait ", title))

  
  #--------------------------------ADH-------------------------------
  adh_test_raw <-subset(test_data, test_data$Accession %in% adh_list)
  adh_test_raw <- adh_test_raw[(info_columns+1):ncol(data_ready)]
  adh_raw_bar <- colSums(adh_test_raw)
  barplot_gw(adh_raw_bar, str_c("ADH ",title))

}  



# create final excel documents
Final_Excel_gw <- function(df, filename) {
  require(openxlsx)
  tempfile <- str_c(file_prefix, filename)
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet =1, df)  
  
  z=1
  for(i in 1:comp_number)  {
    comp_string <- str_c(group_list[group_comp[z]],"_v_",group_list[group_comp[z+1]])
    assign(comp_string, subset(df, get(comp_pval_groups[i])<=pvalue_cutoff &  
                                 (get(comp_fc_groups[i])>=fc_cutoff | get(comp_fc_groups[i])<= -fc_cutoff)) )
    addWorksheet(wb, comp_string)
    writeData(wb, sheet = i+1, get(comp_string))
    z <- z+2
  }
  saveWorkbook(wb, tempfile, overwrite = TRUE)
  
}


cv_stats <- function() {
  # create summary table for CV's for groups under different normalize conditions
  cv_start <- info_columns+sample_number+sample_number+1
  #raw_avg_cv <- colMeans(data_ready_final[cv_start:(cv_start+group_number-1)])
  fill_avg_cv <- colMeans(data_ready_fill_final[cv_start:(cv_start+group_number-1)])
  sl_avg_cv <- colMeans(data_ready_sl_final[cv_start:(cv_start+group_number-1)])
  tmm_avg_cv <- colMeans(data_ready_tmm_final[cv_start:(cv_start+group_number-1)])
  sl_tmm_avg_cv <- colMeans(data_ready_sl_tmm_final[cv_start:(cv_start+group_number-1)])
  if (normalize_to_protein) {
    protein_avg_cv <- colMeans(data_ready_protein_norm_final[cv_start:(cv_start+group_number-1)])
    final_cv_avg <- rbind(fill_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv, protein_avg_cv )
  }else{
    final_cv_avg <- rbind(fill_avg_cv, sl_avg_cv, tmm_avg_cv, sl_tmm_avg_cv)
  }
  
  png(filename=str_c(output_dir, "Final_CV_Avg.png"), width = 888, height = 571)
  p<-grid.table(format(final_cv_avg, digits=3))
  dev.off()
  
  
  fill_cv <- data_ready_fill_final[cv_start:(cv_start+group_number-1)]
  colnames(fill_cv) <- paste("Fill", colnames(fill_cv), sep = "_")
  sl_cv <- data_ready_sl_final[cv_start:(cv_start+group_number-1)]
  colnames(sl_cv) <- paste("SL", colnames(sl_cv), sep = "_")
  tmm_cv <- data_ready_tmm_final[cv_start:(cv_start+group_number-1)]
  colnames(tmm_cv) <- paste("TMM", colnames(tmm_cv), sep = "_")
  sl_tmm_cv <- data_ready_sl_tmm_final[cv_start:(cv_start+group_number-1)]
  colnames(sl_tmm_cv) <- paste("SLTMM", colnames(sl_tmm_cv), sep = "_")
  if (normalize_to_protein) {
    protein_cv <- data_ready_protein_norm_final[cv_start:(cv_start+group_number-1)]
    colnames(protein_cv) <- paste("Protein", colnames(protein_cv), sep = "_")
    final_cv_all <- cbind(fill_cv, sl_cv, tmm_cv, sl_tmm_cv, protein_cv )
  }else{
    final_cv_all <- cbind(fill_cv, sl_cv, tmm_cv, sl_tmm_cv)
  }
  
  png(filename=str_c(output_dir, "CV_Comparison", "_boxplot.png"), width = 800, height = 600)
  data_box <- log2(final_cv_all)
  data_box[data_box ==-Inf ] <- NA
  boxplot(data_box, 
          col = color_choices[1:group_number], 
          notch = TRUE, 
          boxwex = 0.8,
          main = c("CV Comparison", group_title),
          axes=TRUE,
          horizontal = TRUE,
          las=1,
          par(mar=c(8,10,4,2)))
  dev.off()
}



lr_normalize <- function(data_normalize, data_ready) {
  data_lr <- data_normalize
  data_out <- data_ready
  data_out[data_out==0] <-NA
  
  for(i in 1:sample_number){
    temp <- data_lr[,i]
    colnames(temp) <- "test"
    temp <- arrange(temp, test)
    data_lr[,i] <- temp
  }
  
  colnames(data_lr) <- seq(from=1, to=sample_number)
  data_lr$avg <- apply(data_lr, 1, FUN = function(x) {mean(x[x > 0])})
  
  for(i in 1:sample_number){
    data_test <- cbind(data_lr[,i], data_lr$avg)
    colnames(data_test) <- c("x", "y")
    LMfit <- rlm(y~x, data_test, na.action=na.exclude)
    Coeffs <- LMfit$coefficients
    m <- Coeffs[2] # y = mX + b
    b <- Coeffs[1] 
    normdata <- (data_out[,i] - b) / m
    data_out[,i] <- normdata
  }
  data_out[is.na(data_out)] <- 0.0
  return(data_out)
}