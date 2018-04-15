#Function List for Quant


# Rearrange columns if raw data is psm, PD does not organize
order_columns <- function(data_in){
  #use setup data to reorder columns
  order_list <- seq(1:sample_number)
  order_frame <- data.frame(order_list, excel_order)
  order_frame <- order_frame[order(order_frame$excel_order),]
  
  data_out <- data_in[, c(order_frame$order_list)]
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


#--- collapse psm to peptide-------------------------------------------------------------
collapse_psm <- function(psm_data){

  psm_data <- data_ready
  
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
  return_list <- list(final_data2, annotate_df)
  
  return(return_list)
}
 








#t.test ---------------------------------
ttest_gw <- function(x,y) {
  ttest_pvalue = try(t.test(x,y, 
                            alternative="two.sided",
                            var.equal = FALSE,
                            paired = FALSE), silent=TRUE)
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
  ave_x = rowMeans(x)
  ave_y = rowMeans(y)
  test = ave_x/ave_y
  fc <- ifelse ((test >= 1), test, -1/test)
  return(signif(fc, digits = 3))
}

#fold change decimal ---------------------------------
foldchange_decimal_gw <- function(x,y) {
  ave_x = rowMeans(x)
  ave_y = rowMeans(y)
  test = ave_x/ave_y
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

#Box plot-------------------------------------------------
boxplot_gw <- function(x,y) {
  png(filename=paste(output_dir, y, "_boxplot.png"), width = 888, height = 571)
  boxplot(log2(x), 
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
  png(filename=paste(output_dir, y, "_MDS.png"), width = 888, height = 571)  
  plotMDS(log2(x), 
        col = color_list,
        main = y)
 dev.off()
}

#Bar plot-------------------------------------------------
barplot_gw <- function(x,y) {
  png(filename=paste(output_dir, y, "_barplot.png"), width = 888, height = 571)  
  barplot(x, 
          col = color_list,
          main = y)
  dev.off()
}

plotDensities_gw <- function(x,y) {
  png(filename=paste(output_dir, y, "_Density.png"), width = 888, height = 571)  
  plotDensities(log2(x), 
                group = treatment_groups,   
                col = group_color, 
                main = y)
  dev.off()
}



#PCA 2D 3D-------------------------------------------------
PCA_gw <- function(x,y) {
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
    snapshotPCA3d(file=paste(output_dir, y, "_PCA3d", ".png"))
    
    png(filename=paste(output_dir, y, "_PCA2D.png"), width = 888, height = 571)    
      pca2d(x_pca, 
            group=x_gr, 
            legend="topleft")
    dev.off()
    
    return("done")
}

#Volcano-------------------------------------------------
volcano_gw <- function(x,y,z,title)
{
  file_name <- str_c(output_dir, title, "_volcano.png")
  ggplot(x, aes(log2(y), -log10(z))  ) +
    geom_point(alpha=0.4, size=2, aes(color = z)) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    scale_colour_gradient(low = "blue", high = "black") +
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")  
  #legend.text=element_text(size=10, vjust=-1), 
  #legend.title = element_text(size=10, vjust = 0.2),
  #legend.title = element_blank(),
  #legend.position=c(.9,0.8))
  ggsave(file_name, width=5, height=4)
}


#Percent CV-------------------------------------------------
stat_cv_gw <- function(x) {
  
  # separate by treatment
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(x[c(group_startcol[i]:group_endcol[i])]))
    temp_string <- c(group_list[i], "_CV")
    assign(temp_string, percentCV_gw(get(group_list[i])))
  }   
}

#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_test_gw <- function(x, title) {
  title <- title
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(x[c(group_startcol[i]:group_endcol[i])]))
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
      temp_pval[y] <- ttest_gw(temp_comp1[y,], temp_comp2[y,])    
    }
    assign(comp_pval_groups[i], temp_pval)
    z <- z +2
  }
  
  # volcano
  for(i in 1:comp_number)
  {
    volcano_gw(x, get(comp_fc2_groups[i]), get(comp_pval_groups[i]), title)
  }
  
  
  
  # Create tables for excel--------------------------------------------------
  data_table <- data.frame(annotate_df, data_high, x)
  
  for(i in 1:group_number) 
  {
    data_table <- data.frame(data_table, get(group_cv[i]))
  }
  
  for(i in 1:comp_number) 
  {
    data_table <- data.frame(data_table, get(comp_fc_groups[i]), get(comp_fc2_groups[i]), get(comp_pval_groups[i]))
  }

  
  return(data_table)
}

