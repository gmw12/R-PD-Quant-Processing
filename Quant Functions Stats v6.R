#data_in <- data_sl
#data_title <- "erasemenow"
#title <- "eraseme2"
#i=1
  
dostats <- function(data_in, data_title) {
  #generate plots for proteins of interest
  plot_dir <- str_c(output_dir, data_title, "//")
  BioID_normalize_gw(data_in, data_title, plot_dir)
  PCA_gw(data_in, data_title, plot_dir)
  heatmap_gw(data_in, data_title, plot_dir)
  data_out <- stat_test_gw(data_in, data_title, plot_dir)
  cv_stats_gw(data_out, data_title)
  #final excel output
  Final_Excel_gw(data_out, str_c(file_prefix1, "_", data_title, "_final.xlsx"))
  return(data_out)
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


pvalue_gw <- function(x,y){
  x <- log2(x)
  y <- log2(y)
  temp_pval <- rep(NA, nrow(x))
  for(i in 1:nrow(x)) 
  {
    temp_pval[i] <- ttest_gw(as.numeric(x[i,]), as.numeric(y[i,]) )   
  }
  return(temp_pval)
}

exactTest_gw <- function(x,y){
  #x <- log2(x)
  #y <- log2(y)
  et <- exactTestDoubleTail(x,y) 
  return(et)
} 
  


#x<-comp_N_data
#y<-comp_D_data
limma_gw <- function(x,y, comp_name, plot_dir){
  xy <- cbind(x,y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- model.matrix(~ 0+factor(c(rep(1,n), rep(0,d))))
  colnames(design) <- c("group1", "group2")
  contrast.matrix <- makeContrasts(group2-group1, levels=design)
  fit <- lmFit(xy,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topfit <- topTable(fit2,coef=1, sort="none", number=Inf)
  data_out <-topfit$P.Value
  limma_qq(fit2$t, comp_name, plot_dir)
  limma_ma(topfit, comp_name, plot_dir)
  limma_volcano(fit2, comp_name, plot_dir)
  return(data_out)
}


old_limma_gw <- function(x,y){
  xy <- cbind(x,y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- cbind(Grp1=1,Grp2vs1=c(rep(1,n), rep(0,d)))
  #design <- c(0,0,0,1,1,1)
  # Ordinary fit
  fit <- lmFit(xy,design)
  fit <- eBayes(fit)
  topfit <- topTable(fit,coef=2, adjust.method="BH", sort="none", number=Inf)
  data_out <-topfit$P.Value
  return(data_out)
}


#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_test_gw <- function(data_in, title, plot_dir) {
  annotate_in <- data_in[1:info_columns_final]
  data_in <- data_in[(info_columns_final+1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for(i in 1:group_number) 
  {
    stat_df[ , str_c(sample_groups$Group[i], "_CV")] <- percentCV_gw(data_in[sample_groups$start[i]:sample_groups$end[i]])
  } 
  #generate pvalue and FC for each comparison
  for(i in 1:comp_number){
    comp_N_data <- data_in[comp_groups$N_start[i]:comp_groups$N_end[i]]
    comp_D_data <- data_in[comp_groups$D_start[i]:comp_groups$D_end[i]]
    stat_df[ , comp_groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ , comp_groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ , comp_groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    #volcano plot
    volcano_gw(stat_df, comp_groups$comp_name[i], title, plot_dir)
  } 
  
  # Create tables for excel--------------------------------------------------
  stat_df <- stat_df[2:ncol(stat_df)]
  #stat_df <- stat_df[,!grepl("_FC2",names(stat_df))] #removes decimal fold change, not needed after volcano plot
  #stat_df <- stat_df[,!grepl("LimmaPval",names(stat_df))] #removes limma pval 
  data_table <- cbind(data_raw, data_in, stat_df)
  return(data_table)
}



BioID_normalize_gw <- function(data_in, data_title, plot_dir) {
  
  #--------------------------------BirA-----O66837-------------------------------
  bira_test_raw <-subset(data_in, data_in$Accession %in% bira_list)
  bira_test_raw <- bira_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(bira_test_raw, str_c("BirA_", data_title), plot_dir)
  
  #--------------------------------Casein-----P02662 P02663-------------------------------
  casein_test_raw <-subset(data_in, data_in$Accession %in% casein_list)
  casein_test_raw <- casein_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(casein_test_raw, str_c("Casein_", data_title), plot_dir)
  
  #--------------------------------Carboxylase-----Q05920, Q91ZA3-------------------------------
  carbox_test_raw <-subset(data_in, data_in$Accession %in% carbox_list)
  carbox_test_raw <- carbox_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(carbox_test_raw, str_c("Carbox ", data_title), plot_dir)
  
  #--------------------------------Avidin------------------------------------
  avidin_test_raw <-subset(data_in, data_in$Accession %in% avidin_list)
  avidin_test_raw <- avidin_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(avidin_test_raw, str_c("Avidin ", data_title), plot_dir)
  
  #--------------------------------Bait-------------------------------
  bait_test_raw <-subset(data_in, data_in$Accession %in% bait_list)
  bait_test_raw <- bait_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(bait_test_raw, str_c("Bait ", data_title), plot_dir)
  
  #--------------------------------Bait-------------------------------
  ko_test_raw <-subset(data_in, data_in$Accession %in% protein1_list)
  ko_test_raw <- ko_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(ko_test_raw, str_c(protein1_list, "_", data_title), plot_dir)
  
  #--------------------------------Bait-------------------------------
  ko_test_raw <-subset(data_in, data_in$Accession %in% protein2_list)
  ko_test_raw <- ko_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(ko_test_raw, str_c(protein2_list, "_", data_title), plot_dir)
  #--------------------------------Bait-------------------------------
  ko_test_raw <-subset(data_in, data_in$Accession %in% protein3_list)
  ko_test_raw <- ko_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(ko_test_raw, str_c(protein3_list, "_", data_title), plot_dir)
  #--------------------------------Bait-------------------------------
  ko_test_raw <-subset(data_in, data_in$Accession %in% protein4_list)
  ko_test_raw <- ko_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(ko_test_raw, str_c(protein4_list, "_", data_title), plot_dir)
  
  
  #--------------------------------ADH-------------------------------
  adh_test_raw <-subset(data_in, data_in$Accession %in% adh_list)
  adh_test_raw <- adh_test_raw[(info_columns_final+1):ncol(data_in)]
  barplot_gw(adh_test_raw, str_c("ADH ",data_title), plot_dir)
  
}  



cv_stats <- function(summary_cv, total_cv) {
  # create summary table for CV's for groups under different normalize conditions
  Simple_Excel(summary_cv, str_c(file_prefix1, "_Average_CV.xlsx", collapse = " "))
  
  total_cv <- total_cv[2:ncol(total_cv)]
  png(filename=str_c(output_dir, "CV_Comparison", "_boxplot.png"), width = 800, height = 600)
  data_box <- log2(total_cv)
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



#collect CV statistics
cv_stats_gw <- function(data_out, data_title){
  cv_start <- info_columns_final+sample_number+sample_number+1
  avg_cv <- colMeans(data_out[cv_start:(cv_start+group_number-1)])
  summary_cv <-get("summary_cv",.GlobalEnv)
  summary_cv <- cbind(summary_cv, avg_cv)
  colnames(summary_cv)[ncol(summary_cv)] <- data_title
  assign("summary_cv",summary_cv,.GlobalEnv)
  all_cv <- data_out[cv_start:(cv_start+group_number-1)]
  column_names <- colnames(all_cv)
  column_names <- str_c(data_title, column_names)
  colnames(all_cv) <- column_names
  total_cv <-get("total_cv",.GlobalEnv)
  total_cv <- cbind(total_cv, all_cv)
  assign("total_cv",total_cv,.GlobalEnv)
}
