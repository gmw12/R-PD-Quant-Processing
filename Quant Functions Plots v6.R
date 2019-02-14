
# plot all at once
Plot_All_gw <- function(df, y) {
  plot_dir <- create_dir(str_c(output_dir, y))
  boxplot_gw(df, y, plot_dir)
  plotMDS_gw(df, y, plot_dir) #str_c(y, "Multidimension Scaling"))
  barplot_gw(df, y, plot_dir)
  barplot_old_gw(df, y, plot_dir)
  #barplot_log2_gw(df_bar, y)
  plotDensities_gw(df, y, plot_dir)
}


#Student t qq plot -------------------------------------------------
limma_qq <- function(x, comp_name, plot_dir) {
  png(filename=str_c(plot_dir, comp_name, "_StudentTQQ.png"), width = 888, height = 571)  
    qqt(x, df=Inf, pch=16,cex=0.2, main=str_c("Student's t Q-Q Plot", " ", comp_name))
    abline(0,1) 
  dev.off()
}

#Student t qq plot -------------------------------------------------
limma_volcano <- function(x, comp_name, plot_dir) {
  png(filename=str_c(plot_dir, comp_name, "_LimmaVolcano.png"), width = 888, height = 571)  
  volcanoplot(x, coef = 1, highlight = 8, names=data_raw$Accession, main= comp_name)
  dev.off()
}

#MA plot-------------------------------------------------
limma_ma <- function(x, comp_name, plot_dir) {
  df <- data.frame(data_raw$Accession)
  colnames(df) <- "Accession"
  df$baseMean <- x$AveExpr
  df$log2FoldChange <- x$logFC
  df$padj <- x$P.Value
  file_name <- str_c(plot_dir, comp_name, "_MAplot.png")
  ggmaplot(df, main = comp_name,
           fdr = pvalue_cutoff, fc = fc_cutoff, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = as.vector(df$Accession),
           legend = "top", top = 20,
           font.label = c("bold", 11), label.rectangle = TRUE,
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())
  ggsave(file_name, width=5, height=4)
}






#Box plot-------------------------------------------------
boxplot_gw <- function(x,y, plot_dir) {
  png(filename=str_c(plot_dir, y, "_boxplot.png"), width = 888, height = 571)
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
plotMDS_gw <- function(x,y,plot_dir) {
  png(filename=str_c(plot_dir, y, "_MDS.png"), width = 888, height = 571)  
  plotMDS(log2(x), 
          col = color_list,
          main = y)
  dev.off()
}

#Bar plot-------------------------------------------------
barplot_old_gw <- function(df,y,plot_dir) {
  namex <- sample_info$Label
  datay <- colSums(df)
  df2 <- data.frame(cbind(namex, datay))
  colnames(df2) <- c("Sample", "Total_Intensity")
  df2$Sample <- factor(df2$Sample, levels = df2$Sample)
  file_name <- str_c(plot_dir, y, "_barplot.png")
  
  ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
    geom_bar(stat="identity", fill=sample_info$colorlist)+ theme_classic() + 
    ggtitle(y) + 
    scale_y_discrete(labels = NULL) +
    coord_cartesian(ylim=NULL, expand=FALSE) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=6)) 
  ggsave(file_name, width=5, height=4)
  return("done")
 }


#Bar plot-------------------------------------------------
barplot_gw <- function(x,y,plot_dir) {
  x <- colSums(x)
  png(filename=str_c(plot_dir, y, "_barplot2.png"), width = 888, height = 571)  
  barplot(x, 
          col = color_list,
          names.arg = treatment_groups,
          cex.names = .8,
          las=2,
          main = y)
  dev.off()
}


#Bar plot-------------------------------------------------
barplot_adh <- function(x,name,y,plot_dir) {
  png(filename=str_c(plot_dir, y, "_barplot2.png"), width = 888, height = 571)  
  barplot(x, 
          col = rainbow(n=12),
          names.arg = name,
          cex.names = 0.8,
          las=2,
          main = y)
  dev.off()
}

#Bar plot-------------------------------------------------
barplot_log2_gw <- function(x,y,plot_dir) {
  x <- colSums(x)
  x<-log(x,2)
  png(filename=str_c(plot_dir, y, "_barplot_log2.png"), width = 888, height = 571)  
  barplot(x, 
          col = color_list,
          names.arg = treatment_groups,
          cex.names = .8,
          las=2,
          main = y)
  dev.off()
}

plotDensities_gw <- function(x,y,plot_dir) {
  png(filename=str_c(plot_dir, y, "_Density.png"), width = 888, height = 571)  
  plotDensities(log2(x), 
                group = treatment_groups,   
                col = group_color, 
                main = y)
  dev.off()
}

?pca3d
#PCA 2D 3D-------------------------------------------------
PCA_gw <- function(x,y, plot_dir) {
  x <- x[(info_columns_final+1):ncol(x)] # strip off info columns
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
  snapshotPCA3d(file=str_c(plot_dir, y, "_PCA3d", ".png"))
  
  df_out <- as.data.frame(x_pca$x)
  file_name <- str_c(plot_dir, y, "_PCA2D.png")
  ggplot(df_out,aes(x=PC1,y=PC2,color=x_gr )) +
    geom_point(size =3) +
    theme(legend.title=element_blank()) +
    ggtitle(y) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file_name, width=5, height=4)
  return("done")
}


#Volcano-------------------------------------------------
volcano_gw <- function(x, comp_name, title, plot_dir)
{
  newtitle <- str_c(title,"_", i)
  plottitle <- str_c(title,"_", comp_name)
  y <- unlist(x[ , str_c(comp_name,"_FC2")])
  z <- unlist(x[ , str_c(comp_name, "_Pval")])
  file_name <- str_c(plot_dir, title, "_", comp_name, "_volcano.png")
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



#Histogram for total intensity-------------------------------------------------
histogram_gw <- function(x,title,plottitle)
{
  testthis1 <- as.matrix(log2(x))
  x_mean <- mean(testthis1, na.rm=TRUE)
  x_stdev <- sd(testthis1, na.rm=TRUE)
  
  file_name <- str_c(output_dir, title, "_histogram.png")
  testgather <- gather(x)
  testgather <- subset(testgather, testgather$value>0)
  testgather$value <- log2(testgather$value)
  ggplot(testgather, aes(value))+
    geom_histogram(bins=100, fill="blue") +
    geom_vline(aes(xintercept = log2(intensity_cutoff))) +
    geom_vline(aes(xintercept = x_mean)) +
    geom_vline(aes(xintercept = (x_mean + x_stdev))) +
    geom_vline(aes(xintercept = (x_mean - x_stdev))) +
    ggtitle(plottitle)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") 
  ggsave(file_name, width=5, height=4)
  x_mean <- 2^x_mean
  if (x_mean < intensity_cutoff){intensity_cutoff<-x_mean}
  return(intensity_cutoff)
}


#Heat map-------------------------------------------------
heatmap_gw <- function(y,plottitle,plot_dir)
{
  y <- y[(info_columns_final+1):ncol(y)] # strip off info columns
  y <- log2(y)
  y <- data.matrix(y)
  ## Row- and column-wise clustering 
  hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
  hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
  ## Tree cutting
  mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
  ## Plot heatmap 
  mycol <- redgreen(75) #colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
  png(filename=str_c(plot_dir, plottitle, "_heatmap.png"), width = 888, height = 571)  
  heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, labCol= treatment_groups, 
            scale="row", density.info="none", trace="none", RowSideColors=mycolhc, main = plottitle) 
  dev.off()
}

