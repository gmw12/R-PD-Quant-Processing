

# fix missings - by replicate group if more than half of the values exist, replace missings with average, 
#if less than half values exist all values go to area_floor
missing_average <- function(data_in){
  data_in[data_in==0] <- NA
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(data_in[c(group_startcol[i]:group_endcol[i])]))
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- sample_groups$Count[i]
    df$min <- sample_groups$Count[i]/2
    df$missings <- rowSums(df[1:sample_groups$Count[i]] == "0")
    df$average <- df$sum / (df$rep - df$missings)
    
    for (j in 1:nrow(data_in)){
      for (k in 1:sample_groups$Count[i]){
        if (df[j,k] == "0" && df$missings[j] <= df$min[j]) { df[j,k] = df$average[j]}
        else if (df$missings[j] > df$min[j]) { df[j,k] = area_floor}
      }
    }
    assign(group_list[i], df[1:sample_groups$Count[i]])
  }
  for(i in 1:group_number)  {data_ready <- cbind(data_ready, get(group_list[i]))}
  data_ready2 <- data_ready[(sample_number+1):(sample_number*2)]
  return(data_ready2)
}


# fix missings - by replicate group if more than half of the values exist, replace missings with average, if less than half values exist all values go to area_floor
missing_minimum <- function(df){
  df[df==0] <- NA
  df$minimum <- apply(df, 1, FUN = function(x) {min(x[x > 0])})
  for (j in 1:nrow(df)){
    for (k in 1:sample_number){
      if (df[j,k] == "0") {df[j,k] = df$minimum[j]}
    }
  }
  return(df[1:sample_number])
}

# Imputing missing values using the EM algorithm proposed in section 5.4.1 of Schafer (1997).
missing_mle <- function(df){
  require(imp4p)
  df[df==0] <- NA
  df <- log2(df)
  df_mle <- impute.mle(df, group_factor) 
  df_mle <- data.frame(df_mle)
  return(df_mle)
}

# imputation of missing data
missing_fill <- function(data_in){
  data_in <- log(data_in,2)
  data_in[data_in==-Inf] = 0
  count_align <-0
  count_onemissing <- 0
  count_censor <- 0
  for(i in 1:group_number){
    # calculate stats for each sample group
    assign(group_list[i], data.frame(data_in[c(sample_groups$start[i]:sample_groups$end[i])]))
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- sample_groups$Count[i]
    df$min <- sample_groups$Count[i]/2
    df$missings <- rowSums(df[1:sample_groups$Count[i]] == 0.0)
    df$average <- apply(df[1:sample_groups$Count[i]], 1, FUN = function(x) {mean(x[x > 0])})
    df$sd <- apply(df[1:sample_groups$Count[i]], 1, FUN = function(x) {sd(x[x > 0])})
    df$bin <- ntile(df$average, 20)  
    sd_info <- subset(df, missings ==0) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
    for (x in 1:19){sd_info$max[x] <- sd_info$min[x+1]}
    sd_info$max[nrow(sd_info)] <- 100
    sd_info$min2 <- sd_info$min
    sd_info$min2[1] <- 0
    sd_info <- sd_info[-21,]
    
    # if the number of missing values <= minimum then will impute based on normal dist of measured values
    if (impute_method == "Duke"){  
      for (j in 1:nrow(data_in)){
        for (k in 1:sample_groups$Count[i]){
          if (df[j,k] == 0.0 && df$missings[j] <= df$min[j]) {
            findsd <- sd_info %>% filter(df$average[j] >= min2, df$average[j]<= max)
            nf <-  rnorm(1, 0, 1)
            testsd <- findsd$sd 
            df[j,k] = df$average[j] + (nf*findsd$sd)
            count_onemissing <- count_onemissing+1
           }
         }
       }
    }
    
    # if number of missing greater than minimum and measured value is above intensity cuttoff then remove measured value
    if (misaligned_filter){
      for (j in 1:nrow(data_in)){
        for (k in 1:sample_groups$Count[i]){
          if (df[j,k] > 1 && df$missings[j] > df$min[j]) {
            if (df[j,k] > log(intensity_cutoff,2)) {
              df[j,k] = 0.0
              count_align <- count_align+1}
            }
          }
        }
      }

    # missing > minimum will impute at background level, bottom 5%, random normal distribution
    for (j in 1:nrow(data_in)){
      for (k in 1:sample_groups$Count[i]){
        if (df[j,k] == 0.0 && df$missings[j] > df$min[j]) { 
          df[j,k] = runif(1, sd_info$min[1], sd_info$max[1])
          count_censor <- count_censor +1}
      }
    }
    
    #save group dataframe
    assign(group_list[i], df[1:sample_groups$Count[i]])
  }
  
  #saving original data with imputed data frame for trouble shooting
  df2 <- data_in
  for(i in 1:group_number)  {df2 <- cbind(df2, get(group_list[i]))}
  df3 <- df2[(sample_number+1):(sample_number*2)]
  
  if (impute_method == "LocalLeastSquares"){df3 <- lls_fill(df3)}
  if (impute_method == "KNN"){df3 <- knn_fill(df3)}  
  df3 <- data.frame(2^df3)
  df3[df3 ==1 ] <- 0
  return(df3)
}


# Local least squares imputation (lls)
lls_fill <- function(data_in){
  column_names <- colnames(data_in)
  data_in[data_in==0] <- NA
  tdata_in <-as.data.frame(t(data_in))
  tdata_in_impute <- llsImpute(tdata_in, k=150, allVariables = TRUE)
  data_out <- tdata_in_impute@completeObs
  data_out <- as.data.frame(t(data_out))
  colnames(data_out) <- column_names
  return(data_out)
}


# Local least squares imputation (lls)
knn_fill <- function(data_in){
  require(impute)
  knndata <- log2(data_in)
  knndata[knndata==-Inf] = NA
  knndata[knndata==0] <- NA
  knndata <- data.matrix(knndata)
  knnimpute <- impute.knn(knndata, k=10)
  data_out <- knnimpute$data
  return(data_out)
}