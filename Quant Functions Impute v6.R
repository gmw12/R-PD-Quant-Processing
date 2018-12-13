

# fix holes - by replicate group if more than half of the values exist, replace holes with average, 
#if less than half values exist all values go to area_floor
hole_average <- function(data_in){
  for(i in 1:group_number) 
  {
    assign(group_list[i], data.frame(data_in[c(group_startcol[i]:group_endcol[i])]))
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- sample_groups$Count[i]
    df$min <- sample_groups$Count[i]/2
    df$holes <- rowSums(df[1:sample_groups$Count[i]] == "0")
    df$average <- df$sum / (df$rep - df$holes)
    
    for (j in 1:nrow(data_in)){
      for (k in 1:sample_groups$Count[i]){
        if (df[j,k] == "0" && df$holes[j] <= df$min[j]) { df[j,k] = df$average[j]}
        else if (df$holes[j] > df$min[j]) { df[j,k] = area_floor}
      }
    }
    assign(group_list[i], df[1:sample_groups$Count[i]])
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
  data_in <- log(data_in,2)
  data_in[data_in==-Inf] = 0
  count_align <-0
  count_onehole <- 0
  count_censor <- 0
  for(i in 1:group_number){
    assign(group_list[i], data.frame(data_in[c(sample_groups$start[i]:sample_groups$end[i])]))
    df <- get(group_list[i])
    df$sum <- rowSums(df)
    df$rep <- sample_groups$Count[i]
    df$min <- sample_groups$Count[i]/2
    df$holes <- rowSums(df[1:sample_groups$Count[i]] == 0.0)
    df$average <- apply(df[1:sample_groups$Count[i]], 1, FUN = function(x) {mean(x[x > 0])})
    df$sd <- apply(df[1:sample_groups$Count[i]], 1, FUN = function(x) {sd(x[x > 0])})
    df$bin <- ntile(df$average, 20)  
    sd_info <- subset(df, holes ==0) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
    for (x in 1:19){sd_info$max[x] <- sd_info$min[x+1]}
    sd_info$max[nrow(sd_info)] <- 100
    sd_info$min2 <- sd_info$min
    sd_info$min2[1] <- 0
    sd_info <- sd_info[-21,]
    
    for (j in 1:nrow(data_in)){
      for (k in 1:sample_groups$Count[i]){
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
      for (k in 1:sample_groups$Count[i]){
        if (df[j,k] > 1 && df$holes[j] > df$min[j]) {
          if (df[j,k] > log(intensity_cutoff,2)) {
            df[j,k] = 0.0
            count_align <- count_align+1}
        }
      }
    }    
    for (j in 1:nrow(data_in)){
      for (k in 1:sample_groups$Count[i]){
        if (df[j,k] == 0.0) {
          df[j,k] = runif(1, sd_info$min[1], sd_info$max[1])
          count_censor <- count_censor +1}
      }
    }    
    assign(group_list[i], df[1:sample_groups$Count[i]])
  }
  #saving original data with imputed data frame for trouble shooting
  df2 <- data_in
  for(i in 1:group_number)  {df2 <- cbind(df2, get(group_list[i]))}
  df3 <- df2[(sample_number+1):(sample_number*2)]
  df3 <- data.frame(2^df3)
  df3[df3 ==1 ] <- 0
  return(df3)
}

