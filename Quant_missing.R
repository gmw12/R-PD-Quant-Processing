

# fix holes - by replicate group if more than half of the values exist, replace holes with average, 
#if less than half values exist all values go to area_floor


hole_fill <- function(data_in){
   
  data_in <- data_ready
  count_align <-0
  count_onehole <- 0
  count_censor <- 0

  for(i in 1:group_number){
    i=2
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
          if (df[j,k] > log(5000000,2)) {
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
  



 
  
prayer <- subset(df, holes==0)
test_d <- density(prayer$average)
plot(test_d)   

prayer2 <- subset(df, holes==1)
test_d <- density(prayer2$average)
plot(test_d)  

prayer3 <- subset(df, holes==2)
test_d <- density(prayer3$average)
plot(test_d)



prayer_sd <-sd(prayer$average)
