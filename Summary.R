
for(i in 1:group_number) 
{
  test_average <- cbind(mean(data_ready_final[, group_cv[i]]))
}

iterations = 1
variables = group_number

output <- matrix(ncol=1, nrow=3)

for(i in 1:group_number){
  output[i,] <- mean(data_ready_final[, group_cv[i])
  
}

output

output <- data.frame(output)
class(output)