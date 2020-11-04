gfs <- function(vec)
{
  #alpha_1 <- x * length(vec)
  #alpha_2 <- y * length(vec)
  alpha_1 <- 0.1 * length(vec)
  alpha_2 <- 0.2 * length(vec)
  
  ranks <- rank(-vec)
  
  new_vec <- vec
  
  new_vec[ranks[1:round(alpha_1)]] <- 1
  new_vec[ranks[round(alpha_2)+1:length(vec)]] <- 0
  new_vec[ranks[round(alpha_1)+1:(round((alpha_2 - alpha_1)/5))]] <- 0.8 #bin 1
  new_vec[ranks[(round(alpha_1 + (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.6 #bin 2
  new_vec[ranks[(round(alpha_1 + 2 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.4 #bin 3
  new_vec[ranks[(round(alpha_1 + 3 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.2 #bin 4
  new_vec[ranks[(round(alpha_1 + 4 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.1 #bin 5
  
  return(new_vec)
}

















