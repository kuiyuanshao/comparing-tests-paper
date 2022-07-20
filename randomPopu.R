#######################################################################################################################
################################# Population Simulation Under the Null Hypothesis #####################################
#######################################################################################################################
library(tidyverse)
#population >>> size of the population.
#nstrata >>> number of strata possible to make stratification.
#ncluster >>> number of clusters possible to make clustering.
#nested for loops to make differences between the strata and clusters.
#nesetd for loops to make similarity within the strata and clusters.

random_population <- function(population = 1e7, nstrata = 10, ncluster = 200, mean){
  personID <- 1:population
  strata <- 1:nstrata
  
  p_cluster1 <- abs(rnorm(ncluster))
  p_cluster1 <- p_cluster1 / sum(p_cluster1)
  cluster1 <- sample(1:ncluster, population, replace = T, prob = p_cluster1)
  
  df <- as.data.frame(cbind(personID, cluster1))
  
  df$cluster2 <- 0
  df$strata <- 0
  
  df$out1 <- 0
  df$out2 <- 0
  
  df$V1 <- 0
  df$V2 <- 0
  df$V3 <- 0
  df$V4 <- 0
  
  for (i in 1:ncluster){
    index <- which(df$cluster1 == i)
    
    df$strata[index] <- rep(sample(1:nstrata, 1), length(index))
    
    ncluster2 <- sample(10:20, 1)
    cluster2 <- paste(i, 1:ncluster2, sep = "-")
    p_cluster2 <- abs(rnorm(ncluster2))
    p_cluster2 <- p_cluster2 / sum(p_cluster2)
    
    df$cluster2[index] <- sample(cluster2, length(index), replace = T, prob = p_cluster2)
  }
  
  for (i in 1:nstrata){
    index1 <- which(df$strata == i)
    
    cluster1 <- unique(df$cluster1[index1])
    for (j in 1:length(cluster1)){
      index2 <- which(df$cluster1[index1] == cluster1[j])
      
      cluster2 <- unique(df$cluster2[index1][index2])
      for (k in 1:length(cluster2)){
        index <- which(df$cluster2[index1][index2] == cluster2[k])
        
        V1 <- abs(rnorm(2, mean = 10 * i + j + k, sd = k))
        V1 <- V1 / sum(V1)
        V1 <- rbinom(length(index), 1, V1[1])
        
        df$V1[index1][index2][index] <- V1
        
        V2 <- rnorm(length(index), mean = 10 * i + j + k, 
                    sd = sqrt(i * 100))
        V2[V1 == 1] <- V2 + rnorm(sum(V1 == 1), 
                                  mean = i + j + k, 
                                  sd = (i + j + k) / 2)
        df$V2[index1][index2][index] <- V2
        
        V3 <- rnorm(length(index), mean = 15 * i + j + k, 
                    sd = sqrt(i * 100)) + V2
        df$V3[index1][index2][index] <- V3
        
        V4 <- rnorm(length(index), mean = 15 * i + j + k, 
                    sd = sqrt(i * 100)) + V3
        df$V4[index1][index2][index] <- V4
        
        
        sd <- sqrt(i + j + k)
        out1 <- rnorm(length(index), mean = 5 * i + 2 * j + k, 
                      sd = sqrt((5 * i + 2 * j + k)))
        if (mean == 0){
          out1 <- rpois(length(index), 1)
        }else{
          out1 <- scale(out1) * sd + mean
          out1[which(out1 < 0)] <- 0
        }
        df$out1[index1][index2][index] <- as.integer(out1)
        
        p <- abs(rnorm(2, mean = i + sqrt(j) + k^0.25, sd = k))
        p <- p / sum(p)
        df$out2[index1][index2][index] <- rbinom(length(index), 1, p[1])
      }
    }
  }
  return (df)
}