#######################################################################################################################
###################################### Stratified Multistage Sampling #################################################
#######################################################################################################################

library(sampling)
#population >>> population to sample with.
#ncluster1 >>> total number of PSU in the first sampling stage.
#ncluster2 >>> total number of SSU in the second sampling stage
#cluster_observ >>> how many individuals being sampled in each cluster in the third sampling stage.
sampling <- function(population, ncluster1 = 35, 
                     ncluster2 = 110, cluster_observ = 200){
  
  #A function to randomly allocate the number of clusters being sampled in each stratum or each strata.
  #Limited to not less than 2.
  random_cluster <- function(n, m){
    cluster <- rnorm(n, m/n, 0)
    if (abs(sum(cluster)) < 0.01){
      cluster <- cluster + 1
    }
    cluster <- round(cluster / sum(cluster) * m)
    deviation <- m - sum(cluster)
    for (. in seq_len(abs(deviation))) {
      cluster[i] <- cluster[i <- sample(n, 1)] + sign(deviation)
    }
    while (any(cluster < 2)) {
      sm2 <- cluster < 2
      gr2  <- cluster > 2
      cluster[sm2][i] <- cluster[sm2][i <- 
                                        sample(sum(sm2), 1)] + 1
      cluster[gr2][i]  <- cluster[gr2][i <- 
                                         sample(sum(gr2), 1)] - 1
    }
    cluster
  }
  
  #Filtering the simulated population to avoid insufficent clusters or individuals in the sampling stages.
  cluster2count <- population %>% group_by(cluster2) %>% 
    summarise(n = n()) %>% 
    filter(n > 200)
  popu<- population %>% 
    filter(cluster2 %in% cluster2count$cluster2)
  cluster1count <- popu %>% 
    group_by(cluster1) %>% 
    summarise(n = length(unique(cluster2))) %>% 
    filter(n > (ncluster2 / ncluster1))
  
  popu <- popu %>% 
    filter(cluster1 %in% cluster1count$cluster1) %>%
    mutate(strata = as.factor(strata),
           cluster1 = as.factor(cluster1),
           cluster2 = as.factor(cluster2)) %>%
    arrange(strata)
  
  clus_size1 <- random_cluster(length(unique(popu$strata)), ncluster1)
  clus_size2 <- random_cluster(ncluster1, ncluster2)
  
  ind <- mstage(popu, stage = list("stratified", "cluster", "cluster"),
                      varnames = list("strata", "cluster1", "cluster2"),
                      size = list(size1 = table(popu$strata), size2 = clus_size1, 
                                  size3 = clus_size2),
                      method = list("", "srswor", "srswor", "srswor"))
  sample_2nd <- getdata(popu, ind)[[3]] %>% mutate(Prob_2nd = Prob)
  ind_3rd <- sampling::strata(sample_2nd, stratanames = "cluster2", 
                              size = rep(cluster_observ, ncluster2), method = "srswor")
  
  sample <- getdata(sample_2nd, ind_3rd) %>% mutate(pw = 1 / (Prob_2nd * Prob))
  return (sample)
}

