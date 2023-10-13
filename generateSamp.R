pacman::p_load(sampling, dplyr, hash, stringr)

#A function to randomly allocate the number of clusters being sampled in each stratum or each info_strata.
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

# A function for randomly selecting the clusters without replacement
phaseSamp <- function(hash, num){
  res <- NULL
  prob <- NULL
  for (i in 1:length(hash)){
    res <- c(res, sample(hash[[as.character(i)]], num[i]))
    prob <- c(prob, num[i] / length(hash[[as.character(i)]]))
  }
  return (list(res, prob))
}

generateSamp <- function(data, nstratum, ninfo_cluster1, ninfo_cluster2, nobserv){
  #Filtering the simulated population to avoid insufficent clusters or individuals in the sampling stages.
  info_cluster2count <- data %>% group_by(info_cluster2) %>% 
    summarise(n = n()) %>% 
    filter(n > 200)
  popu <- data %>% 
    filter(info_cluster2 %in% info_cluster2count$info_cluster2)
  info_cluster1count <- popu %>% 
    group_by(info_cluster1) %>% 
    summarise(n = length(unique(info_cluster2))) %>% 
    filter(n > (ninfo_cluster2 / ninfo_cluster1))
  
  popu <- popu %>% 
    filter(info_cluster1 %in% info_cluster1count$info_cluster1) %>%
    arrange(info_strata)
  
  #Randomly allocate the number of clusters to be sampled in different stages
  clus_size1 <- random_cluster(length(unique(popu$info_strata)), ninfo_cluster1)
  clus_size2 <- random_cluster(ninfo_cluster1, ninfo_cluster2)
  
  #Grab some population information and join into the data
  allostrata <- popu %>% 
    group_by(info_strata) %>%
    summarise(n_strata = n())
  
  alloclus1 <- popu %>% 
    group_by(info_strata, info_cluster1) %>% 
    summarise(n_clus1 = n())
  
  alloclus2 <- popu %>% 
    group_by(info_strata, info_cluster1, info_cluster2) %>% 
    summarise(n_clus2 = n())
  
  popu <- popu %>% 
    left_join(allostrata) %>%
    left_join(alloclus1) %>% 
    left_join(alloclus2)
  
  #Randomly select allocated number of cluster1 in each strata
  strata_hash <- hash()
  for (i in 1:nstratum){
    strata_hash[i] <- as.numeric(alloclus1$info_cluster1[alloclus1$info_strata == i])
  }
  #Get the result and probability corresponds to
  selstage1 <- phaseSamp(strata_hash, clus_size1)
  selstage1_prob <- data.frame(info_strata = 1:10, 
                               prob_stage1 = selstage1[[2]])
  selstage1 <- selstage1[[1]]

  #Randomly select allocated number of cluster2 in each selected cluster1
  clus1_hash <- hash()
  for (i in 1:length(selstage1)){
    clus1_hash[i] <- as.numeric(alloclus2$info_cluster2[alloclus2$info_cluster1 == selstage1[i]])
  }
  #Get the result and probability corresponds to
  selstage2 <- phaseSamp(clus1_hash, clus_size2)
  selstage2_prob <- data.frame(info_cluster1 = selstage1, 
                               prob_stage2 = selstage2[[2]])
  selstage2 <- selstage2[[1]]
  #Randomly select nobserv in each selected cluster2
  clus2_hash <- hash()
  for (i in 1:length(selstage2)){
    clus2_hash[i] <- as.numeric(popu$id[popu$info_cluster2 == selstage2[i]])
  }
  #Get the result and probability corresponds to
  selstage3 <- phaseSamp(clus2_hash, rep(nobserv, length(selstage2)))
  selstage3_prob <- data.frame(info_cluster2 = selstage2, 
                               prob_stage3 = selstage3[[2]])
  selstage3 <- selstage3[[1]]
  #Filtering the data and compute the probability to be selected for each individual
  sample <- popu[popu$id %in% selstage3, ] %>%
    left_join(selstage1_prob) %>%
    left_join(selstage2_prob) %>%
    left_join(selstage3_prob) %>%
    mutate(prob = prob_stage1 * prob_stage2 * prob_stage3) %>%
    mutate(pw = 1 / prob)

  return (sample)
}

#generate samples for 25 design degrees of freedom
currdir <- getwd()
dir.create(file.path(currdir, "popu1_samp_25df"), showWarnings = FALSE)
load(paste0(currdir, "/popu1/popu1.RData"))

for (i in 1:1000){
  set.seed(i * 10)
  samp <- generateSamp(popu1, 10, 35, 100, 200)
  save(samp, file = 
         paste0(currdir, "/popu1_samp_25df/Sample_", 
                str_pad(i, nchar(1000), pad = 0), ".RData"), 
       compress = 'xz')
}

#generate samples for 30 design degrees of freedom

#generate samples for 35 design degrees of freedom
currdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dir.create(file.path(currdir, "popu1_samp_35df"), showWarnings = FALSE)
load(paste0(currdir, "/popu1/popu1.RData"))

for (i in 1:1000){
  set.seed(i * 5)
  samp <- generateSamp(popu1, 10, 45, 200, 100)
  save(samp, file = 
         paste0(currdir, "/popu1_samp_35df/Sample_", 
                str_pad(i, nchar(1000), pad = 0), ".RData"), 
       compress = 'xz')
}

# Popu2

#generate samples for 25 design degrees of freedom
currdir <- getwd()
dir.create(file.path(currdir, "popu2_samp_25df"), showWarnings = FALSE)
load(paste0(currdir, "/popu2/popu2.RData"))

for (i in 1:1000){
  set.seed(i * 10)
  samp <- generateSamp(popu2, 10, 35, 100, 200)
  save(samp, file = 
         paste0(currdir, "/popu2_samp_25df/Sample_", 
                str_pad(i, nchar(1000), pad = 0), ".RData"), 
       compress = 'xz')
}

#generate samples for 30 design degrees of freedom

#generate samples for 35 design degrees of freedom
currdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dir.create(file.path(currdir, "popu2_samp_35df"), showWarnings = FALSE)
load(paste0(currdir, "/popu2/popu2.RData"))

for (i in 1:1000){
  set.seed(i * 5)
  samp <- generateSamp(popu2, 10, 45, 200, 100)
  save(samp, file = 
         paste0(currdir, "/popu2_samp_35df/Sample_", 
                str_pad(i, nchar(1000), pad = 0), ".RData"), 
       compress = 'xz')
}
