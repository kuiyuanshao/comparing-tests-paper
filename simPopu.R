pacman::p_load(dplyr, mvtnorm, scales, future.apply)
#Simulation of Person ID, Strata, and Clusters
sim_idsc <- function(npopu, nstrata, nclsuter1, ncluster2){
  num_strata <- 1:nstrata
  num_cluster1 <- split(1:ncluster1, 
                        rep(1:nstrata,
                            as.integer(rnorm(nstrata, 
                                             mean = ncluster1/nstrata))))
  
  num_cluster2 <- split(1:ncluster2,
                        rep(1:ncluster1,
                            as.integer(rnorm(ncluster1,
                                             mean = ncluster2/ncluster1))))
  num_id <- split(1:npopu,
                  rep(1:ncluster2,
                      as.integer(rnorm(ncluster2,
                                       mean = npopu/ncluster2))))
  
  data <- data.frame(id = 1:npopu)
  length_cluster2 <- as.vector(unlist(lapply(num_cluster2, length)))
  length_cluster1 <- as.vector(unlist(lapply(num_cluster1, length)))
  length_id <- as.vector(unlist(lapply(num_id, length)))
  data$info_cluster1 <- rep(rep(1:ncluster1, length_cluster2), length_id)
  data$info_cluster2 <- rep(1:ncluster2, length_id)
  data$info_strata <- rep(rep(rep(1:nstrata, length_cluster1), 
                              length_cluster2), length_id)
  return (data)
}

#Single iteration function prepared for parallel computation
sim_loop <- function(k, data, prob, mean_matrix, vcov_list, tb, scale, lambda_mean){
  ind <- which(data$info_cluster2 == k)
  i <- unique(data[ind, ]$info_strata)
  j <- unique(data[ind, ]$info_cluster1)
  
  prob_param <- rnorm(100, mean = (10 * i + j / 10 + k / 500) / scale, sd = 0.1)
  probc <- prob + sample(prob_param, 12)
  probc <- probc / sum(probc)
  tbc <- tb
  tbc$prob <- prob
  tbc$sexcat_ind <- 1:nrow(tbc)
  sexcat_ind <- apply(rmultinom(length(ind), size = 1,
                                prob = tbc$prob), 2, which.max)
  data$sexcat_ind[ind] <- sexcat_ind
  data$sex[ind] <- tbc$sex[data$sexcat_ind[ind]]
  data$category[ind] <- tbc$category[data$sexcat_ind[ind]]
  
  curr_scale <- 10 * i + j / 10 + k / 500
  #Each mean of the variable is added some random effets
  #Therefore the mean changes in each stratum, cluster.
  curr_sd <- abs(rnorm(1, mean = i, sd = 10)) + 1
  salary_param <- round(rnorm(100, mean = curr_scale * 10 / scale, 
                              sd = sample(curr_sd, 1)), 2)
  age_param <- as.integer(rnorm(100, mean = curr_scale * 10 / scale, 
                                sd = sample(curr_sd, 1)))
  hour_param <- as.integer(rnorm(100, mean = curr_scale * 15 / scale, 
                                 sd = sample(curr_sd, 1)))
  al_param <- as.integer(rnorm(100, mean = curr_scale * 10 / scale, 
                               sd = sample(curr_sd, 1)))
  sa_param <- rnorm(100, mean = curr_scale * 2 / scale, 
                    sd = sample(curr_sd, 1) / 100)
  ot_param <- abs(rnorm(100, mean = curr_scale / scale, 
                        sd = sample(curr_sd, 1) / 100))
  
  salary_muc <- mean_matrix[, 1] + sample(salary_param, 12)
  age_muc <- mean_matrix[, 2] + sample(age_param, 12)
  hour_muc <- mean_matrix[, 3] + sample(hour_param, 12)
  annualleave_muc <- mean_matrix[, 4] + sample(al_param, 12)
  superannuation_muc <- mean_matrix[, 5] + sample(sa_param, 12)
  overtime_muc <- mean_matrix[, 6] + sample(ot_param, 12)
  for (z in unique(sexcat_ind)){
    sec_ind <- which(data$sexcat_ind[ind] == z)
    
    current_vcov <- vcov_list[[z]] + 
      matrix(rnorm(36, 
                   mean = curr_scale / (scale * 50), 
                   sd = sample(curr_sd, 1) / 100), 
             ncol = 6, nrow = 6)
    current_vcov[lower.tri(current_vcov)] <- t(current_vcov)[lower.tri(current_vcov)]
    current_means <- c(salary_muc[z], 
                       age_muc[z], 
                       hour_muc[z], 
                       annualleave_muc[z], 
                       superannuation_muc[z],
                       overtime_muc[z])
    col_ind <- which(colSums(is.na(data[ind, ])) > 0)
    data[ind, ][sec_ind, col_ind] <- abs(rmvnorm(length(sec_ind), 
                                                 mean = current_means, 
                                                 sigma = current_vcov, 
                                                 method = "svd"))

    
    data$overtime[ind][sec_ind] <- ifelse(data$overtime[ind][sec_ind] > 0.5, 1, 0)
    data$superannuation[ind][sec_ind] <- round(data$superannuation[ind][sec_ind] / (i * 50), 4) + 
      ((10 - i) * 0.002) + 0.002 * i
    data$salary[ind][sec_ind] <- round(data$salary[ind][sec_ind], 2)
    data$age[ind][sec_ind] <- round(data$age[ind][sec_ind], 2)
    data$working_hours[ind][sec_ind] <- round(data$working_hours[ind][sec_ind])
    data$days_annual_leave[ind][sec_ind] <- round(data$days_annual_leave[ind][sec_ind])
  }
  #Outcomes
  outcome_lambda <- as.integer(abs(rnorm(1, sd = curr_scale))) + lambda_mean
  outcome_prob <- as.integer(abs(rnorm(2, sd = curr_scale))) + 1
  data$poissY[ind] <- rpois(length(ind), 
                            lambda = outcome_lambda)
  outcome_prob <- round(outcome_prob / sum(outcome_prob), 3)
  data$binomY[ind] <- abs(rbinom(length(ind), 
                                 size = 1, 
                                 prob = outcome_prob) - 
                            sample(0:1, length(ind), replace = T))
  return (data[ind, ])
}

#Age, Salary per Hour, Working Hours per Week, Whether Overtime, Superannuation, Annual Leave per Year
sim_covariates <- function(data, prob, mean_matrix, vcov_list, scale, lambda_mean){
  sex = c("M", "F", "Not to Say")
  category = c("A", "B", "C", "D") 
  tb <- expand.grid(sex, category)
  colnames(tb) <- c("sex", "category")
  
  data$overtime <- data$superannuation <- data$days_annual_leave <- 
    data$working_hours <- data$age <- data$salary <- 
    data$category <- data$sex <- data$sexcat_ind <- NA
  data$poissY <- data$binomY <- 1
  
  data$sex <- factor(NA, levels = sex)
  data$category <- factor(NA, levels = category)
  
  plan(multisession, workers = 5)
  data_full <- future_lapply(1:length(unique(data$info_cluster2)),
                             sim_loop, 
                             data, prob, 
                             mean_matrix, 
                             vcov_list, 
                             tb, scale, 
                             lambda_mean,
                             future.seed = T)
  
  data_full <- as.data.frame(bind_rows(data_full))
  data_full$sexcat_ind <- NULL
  data_full$age[data_full$age > 75] <- 75
  data_full$age[data_full$age < 16] <- 16
  data_full$days_annual_leave[data_full$days_annual_leave < 21] <- 21
  data_full$days_annual_leave[data_full$days_annual_leave > 50] <- 50
  data_full$salary[data_full$salary < 8] <- 8
  data_full$working_hours[data_full$working_hours > 75] <- 75
  return (data_full)
}


