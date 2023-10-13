pacman::p_load(sampling, dplyr, hash, stringr)

#### Simulates the population of interest 

currdir <- getwd()

source(paste0(currdir, "/simPopu.R"))
npopu <- 1e7
nstrata <- 10
ncluster1 <- 200
ncluster2 <- 8000


prob <- c(0.0575, 0.0614, 0.0089, 
          0.1406, 0.1023, 0.0128, 
          0.1074, 0.1023, 0.0233, 
          0.1751, 0.1636, 0.0447)

mean_matrix <- 
  cbind(c(50, 45, 53, 33, 31, 35,
          20, 21, 17, 25, 27, 29),
        c(35, 34, 41, 45, 32, 28,
          46, 26, 30, 53, 31, 49),
        c(30, 32, 40, 50, 41, 54, 
          34, 40, 31, 43, 31, 39),
        c(35, 35, 28, 28, 28, 28,
          21, 28, 28, 21, 35, 28),
        c(0.05, 0.06, 0.035, 0.034, 0.04, 0.03,
          0.04, 0.035, 0.045, 0.035, 0.035, 0.041),
        c(0.01, 0.01, 0.05, 0.02, 0.06, 0.02,
          0.10, 0.03, 0.12, 0.03, 0.15, 0.03))

set.seed(100)
vcov_list <- NULL
for (i in 1:12){
  sigma_mat <- diag(ncol = 6, nrow = 6)
  diag(sigma_mat) <- (1 / 5) * mean_matrix[i, ]
  sigma_mat[sigma_mat == 0] <- rep(diag(sigma_mat) / 2, each = 5) + 
    rnorm(30, mean = 0, sd = 0.1)
  sigma_mat[lower.tri(sigma_mat)] <- t(sigma_mat)[lower.tri(sigma_mat)]
  vcov_list[[i]] <- sigma_mat
  
}
scale <- 10 * nstrata + ncluster1 / 10 + ncluster2 / 500

initial_data <- sim_idsc(npopu, nstrata, ncluster1, ncluster2)

set.seed(12345)
system.time({popu1 <- sim_covariates(initial_data, prob, 
                                     mean_matrix, vcov_list, 
                                     scale, lambda_mean = 40)})

dir.create(file.path(currdir, "popu1"), showWarnings = FALSE)
save(popu1, file = paste0(currdir, "/popu1/popu1.RData"))

set.seed(54321)
system.time({popu2 <- sim_covariates(initial_data, prob, 
                                     mean_matrix, vcov_list, 
                                     scale, lambda_mean = 0)})

dir.create(file.path(currdir, "popu2"), showWarnings = FALSE)
save(popu2, file = paste0(currdir, "/popu2/popu2.RData"))



