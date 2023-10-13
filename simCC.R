pacman::p_load(MASS, survey, mvtnorm, stringr)

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

sim_CC <- function(n, prob, mean_matrix, vcov_list){
  currdir <- getwd()
  dir.create(file.path(currdir, "sample"), showWarnings = FALSE)
  for (i in 1:1000){
    set.seed(i * 10)
    sex = c("M", "F", "Not to Say")
    category = c("A", "B", "C", "D") 
    tb <- expand.grid(sex, category)
    colnames(tb) <- c("sex", "category")
    
    data <- data.frame(id = 1:n)
    
    tbc <- tb
    tbc$prob <- prob
    tbc$sexcat_ind <- 1:nrow(tbc)
    sexcat_ind <- apply(rmultinom(n, size = 1,
                                  prob = tbc$prob), 2, which.max)
    data$sexcat_ind <- sexcat_ind
    data$sex <- tbc$sex[data$sexcat_ind]
    data$category <- tbc$category[data$sexcat_ind]
    
    data$x6 <- data$x5 <- data$x4 <- 
      data$x3 <- data$x2 <- data$x1 <- NA
    
    for (z in unique(sexcat_ind)){
      sec_ind <- which(data$sexcat_ind == z)
      cuurent_mean <- mean_matrix[z, ]
      current_vcov <- vcov_list[[z]]
      x <- abs(rmvnorm(length(sec_ind), 
                     mean = cuurent_mean, 
                     sigma = current_vcov, 
                     method = "svd"))
      data[sec_ind, 5:10] <- x
    }
    data$sexcat_ind <- NULL
    data$x2[data$x2 > 75] <- 75
    data$x2[data$x2 < 16] <- 16
    data$x4[data$x4 < 21] <- 21
    data$x4[data$x4 > 50] <- 50
    data$x1[data$x1 < 8] <- 8
    data$x3[data$x3 > 75] <- 75
    
    y1 <- rexp(n, rate = 1)
    y1 <- as.numeric(y1 < 0.1)
    xstrats <- interaction(data$x1 > quantile(data$x1)[4], 
                           data$x2 > quantile(data$x2)[2])
    
    tf <- sum(xstrats == "TRUE.FALSE")
    
    set.seed(i * 10)
    samplephase2 <- stratsample(xstrats, 
                                  c(FALSE.FALSE = 100,
                                    TRUE.FALSE = tf,
                                    FALSE.TRUE = 100,
                                    TRUE.TRUE = 100))
    samplephase2 <- unique(c(samplephase2, which(y1 == 1)))
      
    id <- (1:n) %in% samplephase2
      
    info_strata <- as.factor(ifelse(y1 == 1, "CASE", as.character(xstrats)))
      
    data$info_strata <- info_strata
    data$sample.phase2 <- id
    data$y1 <- y1
    data$pw <- 0
    weights <- table(info_strata) / table(info_strata[data$sample.phase2])
    data$pw[data$sample.phase2 == T & data$info_strata == "CASE"] <- weights[1]
    data$pw[data$sample.phase2 == T & data$info_strata == "FALSE.FALSE"] <- weights[2]
    data$pw[data$sample.phase2 == T & data$info_strata == "FALSE.TRUE"] <- weights[3]
    data$pw[data$sample.phase2 == T & data$info_strata == "TRUE.FALSE"] <- weights[4]
    data$pw[data$sample.phase2 == T & data$info_strata == "TRUE.TRUE"] <- weights[5]
      
    samp <- data
    save(samp, 
         file = paste0(currdir, "/sample/Sample_", 
                       str_pad(i, nchar(1000), pad = 0), ".RData"), 
                       compress = 'xz')
  }
}

sim_CC(5000, prob, mean_matrix, vcov_list)
