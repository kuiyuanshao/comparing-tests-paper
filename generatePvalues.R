pacman::p_load(tidyverse, survey, stringr, progress)

pValues <- function(samp){
  multidesign <- svydesign(ids = ~ info_cluster1 + info_cluster2 + id,
                           strata = ~ info_strata,
                           data = samp, weights = ~ pw)
  
  pmod <- svyglm(poissY ~ category + sex + 
                   salary + age + 
                   working_hours + 
                   days_annual_leave,
                 multidesign, 
                 family = poisson())
  
  bmod <- svyglm(binomY ~ category + sex +  
                   salary + age + 
                   working_hours + 
                   days_annual_leave,
                 multidesign, 
                 family = quasibinomial())
  
  single_p <- multiTest(pmod, term = 1) 
  single_b <- multiTest(bmod, term = 1) 
  
  multi_p <- multiTest(pmod, term = 3) 
  multi_b <- multiTest(bmod, term = 3) 
  
  
  result_F <- rbind(single_p[[1]], single_b[[1]], multi_p[[1]], multi_b[[1]])
  result_C <- rbind(single_p[[2]], single_b[[2]], multi_p[[2]], multi_b[[2]])
  result_F <- cbind(result_F, c("Poisson", "Binomial", "Poisson", "Binomial"))
  result_C <- cbind(result_C, c("Poisson", "Binomial", "Poisson", "Binomial"))
  result_F <- cbind(result_F, c("1 term", "1 term", "5 terms", "5 terms"))
  result_C <- cbind(result_C, c("1 term", "1 term", "5 terms", "5 terms"))
  
  result_F <- cbind(result_F, rep("F", 4))
  result_C <- cbind(result_C, rep("C", 4))
  
  result <- rbind(result_F, result_C)
  
  colnames(result) <- c("Wald", "Quasi-score", "Likelihood ratio", 
                        "Working Wald", "Working score",
                        "Outcome Distribution", "q", "Reference Distribution")
  
  result <- as.data.frame(result) %>% 
    pivot_longer(cols = 1:5, names_to = "Test", values_to = "p-value")
  return (result)
}

multiTest <- function(mod, term){
  if (term == 1){
    waldF <- regTermTest(mod, ~salary, method = "Wald") 
    scoreF <- svyscoretest(mod, drop.terms= ~salary, method = "pseudoscore")
    lrtF <- regTermTest(mod, ~salary, method = "LRT") 
    wwaldF <- regTermTest(mod, ~salary, method = "WorkingWald") 
    wscoreF <- svyscoretest(mod, drop.terms= ~salary, method = "working")
    
    waldC <- regTermTest(mod, ~salary, method = "Wald", df = Inf) 
    scoreC <- svyscoretest(mod, drop.terms= ~salary, ddf = Inf, method = "pseudoscore")
    lrtC <- regTermTest(mod, ~salary, method = "LRT", df = Inf) 
    wwaldC <- regTermTest(mod, ~salary, method = "WorkingWald", df = Inf) 
    wscoreC <- svyscoretest(mod, drop.terms= ~salary, ddf = Inf, method = "working")
    
    return (list(c(waldF$p, as.numeric(scoreF[4]), lrtF$p, wwaldF$p, wscoreF[3]), 
                 c(waldC$p, as.numeric(scoreC[4]), lrtC$p, wwaldC$p, wscoreC[3])))
  }else{
    waldF <- regTermTest(mod, ~category + sex, method = "Wald") 
    scoreF <- svyscoretest(mod, drop.terms= ~category + sex, method = "pseudoscore")
    lrtF <- regTermTest(mod, ~category + sex, method = "LRT") 
    wwaldF <- regTermTest(mod, ~category + sex, method = "WorkingWald") 
    wscoreF <- svyscoretest(mod, drop.terms= ~category + sex, method = "working")
    
    waldC <- regTermTest(mod, ~category + sex, method = "Wald", df = Inf) 
    scoreC <- svyscoretest(mod, drop.terms= ~category + sex, ddf = Inf, method = "pseudoscore")
    lrtC <- regTermTest(mod, ~category + sex, method = "LRT", df = Inf) 
    wwaldC <- regTermTest(mod, ~category + sex, method = "WorkingWald", df = Inf) 
    wscoreC <- svyscoretest(mod, drop.terms= ~category + sex, ddf = Inf, method = "working")
    
    return (list(c(waldF$p, scoreF[4], lrtF$p, wwaldF$p, wscoreF[3]), 
                 c(waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3])))
  }
}

#transferring the p-values to uniform.
exp_p <- function(pmatrix){
  exp_mat <- pmatrix
  exp_mat$`p-value` <- as.numeric(exp_mat$`p-value`)
  test <- c("Wald", "Quasi-score", "Likelihood ratio", 
            "Working Wald", "Working score")
  ref <- c("F", "C")
  q <- c("1 term", "5 terms")
  out <- c("Poisson", "Binomial")
  for (i in 1:5){
    for (j in 1:2){
      for (k in 1:2){
        for (z in 1:2){
          ind <- which((pmatrix$Test == test[i]) & 
                         (pmatrix$`Reference Distribution` == ref[j]) &
                         (pmatrix$q == q[k]) & (pmatrix$`Outcome Distribution` == out[z]))
          currp <- as.numeric(pmatrix$`p-value`[ind])
          exp_mat[ind, 5] <- (rank(currp, ties.method="first")) / (length(currp) + 1)  
        }
      }
    }
  }
  exp_mat
}

currdir <- getwd()

#### Population 1

#25 Design degrees of freedom
pb25 <- progress_bar$new(
  format = "Running :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 1000, width = 60)
files <- paste0(currdir, '/popu1_samp_25df/Sample_', str_pad(1:1000, nchar(1000), pad = 0), ".RData")
df25_mat <- NULL
for (i in 1:1000){
  load(files[i])
  result <- pValues(samp)
  df25_mat <- rbind(df25_mat, result)
  pb25$tick(tokens = list(what = "pvalues   "))
  Sys.sleep(2 / 100)
}


#35 Design degrees of freedom
pb35 <- progress_bar$new(
  format = "Running :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 1000, width = 60)
files <- paste0(currdir, '/popu1_samp_35df/Sample_', str_pad(1:1000, nchar(1000), pad = 0), ".RData")
df35_mat <- NULL
for (i in 1:1000){
  load(files[i])
  result <- pValues(samp)
  df35_mat <- rbind(df35_mat, result)
  pb35$tick(tokens = list(what = "pvalues   "))
}


df25_mat_unif <- exp_p(df25_mat)
df35_mat_unif <- exp_p(df35_mat)

df25_mat_log <- df25_mat
df35_mat_log <- df35_mat
df25_mat_log$expected <- -log10(df25_mat_unif$`p-value`)
df35_mat_log$expected <- -log10(df35_mat_unif$`p-value`)
df25_mat_log$observed <- -log10(as.numeric(df25_mat$`p-value`))
df35_mat_log$observed <- -log10(as.numeric(df35_mat$`p-value`))

df25_mat_log$`p-value` <- NULL
df35_mat_log$`p-value` <- NULL

df25_mat_log$observed[df25_mat_log$observed > 2.5] <- 2.5
df35_mat_log$observed[df35_mat_log$observed > 2.5] <- 2.5

df25_mat_log$Test <- factor(df25_mat_log$Test, levels = c("Wald", "Quasi-score", "Likelihood ratio", 
                                                          "Working Wald", "Working score"))
df35_mat_log$Test <- factor(df35_mat_log$Test, levels = c("Wald", "Quasi-score", "Likelihood ratio", 
                                                          "Working Wald", "Working score"))

save(df25_mat_unif, df25_mat, df25_mat_log, file = "popu1_df25_pvalues.RData")
save(df35_mat_unif, df35_mat, df35_mat_log, file = "popu1_df35_pvalues.RData")



#### Population 2

#25 Design degrees of freedom
pb25 <- progress_bar$new(
  format = "Running :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 1000, width = 60)
files <- paste0(currdir, '/popu2_samp_25df/Sample_', str_pad(1:1000, nchar(1000), pad = 0), ".RData")
df25_mat <- NULL
for (i in 1:1000){
  load(files[i])
  result <- pValues(samp)
  df25_mat <- rbind(df25_mat, result)
  pb25$tick(tokens = list(what = "pvalues   "))
  Sys.sleep(2 / 100)
}


#35 Design degrees of freedom
pb35 <- progress_bar$new(
  format = "Running :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 1000, width = 60)
files <- paste0(currdir, '/popu2_samp_35df/Sample_', str_pad(1:1000, nchar(1000), pad = 0), ".RData")
df35_mat <- NULL
for (i in 1:1000){
  load(files[i])
  result <- pValues(samp)
  df35_mat <- rbind(df35_mat, result)
  pb35$tick(tokens = list(what = "pvalues   "))
}


df25_mat_unif <- exp_p(df25_mat)
df35_mat_unif <- exp_p(df35_mat)

df25_mat_log <- df25_mat
df35_mat_log <- df35_mat
df25_mat_log$expected <- -log10(df25_mat_unif$`p-value`)
df35_mat_log$expected <- -log10(df35_mat_unif$`p-value`)
df25_mat_log$observed <- -log10(as.numeric(df25_mat$`p-value`))
df35_mat_log$observed <- -log10(as.numeric(df35_mat$`p-value`))

df25_mat_log$`p-value` <- NULL
df35_mat_log$`p-value` <- NULL

df25_mat_log$observed[df25_mat_log$observed > 2.5] <- 2.5
df35_mat_log$observed[df35_mat_log$observed > 2.5] <- 2.5

df25_mat_log$Test <- factor(df25_mat_log$Test, levels = c("Wald", "Quasi-score", "Likelihood ratio", 
                                                          "Working Wald", "Working score"))
df35_mat_log$Test <- factor(df35_mat_log$Test, levels = c("Wald", "Quasi-score", "Likelihood ratio", 
                                                          "Working Wald", "Working score"))

save(df25_mat_unif, df25_mat, df25_mat_log, file = "popu2_df25_pvalues.RData")
save(df35_mat_unif, df35_mat, df35_mat_log, file = "popu2_df35_pvalues.RData")




