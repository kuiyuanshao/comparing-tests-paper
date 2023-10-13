pacman::p_load(tidyverse, survey, stringr, progress)

pValues <- function(samp){
  samp <- samp[samp$sample.phase2 == T, ]
  multidesign <- svydesign(ids = ~ 1,
                           strata = ~ info_strata,
                           data = samp, weights = ~ pw)
  bmod <- svyglm(y1 ~ category + sex + 
                   x1 + x2 + 
                   x3 + 
                   x4,
                 multidesign, 
                 family = quasibinomial())
  
  multi_b <- multiTest(bmod, term = 5) 
  
  
  result <- matrix(multi_b, nrow = 1)
  colnames(result) <- c("Wald", "Quasi-score", "Likelihood ratio", 
                        "Working Wald", "Working score")
  
  result <- as.data.frame(result) %>% 
    pivot_longer(cols = 1:5, names_to = "Test", values_to = "p-value")
  return (result)
}

multiTest <- function(mod, term){
  waldC <- regTermTest(mod, ~category + sex, method = "Wald", df = Inf) 
  scoreC <- svyscoretest(mod, drop.terms= ~category + sex, ddf = Inf, method = "pseudoscore")
  lrtC <- regTermTest(mod, ~category + sex, method = "LRT", df = Inf) 
  wwaldC <- regTermTest(mod, ~category + sex, method = "WorkingWald", df = Inf) 
  wscoreC <- svyscoretest(mod, drop.terms= ~category + sex, ddf = Inf, method = "working")
    
  return (c(waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3]))
}

#transferring the p-values to uniform.
exp_p <- function(pmatrix){
  exp_mat <- pmatrix
  exp_mat$`p-value` <- as.numeric(exp_mat$`p-value`)
  test <- c("Wald", "Quasi-score", "Likelihood ratio", 
            "Working Wald", "Working score")
  for (i in 1:5){
    ind <- which(pmatrix$Test == test[i])
    currp <- as.numeric(pmatrix$`p-value`[ind])
    exp_mat[ind, 2] <- (rank(currp, ties.method="first")) / (length(currp) + 1)  
  }
  exp_mat
}

currdir <- getwd()

cc <- progress_bar$new(
  format = "Running :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 1000, width = 60)
files <- paste0(currdir, '/sample/Sample_', str_pad(1:1000, nchar(1000), pad = 0), ".RData")
cc_mat <- NULL
for (i in 1:1000){
  load(files[i])
  result <- pValues(samp)
  cc_mat <- rbind(cc_mat, result)
  cc$tick(tokens = list(what = "pvalues   "))
  Sys.sleep(2 / 100)
}

cc_mat_unif <- exp_p(cc_mat)
cc_mat_log <- cc_mat
cc_mat_log$expected <- -log10(cc_mat_unif$`p-value`)
cc_mat_log$observed <- -log10(as.numeric(cc_mat$`p-value`))
cc_mat_log$`p-value` <- NULL
cc_mat_log$observed[cc_mat_log$observed > 2.5] <- 2.5
cc_mat_log$Test <- factor(cc_mat_log$Test, levels = c("Wald", "Quasi-score", "Likelihood ratio", 
                                                          "Working Wald", "Working score"))
save(cc_mat_unif, cc_mat, cc_mat_log, file = "cc_pvalues.RData")




