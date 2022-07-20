#######################################################################################################################
################################# P-Values Simulation For Five Hypothesis Testing Approaches ##########################
#######################################################################################################################

library(survey)
library(sampling)
library(future.apply)

#This function gets the p-value for the an individual sample.
pvalues <- function(population, ncluster1 = 35, ncluster2 = 100, 
                    cluster_observ = 200, regression, term = 1){
  sample <- sampling(population, ncluster1, ncluster2, cluster_observ)
  strat_design <- svydesign(id = ~ cluster1 + cluster2 + personID, 
                            strata = ~ strata, 
                            data = sample, weights = ~ pw)
  if (regression == "poisson"){
    mod <- svyglm(out1 ~ V1 + V2 + V3 + V4 + out2, strat_design, 
                  family = poisson())
  }else if(regression == "linear"){
    mod <- svyglm(out1 ~ V1 + V2 + V3 + V4 + out2, strat_design, 
                  family = gaussian())
  }else if(regression == "binomial"){
    mod <- svyglm(out2 ~ V1 + V2 + V3 + V4 + out1, strat_design, 
                  family = quasibinomial())
  }
  
  if (term == 1){
    waldF <- regTermTest(mod, ~V1, method = "Wald")
    scoreF <- svyscoretest(mod, ~V1, method = "pseudoscore")
    lrtF <- regTermTest(mod, ~V1, method = "LRT")
    wwaldF <- regTermTest(mod, ~V1, method = "WorkingWald")
    wscoreF <- svyscoretest(mod, ~V1, method = "working")
    
    waldC <- regTermTest(mod, ~V1, method = "Wald", df = Inf)
    scoreC <- svyscoretest(mod, ~V1, method = "pseudoscore", ddf = Inf)
    lrtC <- regTermTest(mod, ~V1, method = "LRT", df = Inf)
    wwaldC <- regTermTest(mod, ~V1, method = "WorkingWald", df = Inf)
    wscoreC <- svyscoretest(mod, ~V1, method = "working", ddf = Inf)
    
    return(c(waldF$p, scoreF, lrtF$p, wwaldF$p, wscoreF, 
             waldC$p, scoreC, lrtC$p, wwaldC$p, wscoreC))
  }else{
    waldF <- regTermTest(mod, ~V1 + V2 + V3 + V4, 
                         method = "Wald")
    scoreF <- svyscoretest(mod, ~V1 + V2 + V3 + V4, 
                            method = "pseudoscore")
    lrtF <- regTermTest(mod, ~V1 + V2 + V3 + V4, 
                        method = "LRT")
    wwaldF <- regTermTest(mod, ~V1 + V2 + V3 + V4, 
                          method = "WorkingWald")
    wscoreF <- svyscoretest(mod, ~V1 + V2 + V3 + V4, 
                             method = "working")
    
    waldC <- regTermTest(mod, ~V1 + V2 + V3 + V4, 
                         method = "Wald", df = Inf)
    scoreC <- svyscoretest(mod, ~V1 + V2 + V3 + V4, 
                            method = "pseudoscore", ddf = Inf)
    lrtC <- regTermTest(mod, ~V1 + V2 + V3 + V4,
                        method = "LRT", df = Inf)
    wwaldC <- regTermTest(mod, ~V1 + V2 + V3 + V4, 
                          method = "WorkingWald", df = Inf)
    wscoreC <- svyscoretest(mod, ~V1 + V2 + V3 + V4, 
                             method = "working", ddf = Inf)
    
    return(c(waldF$p, scoreF[4], lrtF$p, wwaldF$p, wscoreF[3], 
             waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3]))
  }
}

#Use of parallel computing to speed up the simulation
plan(multisession, works = 4)
pmatrix_generator <- function(population, iteration = 200, ncluster1 = 2, 
                              ncluster2 = 2, cluster_observ = 200, regression, term = 1){
  
  result <- future_replicate(n = iteration, pvalues(population, ncluster1, 
                                   ncluster2, cluster_observ, 
                                   regression, term))
  return (result)
}

#transferring the p-values to uniformly distributed.
exp_p <- function(pmatrix){
  exp_p <- pmatrix
  for (i in 1:dim(exp_p)[1]){
    exp_p[i, ] <- (rank(pmatrix[i, ], ties.method="first")) / (length(pmatrix[i, ]) + 1)
  }
  exp_p
}