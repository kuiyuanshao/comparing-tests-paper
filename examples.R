
####NHANES

#Code and Data retrieved from https://github.com/tslumley/regression-paper
nhanes <- read.csv("nhanes.csv")
nhanesdesign <- svydesign(id = ~ SDMVPSU,strat = ~ SDMVSTRA,weights = ~ fouryearwt,
                          nest = TRUE, data = subset(nhanes, !is.na(WTDRD1)))
nhanesdesign <- update(nhanesdesign, sodium = DR1TSODI/1000, potassium = DR1TPOTA/1000)
nhanesdesign <- update(nhanesdesign, namol = sodium/23, kmol = potassium/39)

nhanesdesign <- update(nhanesdesign, ish = (BPXSAR>140) & (BPXDAR<90))
nhanesdesign <- update(nhanesdesign, age1 = pmin(RIDAGEYR,50)/10)

mod <- svyglm(DR1TVARA ~ BMXBMI + factor(RIDRETH1)*factor(RIAGENDR), design = nhanesdesign, family = poisson)

waldF <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "Wald") 
scoreF <- svyscoretest(mod, drop.terms= ~factor(RIDRETH1):factor(RIAGENDR), method = "pseudoscore")
lrtF <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "LRT") 
wwaldF <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "WorkingWald") 
wscoreF <- svyscoretest(mod, drop.terms= ~factor(RIDRETH1):factor(RIAGENDR), method = "working")

waldC <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "Wald", df = Inf) 
scoreC <- svyscoretest(mod, drop.terms= ~factor(RIDRETH1):factor(RIAGENDR), ddf = Inf, method = "pseudoscore")
lrtC <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "LRT", df = Inf) 
wwaldC <- regTermTest(mod, ~factor(RIDRETH1):factor(RIAGENDR), method = "WorkingWald", df = Inf) 
wscoreC <- svyscoretest(mod, drop.terms= ~factor(RIDRETH1):factor(RIAGENDR), ddf = Inf, method = "working")

round(c(waldF$p, scoreF[4], lrtF$p, wwaldF$p, wscoreF[3]), 4)
round(c(waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3]), 4)


##### Wilms' tumor
data(nwtco)
set.seed(234)
nwtco$incc2<-as.logical(with(nwtco, ifelse(rel | instit == 2, 1, rbinom(nrow(nwtco), 1, .1))))
dcchs <-twophase(id = list(~ seqno, ~ seqno), strata = list(NULL, ~interaction(rel, instit)),
                data = nwtco, subset = ~ incc2)

mod <- svyglm(rel ~ age + factor(stage)*factor(histol), design = dcchs, family = quasibinomial)
waldC <- regTermTest(mod, ~factor(stage):factor(histol), method = "Wald", df = Inf) 
scoreC <- svyscoretest(mod, drop.terms= ~factor(stage):factor(histol), ddf = Inf, method = "pseudoscore")
lrtC <- regTermTest(mod, ~factor(stage):factor(histol), method = "LRT", df = Inf) 
wwaldC <- regTermTest(mod, ~factor(stage):factor(histol), method = "WorkingWald", df = Inf) 
wscoreC <- svyscoretest(mod, drop.terms= ~factor(stage):factor(histol), ddf = Inf, method = "working")

round(c(waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3]), 4)

gccs8 <- calibrate(dcchs, phase = 2, formula = ~interaction(rel, stage, instit))
mod <- svyglm(rel ~ age + factor(stage)*factor(histol), design = gccs8, family = quasibinomial)
waldC <- regTermTest(mod, ~factor(stage):factor(histol), method = "Wald", df = Inf) 
scoreC <- svyscoretest(mod, drop.terms= ~factor(stage):factor(histol), ddf = Inf, method = "pseudoscore")
lrtC <- regTermTest(mod, ~factor(stage):factor(histol), method = "LRT", df = Inf) 
wwaldC <- regTermTest(mod, ~factor(stage):factor(histol), method = "WorkingWald", df = Inf) 
wscoreC <- svyscoretest(mod, drop.terms= ~factor(stage):factor(histol), ddf = Inf, method = "working")

round(c(waldC$p, scoreC[4], lrtC$p, wwaldC$p, wscoreC[3]), 4)


