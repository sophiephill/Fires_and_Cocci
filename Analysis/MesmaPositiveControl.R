
library(dplyr)
options(stringsAsFactors = F)

prepost20 <- read.csv("covariates/FireCases_20km.csv")
prepost20$total <- prepost20$PreCaseCount + prepost20$PostCaseCount

allids <- prepost20[prepost20$total > 0, "fire_id"]
pools <- read.csv("../gsynth/data/WhichHexesToPool.csv")
keep <- pools[pools$pool == "keep", "ids"]

fixedCovs <- c("clay", "sand", 'silt', "ph", "impervious")

#### MESMA analysis ###############################################################################
library(gsynth)
AllGapsEvenPlacebos <- data.frame()
controlWeights <- data.frame()
pvals <- data.frame()
avgAtts <- data.frame()
counts <- data.frame()

keep19 <- sapply(keep, hexidToFireid, thexes = thexes20)

for(id in keep19) {
  print(id)
  df <- getDf(id, thexes20, chexes20, wind20, pop20, soil20, prism20, 
              impervious20, fcases20, ctcases20, mesma20, controlFires20)
  df <- df[,!(names(df) %in% fixedCovs)]
  firedate <- getFireDate(id)
  df$intvar <- ifelse(grepl("FR", df$HEXID) & (df$date >= firedate), 1, 0)
  hexid <- df$HEXID[grepl("FR", df$HEXID)][1]
  
  g <- hush(gsynth(formula = mesma ~ intvar + NDaysAbove5ms + tmax + tmean + ppt + population, 
                   data = df, 
                   Y = "mesma", 
                   D = "intvar", 
                   X = timeVaryingCovs,
                   na.rm = FALSE,
                   index = c("HEXID", "date"), 
                   weight = NULL, 
                   force = "two-way", 
                   r = 1, CV = T, 
                   ### Options if you use a different algorithm for regularization: 
                   # lambda = NULL, nlambda = 10, k = 5,
                   EM = FALSE, # Makes more precisely estimated coefs but takes longer to run
                   estimator = "ife",
                   se = T, nboots = 10000, alpha = 0.05,
                   inference = "parametric", # use parametric when there's <40 treated units
                   cov.ar = 1, 
                   parallel = TRUE, cores = 2, tol = 0.001, 
                   seed = NULL, min.T0 = 5, normalize = FALSE))
  
  ndates <- length(unique(df$date))
  maxMonth <- ndates - 36 - 1
  months <- -36:maxMonth
  
  ## Effects
  placebogaps <- data.frame(g$res.co)
  placebogaps$month <- months
  placebogaps <- reshape2::melt(placebogaps, "month", variable.name = "id", value.name = "gap")
  
  tmp <- data.frame("month" = months,
                    "id" = hexid,
                    "gap" = g$eff[,1],
                    row.names = c())
  tmp <- rbind(tmp, placebogaps)
  hexid <- gsub("000", "", hexid)
  tmp$hexid <- hexid
  AllGapsEvenPlacebos <- rbind(AllGapsEvenPlacebos, tmp)
  
  ## Monthly att
  att <- data.frame(g$est.att)
  att$month <- months
  att$hexid <- hexid
  pvals <- rbind(pvals, att)
  
  ## Observed vs counterfactual counts to calculate percent change
  c <- data.frame("hexid" = hexid, "month" = months, 
                  "observed" = g$Y.tr[,1], "predicted" = g$Y.ct[,1])
  counts <- rbind(counts, c)
  
  ## Att over different post periods period
  for (p in list(c(0, 12), c(13, 24), c(25, maxMonth), c(0, maxMonth))) {
    sub <- att[(att$month >= p[1]) & (att$month <= p[2]),]
    pooledSe <- base::sqrt(sum(sub$S.E. ** 2))
    effect <- sum(sub$ATT)
    periodName <- paste0(p[1], "-", p[2])
    t <- data.frame("hexid" = hexid, "TE" = effect, "seTE" = pooledSe, 
                    "period" = periodName)
    avgAtts <- rbind(avgAtts, t)
  }
}

setwd("~/Desktop/CocciWildfires/gsynth/")
# write.csv(counts, "MesmaAnalysis/MesmaRawMeasures.csv")
# write.csv(avgAtts, "MesmaAnalysis/MesmaEffectsAndSes.csv")
# write.csv(AllGapsEvenPlacebos, "MesmaAnalysis/MesmaGapsEvenPlacebos.csv")
# write.csv(pvals, "MesmaAnalysis/AttMonthly.csv")

##### Pool ######################################################################################################
setwd("~/Desktop/CocciWildfires/gsynth/")
source("RmPostFireOverlaps.R") # Some fires overlapped with another in the post period

library(meta)
# Here I get a CI for the mean predicted counterfactual for each month
# For the effect, I just average the predictions
# For the SE, I use the SE of the effect since var(observed - predicted) = var(predicted)
mesma <- read.csv("MesmaAnalysis/MesmaRawMeasures.csv")
mesma <- mesma[,2:ncol(mesma)]
atts <- read.csv("MesmaAnalysis/AttMonthly.csv")
mesma <- merge(mesma, atts[,c("hexid", "month", "S.E.")], by = c("hexid", "month"))
mesma <- mesma[mesma$hexid %in% keep, ]

metas <- data.frame()
for (month in -36:36) {
  sub <- mesma[mesma$month == month, ]
  sub <- fixProblems(sub, month)
  
  formatOutput <- function(m, data, allLabel, measure, onlyAll) { # input meta output
    all <- data.frame(hexid = allLabel, measure = m$TE.fixed, CI.lower = m$lower.fixed, CI.upper = m$upper.fixed, month = month)
    if (!onlyAll) {
      individuals <- data.frame(hexid = data$hexid, measure = data$TE, CI.lower = m$lower, CI.upper = m$upper, month = month)
      all <- rbind(individuals, all)
    } 
    names(all)[names(all) == "measure"] <- measure
    return(all)
  }
  
  getPooledEstimate <- function(measure, data, allLabel, onlyAll) {
    names(data)[names(data) == measure] <- "TE"
    m <- metagen(TE = TE, seTE = S.E., data=data, studlab= paste(hexid), 
                 comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
    return(formatOutput(m, data, allLabel, measure, onlyAll))
  }
  
  pred <- getPooledEstimate("predicted", sub, "ALL", F)
  obs <- getPooledEstimate("observed", sub, "ALL", F)
  tmp <- obs %>% select(hexid, month, observed) %>% merge(pred, by = c("hexid","month"))
  metas <- rbind(metas, tmp)
}
head(metas)
# write.csv(metas, "MesmaAnalysis/MonthlyPredictedInterval.csv", row.names=FALSE)


### Quick reg: 
# is decrease in month 11-12 associated with sig fires?
mesma <- read.csv("MesmaAnalysis/MonthlyPredictedInterval.csv")
sigs <- c("FR92", "FR91","FR32","FR13")
m <- (mesma %>% filter(hexid != "ALL", month >= 0) 
      %>% mutate(sig = hexid %in%sigs, decrease = abs(observed - predicted) < 0.01)
      %>% mutate(mesmachange = observed-predicted))
logreg <- glm(decrease ~ sig, data = m[m$month%in%11:12,], family = binomial)
data.frame(coef = logreg$coefficients, exp.coef = exp(logreg$coefficients))
summary(logreg) # more likely to NOT have a decrease, but also much more likely to have minor change
m[(m$hexid %in% sigs)&(m$month %in% 11:12),]

# Is change in vegetative landcover associated with changes in cases?
cases <- read.csv("data/AttMonthly.csv") %>% filter(hexid != "ALL")
d <- merge(cases[,c("hexid","month","ATT")], m %>% select(hexid,month,sig,mesmachange), by = c("hexid", "month"))

dates <- data.frame()
for (hexid in unique(cases$hexid)) {
  fireid <- hexidToFireid(hexid, thexes20)
  tmp <- getCaseData(fireid,fcases20, ctcases20, thexes20) 
  tmp <- tmp[tmp$HEXID==hexid,]
  tmp$month <- -36:(length(unique(tmp$date)) - 37)
  tmp$hexid <- hexid
  dates <- rbind(dates, tmp[,c("month", "date","hexid")])
}

seasons <- c("w", "w", "sp",'sp',"sp","su","su","su","f","f","f","w")
d <-( merge(d, dates, by = c("month", "hexid")) %>% mutate(date = as.Date(date, "%Y-%m-%d"), truemonth = lubridate::month(date))
      %>% mutate(season = seasons[truemonth], sigmonth = month %in% 11:12) %>% filter(month > 0))

linreg <- lm(ATT ~ mesmachange + season + sig, data = d)
confint(linreg)

# tiny but significant negative association (increase in case att ~ -0.2% change in mesma)

mesma %>% filter(month %in% 10:12, hexid %in% sigs) %>% mutate(ATT = observed - predicted) %>% arrange(hexid) %>% head(20)


