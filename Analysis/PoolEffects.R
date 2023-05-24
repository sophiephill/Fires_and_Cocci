
library(meta)
library(dplyr)
library(stringr)
setwd("~/Desktop/CocciWildfires/gsynth/")
setwd("~/Desktop/CocciWildfires/gsynth/")

pools <- read.csv("data/WhichHexesToPool.csv")
keep <- pools[pools$pool == "keep", "ids"]
covs <- read.csv("../Synth_analysis/MonthlyCovData/MonthlyCovariatesWPlacebos20km.csv")
covs <- covs[grepl("FR", covs$id), ]
precounts <- aggregate(CaseCount ~ hexid, covs[covs$month < 0,], base::sum)
over10counts <- precounts[precounts$CaseCount >= 10, "hexid"]

#### Pool ATT estimates ############################################################################
# Pools effects by period
atts <- read.csv("25km/EffectsAndSes.csv")
# atts <- read.csv("25km/EffectsAndSes.csv")
atts$period <- gsub("30", "36", atts$period)
atts$period <- gsub("29", "36", atts$period)
atts <- atts[atts$hexid %in% keep, ]

metas <- data.frame()
for (p in unique(atts$period)) {
  sub <- atts[atts$period == p, ]
  sub <- fixProblems(sub, p, p == "0-36") # only inclue partials for full period
  partials <- sub %>% filter(period != p)
  sub <- sub %>% filter(period == p)
  allFires <- metagen(TE = TE, seTE = seTE, data=sub, studlab= paste(hexid), 
               comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
  over10 <- metagen(TE = TE, seTE = seTE,data=sub[sub$hexid %in% over10counts,],
                studlab= paste(hexid), comb.fixed = TRUE,
                comb.random = FALSE, prediction=F, sm="MD")
  under10 <- metagen(TE = TE, seTE = seTE,data=sub[!(sub$hexid %in% over10counts),],
                studlab= paste(hexid), comb.fixed = TRUE,
                comb.random = FALSE, prediction=F, sm="MD")
  formatOutput <- function(m, label, onlyAll) { # input meta output
    all <- data.frame(hexid = label, ATT = m$TE.fixed, CI.lower = m$lower.fixed, CI.upper = m$upper.fixed, period = p)
    if (!onlyAll) {
      individuals <- data.frame(hexid = sub$hexid, ATT = sub$TE, CI.lower = m$lower, CI.upper = m$upper, period = p)
      all <- rbind(individuals, all)
    }
    return(all)
  }
  tmp <- rbind(formatOutput(allFires, "ALL", onlyAll = F), formatOutput(over10, "ALL Over 10", onlyAll = T))
  tmp <- rbind(tmp, formatOutput(under10, "ALL Under 10", onlyAll = T))
  if (nrow(partials) > 0 ) {
    mp <- metagen(TE = TE, seTE = seTE, data=partials, studlab= paste(hexid), 
                        comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
    
    tmp <- rbind(tmp, data.frame(hexid = partials$hexid, ATT = partials$TE, CI.lower = mp$lower, CI.upper = mp$upper, period = partials$period))
  }
  metas <- rbind(metas, tmp)
}

metas

# write.csv(metas, "15km/AttFullPeriod.csv", row.names = F)

#metas <- read.csv("~/Desktop/CocciWildfires/gsynth/data/AttFullPeriod.csv")
metas$insig <- (metas$CI.lower < 0)&(metas$CI.upper > 0)
metas %>% mutate(name = NameFire(hexid)) %>% filter(period == "0-36")
sig <- metas[!metas$insig,]
sig %>% arrange(hexid)


#### Cases by period #########################
# Pools observed/predicted cases by period
cases <- read.csv("data/ObsSynthCaseCounts.csv")
atts <- read.csv("data/AttMonthly.csv")
atts <- merge(atts, cases, by = c("hexid", "month"))
atts <- atts[atts$hexid %in% keep, c("hexid", "month", "observed", "predicted", "S.E.")]

periodCases <- data.frame()
## Sum cases over many months: also sum the variances
## This is so the percent change takes into account the weighting used to pool the ATTs across all
for (hexid in unique(atts$hexid)) {
  a <- atts[atts$hexid == hexid,]
  for (p in list(c(0, 12), c(13, 24), c(25, 36), c(0, 36))) {
    sub <- a[(a$month >= p[1]) & (a$month <= p[2]),]
    periodName <- paste0(p[1], "-", p[2])
    pooledSe <- base::sqrt(sum(sub$S.E. ** 2))
    obs <- sum(sub$observed)
    pred <- sum(sub$predicted)
    t <- data.frame("hexid" = hexid, "observed" = obs, "predicted" = pred, "seTE" = pooledSe, "period" = periodName)
    periodCases <- rbind(periodCases, t)
  }
}
problems
periodCases

metas <- data.frame()
for (p in unique(periodCases$period)) { 
  sub <- periodCases[periodCases$period == p, ]
  sub <- fixProblems(sub, p)
  
  formatOutput <- function(m, data, allLabel, measure, onlyAll) { # input meta output
    all <- data.frame(hexid = allLabel, measure = m$TE.fixed, CI.lower = m$lower.fixed, CI.upper = m$upper.fixed, period = p)
    if (!onlyAll) {
      individuals <- data.frame(hexid = data$hexid, measure = data$TE, CI.lower = m$lower, CI.upper = m$upper, period = p)
      all <- rbind(individuals, all)
    } 
    names(all)[names(all) == "measure"] <- measure
    return(all)
  }
  
  getPooledEstimate <- function(measure, data, allLabel, onlyAll) {
    names(data)[names(data) == measure] <- "TE"
    m <- metagen(TE = TE, seTE = seTE, data=data, studlab= paste(hexid),
                 comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
    w <- weights(m, comb.fixed = T) 
    w <- w %>% mutate(hexid = row.names(w))
    return( list(df = formatOutput(m, data, allLabel, measure, onlyAll), weights = w) )
  }
  
  allFires <- getPooledEstimate("predicted", sub, "ALL", F)
  over10 <- getPooledEstimate("predicted", sub[sub$hexid %in% over10counts,], "ALL Over 10", T)
  under10 <- getPooledEstimate("predicted", sub[!(sub$hexid %in% over10counts),], "ALL Under 10", T)
  pred <- allFires$df %>% rbind(over10$df) %>% rbind(under10$df)
  
  observedIndividuals <- (getPooledEstimate("observed", sub, "ALL", F)$df %>% filter(grepl("FR", hexid)) 
                          %>% mutate(obsSE = NA) ) # get observed case data, no pools
  poolObsFromWeights <- function(w, label) {
    obs <- (w %>% merge(., observedIndividuals, by = "hexid") %>% mutate(contrib = observed * w.fixed) 
             %>% mutate(contrib = contrib / sum(w$w.fixed)) %>% pull(contrib) )
    return(data.frame(hexid = label, observed = sum(obs), CI.lower = NA, CI.upper = NA, period = p, obsSE = sd(obs)))
  }
  
  observedAll <- (poolObsFromWeights(allFires$weights, "ALL") 
          %>% rbind(poolObsFromWeights(over10$weights, "ALL Over 10")) 
          %>% rbind(poolObsFromWeights(under10$weights, "ALL Under 10"))
          %>% rbind(observedIndividuals))
  
  tmp <- observedAll %>% select(hexid, period, observed, obsSE) %>% merge(pred, by = c("hexid","period"))
  metas <- rbind(metas, tmp)
}

observedSe <- cases %>% group_by(hexid) %>% summarise(obsSE2 = sd(observed))
metas <- metas %>% merge(., observedSe, by = "hexid", all.x = T)
metas <- metas %>% mutate(obsSE = sapply(1:length(obsSE), function(i) ifelse(is.na(obsSE[i]), obsSE2[i], obsSE[i]))) %>% select(-obsSE2)
metas <- metas %>% mutate(predSE = CI.upper - CI.lower) %>% mutate(predSE = predSE/(1.96*2))
# delta refers to the delta method ie how the Census calculates se for percent change
metas <- metas %>% mutate(percentChange = (observed - predicted) / predicted,
                          percChangeSeDelta = abs(observed/predicted) * sqrt((obsSE**2 / observed**2) + (predSE**2 / predicted**2)),
                          DeltaCI.lower = percentChange - 1.96 * percChangeSeDelta, DeltaCI.upper = percentChange + 1.96 * percChangeSeDelta)
head(metas)
 # write.csv(metas, "data/CasesPooledByPeriod.csv", row.names = F)

#### Monthly Pooled Effects #####
# This gets the effect estimates for Month X, the code above handles ranges like 0-12
# This is to make an interval around the case counts plot
# Here I get a CI for the mean predicted counterfactual for each month
# For the effect, I just average the predictions
# For the SE, I use the SE of the effect since var(observed - predicted) = var(predicted)
cases <- read.csv("data/ObsSynthCaseCounts.csv")
atts <- read.csv("data/AttMonthly.csv")
atts <- merge(atts, cases, by = c("hexid", "month"))
atts <- atts[atts$hexid %in% keep, c("hexid", "month", "observed", "predicted", "S.E.")]


metas <- data.frame()
for (month in -36:36) {
  sub <- atts[atts$month == month, ]
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
  
  allFires <- getPooledEstimate("predicted", sub, "ALL", F)
  over10 <- getPooledEstimate("predicted", sub[sub$hexid %in% over10counts,], "ALL Over 10", T)
  under10 <- getPooledEstimate("predicted", sub[!(sub$hexid %in% over10counts),], "ALL Under 10", T)

  pred <- allFires %>% rbind(over10) %>% rbind(under10)
  
  allFires <- getPooledEstimate("observed", sub, "ALL", F)
  over10 <- getPooledEstimate("observed", sub[sub$hexid %in% over10counts,], "ALL Over 10", T)
  under10 <- getPooledEstimate("observed", sub[!(sub$hexid %in% over10counts),], "ALL Under 10", T)
  
  obs <- allFires %>% rbind(over10) %>% rbind(under10)
  
  tmp <- obs %>% select(hexid, month, observed) %>% merge(pred, by = c("hexid","month"))
  metas <- rbind(metas, tmp)
}


# write.csv(metas, "data/MonthlyAverageEffect.csv", row.names = F)

