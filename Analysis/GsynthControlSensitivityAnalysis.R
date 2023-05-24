

options(stringsAsFactors = F)

prepost20 <- read.csv("covariates/FireCases_20km.csv")
prepost20$total <- prepost20$PreCaseCount + prepost20$PostCaseCount

prepost25 <- read.csv("covariates/FireCases_25km.csv")
prepost25$total <- prepost25$PreCaseCount + prepost25$PostCaseCount
prepost15 <- read.csv("covariates/FireCases_15km.csv")
prepost15$total <- prepost15$PreCaseCount + prepost15$PostCaseCount

response <- c("CaseCount")
timeVaryingCovs <- c("NDaysAbove5ms", "tmax", "tmean", "ppt", "mesma", "population")
fixedCovs <- c("clay", "sand", 'silt', "ph", "impervious")

###############################################################################################
##### gsynth ##################################################################################
###############################################################################################

#### Run gsynth #################################################################################
library(gsynth)

attsOmittedControl <- data.frame()

#allids <- prepost20[prepost20$total > 0, "fire_id"]
pools <- read.csv("../gsynth/data/WhichHexesToPool.csv")
keep <- pools[pools$pool == "keep", "ids"]
kept19 <- thexes20@data[thexes20$HEXID %in% keep, "fire_id"]

for(id in kept19) {
  print(id)
  df <- getDf(id, thexes20, chexes20, wind20, pop15, soil20, prism20, 
              impervious20, fcases20, ctcases20, mesma20, controlFires20)
  df <- df[,!(names(df) %in% fixedCovs)]
  firedate <- getFireDate(id)
  df$intvar <- ifelse(grepl("FR", df$HEXID) & (df$date >= firedate), 1, 0)
  hexid <- df$HEXID[grepl("FR", df$HEXID)][1]
  hexid <- gsub("000", "", hexid)
  begin <- Sys.time()
  for(control in unique(df$HEXID[grepl("CT", df$HEXID)])) {
    dfT <- df[df$HEXID != control,]
    g <- hush(gsynth(formula = CaseCount ~ intvar + NDaysAbove5ms + tmax + tmean + ppt + mesma + population, 
                     data = dfT, 
                     Y = "CaseCount", 
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
                     se = T, nboots = 5000, alpha = 0.05,
                     inference = "parametric", # use parametric when there's <40 treated units
                     cov.ar = 1, 
                     parallel = TRUE, cores = 2, tol = 0.001, 
                     seed = 4343434, min.T0 = 5, normalize = FALSE))
    
    ndates <- length(unique(dfT$date))
    maxMonth <- ndates - 36 - 1
    months <- -36:maxMonth
    
    ## Monthly att
    att <- data.frame(g$est.att)
    att$month <- months
    att$hexid <- hexid
    
    for (p in list(c(0, maxMonth), c(0, 12), c(13, 24), c(25, maxMonth) ) ) {
      ## Att over different post periods period
      sub <- att[(att$month >= max(p[1], 1)) & (att$month <= p[2]),]
      pooledSe <- base::sqrt(sum(sub$S.E. ** 2))
      effect <- sum(sub$ATT)
      periodName <- paste0(p[1], "-", p[2])
      t <- data.frame("hexid" = hexid, "TE" = effect, "seTE" = pooledSe, "period" = periodName)
      t$OmittedControl <- control
      t$hexid <- hexid
      attsOmittedControl <- rbind(attsOmittedControl, t)
    }
  }
  end <- Sys.time()
  print(end - begin)
}


# write.csv(attsOmittedControl, "../gsynth/OmittedControlAtts.csv")

### Pool leaving out 2014, 2015 controls
library(dplyr)
setwd("~/Desktop/CocciWildfires/gsynth/")
keep <- read.csv("data/WhichHexesToPool.csv") %>% pull(ids)
covs <- read.csv("../Synth_analysis/MonthlyCovData/MonthlyCovariatesWPlacebos20km.csv")
covs <- covs[grepl("FR", covs$id), ]
precounts <- aggregate(CaseCount ~ hexid, covs[covs$month < 0,], base::sum)
over10counts <- precounts[precounts$CaseCount >= 10, "hexid"]

source("RmPostFireOverlaps.R") # Some fires overlapped with another in the post period

library(meta)
omit <- read.csv("data/OmittedControlAtts.csv")
atts <- read.csv("data/EffectsAndSes.csv")
atts$period <- gsub("30", "36", atts$period)
atts$period <- gsub("29", "36", atts$period)
omit$period <- gsub("29", "36", omit$period)

effectRange <- function(hexid, p = "0-36") {
  formatOutput <- function(m, label, onlyAll) { # input meta output
    all <- data.frame(hexid = label, ATT = m$TE.fixed, CI.lower = m$lower.fixed, CI.upper = m$upper.fixed, period = p)
    if (!onlyAll) {
      individuals <- data.frame(hexid = sub$hexid, ATT = sub$TE, CI.lower = m$lower, CI.upper = m$upper, period = p)
      all <- rbind(individuals, all)
    }
    return(all)
  }
  og <- metagen(TE = TE, seTE = seTE, data = fixProblems(atts[atts$period == p,], p),
                studlab= paste(hexid), comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
  erange <- formatOutput(og, "ALL", T) %>% mutate(OmittedControl = "NONE")
  
  o <- omit[(omit$hexid == hexid) & (omit$period == p),]
  a <- atts[(atts$hexid != hexid) & (atts$period == p),]
  
  for (control in unique(o$OmittedControl)) {
    sub <- a %>% rbind(., o %>% filter(OmittedControl == control) %>% select(-X, -OmittedControl))
    sub <- fixProblems(sub, p)
    allFires <- metagen(TE = TE, seTE = seTE, data=sub, studlab= paste(hexid), 
                        comb.fixed = TRUE, comb.random = FALSE, prediction=F, sm="MD")
    over10 <- metagen(TE = TE, seTE = seTE,data=sub[sub$hexid %in% over10counts,],
                      studlab= paste(hexid), comb.fixed = TRUE,
                      comb.random = FALSE, prediction=F, sm="MD")
    under10 <- metagen(TE = TE, seTE = seTE,data=sub[!(sub$hexid %in% over10counts),],
                       studlab= paste(hexid), comb.fixed = TRUE,
                       comb.random = FALSE, prediction=F, sm="MD")
    
    tmp <- rbind(formatOutput(allFires, "ALL", onlyAll = T), formatOutput(over10, "ALL Over 10", onlyAll = T))
    tmp <- rbind(tmp, formatOutput(under10, "ALL Under 10", onlyAll = T))
    tmp <- tmp %>% mutate("OmittedControl" = control)
    erange <- rbind(erange, tmp)
  }
  erange <- erange %>% mutate(insig = (CI.lower < 0) & (CI.upper > 0))
  return(erange)
}


eff.ranges <- data.frame()
for (p in c("0-36", "0-12","13-24","25-36")) {
  for (h in keep) {
    r <- effectRange(h, p)
    tmp <- r %>% filter(hexid == "ALL") %>% mutate(fire = h)
    eff.ranges <- rbind(eff.ranges, tmp)
  }
}
write.csv(eff.ranges, "data/EffectRanges.csv", row.names = F)


