
setwd("~/Desktop/CocciWildfires/ExtractGaps/")
source("../synthetic_functions.R")
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
AllGapsEvenPlacebos <- data.frame()
controlWeights <- data.frame()
pvals <- data.frame()
avgAtts <- data.frame()
rateRatios <- data.frame()
counts <- data.frame()
highestWeightedControl <- data.frame()

allids <- prepost20[prepost20$total > 0, "fire_id"]
#allids <- prepost25[prepost25$PreCaseCount >= 1, "fire_id"]

# Exclude poorly fit placebos so they don't bias the bootstrap
# Run gsynth twice: first to identify badfits, then remove them
preMse <- function(df) {
  df <- df[grepl("CT", df$id),]
  df <- df[df$month < 0, ]
  df <- aggregate(gap ~ id, df, function(v) mean(v**2))
  return(df)
}
badFits <- data.frame()
badfits <- read.csv("../gsynth/PoorlyFitPlacebosMse5.csv")

## For 15km, some regions too small (few cases) so I didn't extract cov data
haveDataOn <- hexidToFireid(soil15[grepl("FR", soil15$HEXID), 'HEXID'], thexes15)

# If you only need data on the fires we kept
pools <- read.csv("../gsynth/WhichHexesToPool.csv")
keep <- pools[pools$pool == "keep", "ids"]
kept19 <- sapply(keep, hexidToFireid, thexes = thexes20)

for(id in kept19) {
  ### If 15km, some regions were too small so I didn't extract cov data ####
  #if (!(id %in% haveDataOn)) {
   # next
  #}
  ###########################################################################
  print(id)
  df <- getDf(id, thexes20, chexes20, wind20, pop15, soil20, prism20, 
              impervious20, fcases20, ctcases20, mesma20, controlFires20)
  df <- df[,!(names(df) %in% fixedCovs)]
  firedate <- getFireDate(id)
  df$intvar <- ifelse(grepl("FR", df$HEXID) & (df$date >= firedate), 1, 0)
  hexid <- df$HEXID[grepl("FR", df$HEXID)][1]
  
  ### Filter out placebos that were poorly fit so they don't influence the inference ####
  # exclude <- badfits[(badfits$hexid == gsub("000", "", hexid)) & (badfits$gap > 100), "id"]
  # df <- df[!(df$HEXID %in% exclude), ]
  ##################################################################################
  g <- hush(gsynth(formula = CaseCount ~ intvar + NDaysAbove5ms + tmax + tmean + ppt + mesma + population, 
                   data = df, 
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
                   se = T, nboots = 8000, alpha = 0.05,
                   inference = "parametric", # use parametric when there's <40 treated units
                   cov.ar = 1, 
                   parallel = TRUE, cores = 2, tol = 0.001, 
                   seed = 4343434, min.T0 = 5, normalize = FALSE))
  
  ndates <- length(unique(df$date))
  maxMonth <- ndates - 36 - 1
  months <- -36:maxMonth
  
  ## Plain Old Effects
  placebogaps <- data.frame(g$res.co)
  placebogaps$month <- months
  placebogaps2 <- reshape2::melt(placebogaps, "month", variable.name = "id", value.name = "gap")
  
  tmp <- data.frame("month" = months,
                    "id" = hexid,
                    "gap" = g$eff[,1],
                    row.names = c())
  tmp <- rbind(tmp, placebogaps2)
  hexid <- gsub("000", "", hexid)
  tmp$hexid <- hexid
  AllGapsEvenPlacebos <- rbind(AllGapsEvenPlacebos, tmp)
  
  ## Which controls are poorly fit and should be removed to not bias inference?
  mses <-preMse(tmp)
  mses <- mses[mses$gap > 5, ]
  if (nrow(mses) > 0) {
    mses$hexid <- hexid
    badFits <- rbind(badFits, mses)
  }
  
  ##### Weights per control 
  # ct <- g$wgt.implied
  # ct <- data.frame("hexid" = hexid, "id" = row.names(ct), 
  #                  "weight" = ct[,1], row.names = c())
  # controlWeights <- rbind(controlWeights, ct)
  
  
  ## Monthly att + SE + CI 
  att <- data.frame(g$est.att)
  att$month <- months
  att$hexid <- hexid
  pvals <- rbind(pvals, att)
  
  ## Observed vs counterfactual raw counts to calculate percent change
  count <- data.frame("hexid" = hexid, "month" = months, 
                     "observed" = g$Y.tr[,1], "predicted" = g$Y.ct[,1])
  counts <- rbind(counts, count)
  
  ## Find ATT over different periods
  boots <- data.frame(g$att.boot, month = months) # nmonths by 200 sized dataframe of bootstrapped ATTs
  
  for (p in list(c(0, 12), c(13, 24), c(25, maxMonth), c(0, maxMonth))) {
    ## Att over different post periods period
    sub <- att[(att$month >= max(p[1], 1)) & (att$month <= p[2]),]
    pooledSe <- base::sqrt(sum(sub$S.E. ** 2))
    effect <- sum(sub$ATT)
    periodName <- paste0(p[1], "-", p[2])
    t <- data.frame("hexid" = hexid, "TE" = effect, "seTE" = pooledSe, 
                    "period" = periodName)
    avgAtts <- rbind(avgAtts, t)
  }
}

setwd("~/Desktop/CocciWildfires/gsynth/")
# write.csv(counts, "ObsSynthCaseCounts.csv")
# write.csv(avgAtts, "EffectsAndSes.csv")
# write.csv(pvals, "AttMonthly.csv")



# write.csv(AllGapsEvenPlacebos, "AllGapsEvenPlacebos.csv")
# write.csv(badFits, "PoorlyFitPlacebosMse5.csv")
# write.csv(controlWeights, "ControlWeights.csv")
