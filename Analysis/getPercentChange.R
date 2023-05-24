
### Here I bootstrap to obtain uncertainty estimates for percent change
cases <- read.csv("data/CasesPooledByPeriod.csv") %>% filter(grepl("ALL", hexid))
head(cases)

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
