
library(sf)
library(rgdal)
library(tidyverse)
options(stringsAsFactors = F)

#### Fit the ITS model
library(foreign)
library(tsModel)
library(lmtest)
library(Epi)
library(splines)
library(vcd)

keep <- read.csv("data/WhichHexesToPool.csv") %>% pull(ids)
covs <- read.csv("~/Desktop/CocciWildfires/Synth_analysis/MonthlyCovData/MonthlyCovariatesWPlacebos20km.csv")
covs <- covs %>% filter(hexid %in% keep, grepl("FR", id))
covs <- covs %>% mutate(truemonth = lubridate::month(as.Date(date, origin = "1970-01-01")))
head(as.Date(covs$date, origin = "1970-01-01"))

plotEffects <- function(var) {
  postfires <- data.frame()
  for (h in keep) {
    df <- covs %>% filter(hexid == h, month < 13) %>% mutate(postfire = as.integer(month > 0))
    
    mod <- glm(as.formula(paste0(var, "~ postfire + truemonth + harmonic(truemonth, 2, 12)")), family=gaussian, df)
    effects <- round(ci.lin(mod,Exp = T), 3)["postfire",]
    tmp <- data.frame(hexid = df$hexid[1], postfire = effects["exp(Est.)"], ci.lower = effects["2.5%"], 
                      ci.upper = effects["97.5%"], row.names = c())
    postfires <- rbind(postfires, tmp)
  }
  postfires$significant <- postfires$hexid %in% c("FR13", "FR32","FR91","FR92")
  plot <- (ggplot(postfires) 
           + geom_pointrange(aes(x = reorder(hexid, postfire), y = postfire, ymin=ci.lower, ymax=ci.upper, color = significant)) 
           + xlab("") + ylab("exp(postfire)") + geom_abline(intercept = 1, slope = 0))
  print(plot)
  return(postfires)
}

tmp <- plotEffects("tmean")

plotRawVars <- function(var) {
  df <- (covs %>% filter(month %in% 0:12) %>% group_by(hexid) 
         %>% summarise(avg = mean(!!sym(var)), sd = sd(!!sym(var))) 
         %>% mutate(lower = avg - sd, upper = avg + sd) )
  df$significant <- postfires$hexid %in% c("FR13", "FR32","FR91","FR92")
  plot <- (ggplot(df) 
           + geom_pointrange(aes(x = reorder(hexid, avg), y = avg, ymin= lower, ymax=upper, color = significant)) 
           + xlab("") + ylab(paste(var, "average 0-12 months post fire")))
  print(plot)
  return(df)
}

# EFFECTS
ci.lin(mod,Exp=T)["postfire",5:7]

# TREND
exp(coef(mod)["month"]*12)

# We again check the model and autocorrelation functions
res <- residuals(mod,type="deviance")
plot(res,ylim=c(-5,10),pch=19,cex=0.7,col=grey(0.6),main="Residuals over time",
     ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2)
acf(res)
pacf(res)


