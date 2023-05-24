
library(sf)
library(rgdal)
library(tidyverse)
options(stringsAsFactors = F)
setwd("~/Desktop/CocciWildfires/gsynth/")

#### Make fire boundary shapefile that I will extract mesma values to
# 
# fires <- (st_read("wildfire_data/") 
#             %>% st_transform(., CRS("+init=epsg:3488 +units=m"))
#             %>% mutate(igntn_d = paste(YEAR, STARTMONTH,STARTDAY, sep = "/"), igntn_d = as.Date(igntn_d)))
# 
# keep <- (st_read("hexagons/FireHexes20km/") %>% select(HEXID, fire_id) 
#          %>% merge(.,  read.csv("gsynth/data/WhichHexesToPool.csv"), by.x = "HEXID", by.y = "ids")
#          %>% pull(fire_id))
# keepFires <- fires %>% filter(FIRE_ID %in% keep)
# st_write(keepFires, "~/Desktop/FireBoundaries19/", "FireBoundaries19", driver = "ESRI Shapefile", overwrite_layer = T)

#### Clean the data
npv <- (read.csv("data/NPVFires.csv") %>% mutate(month = lubridate::month(date), year = lubridate::year(date))
        %>% group_by(year, month, fire_id) %>% summarise(npv = base::mean(npv))
        %>% mutate(month = str_pad(month, 2, "left", "0"), date = paste(year, month, "01", sep = "-"))
        %>% as.data.frame()%>% select(fire_id, date, npv))
gv <- (read.csv("data/GVFires.csv") %>% mutate(month = lubridate::month(date), year = lubridate::year(date))
        %>% group_by(year, month, fire_id) %>% summarise(gv = base::mean(gv))
        %>% mutate(month = str_pad(month, 2, "left", "0"), date = paste(year, month, "01", sep = "-"))
        %>% as.data.frame()%>% select(fire_id, date, gv))
mesma <- merge(npv, gv, by = c("fire_id", "date")) %>% mutate(mesma = npv + gv)

fireidKey <- st_read("../hexagons/FireHexes20km/") %>% as.data.frame() %>%  select(HEXID, fire_id)
mesma <- (read.csv("../Synth_analysis/MonthlyCovData/MonthlyCovariatesWPlacebos20km.csv")
        %>% filter(grepl("FR", id))
        %>% mutate(id = gsub("000", "", id), date = as.Date(date, origin = "1970-01-01"))
        %>% select(id, date, month, year)
        %>% merge(., fireidKey, by.x = "id", by.y = "HEXID", all.x = T)
        %>% merge(., mesma, by = c("fire_id", "date")))
mesma <- mesma %>% mutate(postfire = as.integer(month > 0))
head(mesma)

#### Fit the ITS model
library(foreign)
library(tsModel)
library(lmtest)
library(Epi)
library(splines)
library(vcd)

keep <- (st_read("../hexagons/FireHexes20km/") %>% select(HEXID, fire_id)
         %>% merge(.,  read.csv("data/WhichHexesToPool.csv"), by.x = "HEXID", by.y = "ids")
         %>% pull(fire_id))
head(mesma)

postfires <- data.frame()
for (fireid in keep) {
  df <- mesma %>% filter(fire_id == fireid, month < 13)
  mod <- glm(mesma ~ postfire + month + harmonic(month, 2, 12), family=gaussian, df)
  effects <- round(ci.lin(mod,Exp = T), 3)["postfire",]
  tmp <- data.frame(fire_id = fireid, hexid = df$id[1], postfire = effects["exp(Est.)"], ci.lower = effects["2.5%"], 
                    ci.upper = effects["97.5%"], row.names = c())
  postfires <- rbind(postfires, tmp)
 
}

postfires$significant <- ifelse(postfires$hexid %in% c("FR13", "FR32","FR91","FR92"),
                                "Followed by Significant\nChanges", "Insignificant Effect")
postfires$FireName <- NameFire(postfires$hexid)

p <- (ggplot(postfires) + geom_pointrange(aes(y = reorder(FireName, postfire), x = postfire, xmin=ci.lower, xmax=ci.upper, color = significant)) 
  + ylab("") + xlab("IRR") + geom_vline(xintercept = 1) 
  + labs(color='') + theme(axis.text.y = element_text(size=12) ) 
  + theme(legend.text=element_text(size=12)))
ggsave("~/Desktop/CocciWildfires/figures/MesmaIts.jpeg", dpi = 600)

(ggplot(mesma) + geom_line(aes(x = month, y = mesma)) + facet_wrap(~id))

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










