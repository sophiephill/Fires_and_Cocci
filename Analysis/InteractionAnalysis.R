library(rgdal)
library(dplyr)

options(stringsAsFactors = F)

covariates <- read.csv("MonthlyCovariatesWPlacebos20km.csv")
thexes <- readOGR("../hexagons/FireHexes20km/")

fires <- thexes@data[,c("HEXID","igntn_d")]
fires$igmonth<- as.integer(stringr::str_match(fires$igntn_d, "-(\\d+)-")[,2])
fires$igseason <- ifelse(fires$igmonth %in% c(12,1,2), "w",
                         ifelse(fires$igmonth%in% c(3,4,5), "sp",
                                ifelse(fires$igmonth%in% c(6,7,8), "su", "f")))
covariates$fall <- ifelse(covariates$season == "f", 1, 0)
covariates$summer <- ifelse(covariates$season == "su", 1, 0)
covariates$spring <- ifelse(covariates$season == "sp", 1, 0)
covariates$InFireYear <-  as.integer(covariates['year']==0)
covariates <- merge(covariates, fires[,c("HEXID", "igseason")], by.x = "hexid",
                    by.y = "HEXID", all.x = T)
keep <- read.csv("data/WhichHexesToPool.csv") %>% pull(ids)
covariates <- covariates[covariates$hexid %in% keep,]
covariates$InFireSeason <- ifelse((covariates$year==0)&(covariates$season==covariates$igseason),1,0)
covariates$InFireMonth <- ifelse((covariates$month==0),1,0)
covariates$IsFireHex <-  as.integer(grepl("FR", covariates$id))

yearorseason <- "InFireYear"
covariates$interaction <- covariates[[yearorseason]]*covariates$IsFireHex

fm <- as.formula(paste0("NDaysAbove5ms~fall+spring+summer+interaction+IsFireHex+",yearorseason ))
model <- lm(fm, covariates)
#poismodel <- glm(NDaysAbove5ms~fall+spring+summer+InFireSeason+IsFireHex+interaction, covariates, family = quasipoisson)
summary(model)
confint(model, c("IsFireHex", "interaction"),level = 0.95)


