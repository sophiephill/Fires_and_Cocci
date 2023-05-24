
library(sf)
library(raster)
library(rgdal)
library(gdalUtils)
library(tidyverse)
setwd("~/Desktop/CocciWildfires/covariates/")
options(stringsAsFactors = F)

## Download tifs
# for (url in readLines("elevationTifs.txt")) {
#   name <- str_match(url, "(n\\w+)\\.tif")[,2]
#   download.file(url, destfile = paste0("elevation_tifs/", name))
# }

## Get coordinates of center of fire hex
occpoints <- st_read("../hexagons/FireHexes20km/") %>% st_transform(., 3488) %>% st_centroid() %>% st_transform(., 4326) 

## For each point, check which raster it overlaps with
individualTifs <- paste0("elevation_tifs/", list.files("elevation_tifs/")) 
individualTifs <- lapply(individualTifs, raster)

elevs <- data.frame()
for (i in 1:nrow(occpoints)) {
  elev <- c()
  for (r in individualTifs) {
    e <- raster::extract(r, occpoints[i,], method = "simple")
    if (!is.na(e)) {
      elev <- c(elev, e) 
    }
  }
  summ <- list("hexid" = occpoints$HEXID[i], NumMatches = length(elev),
               "Elev" = mean(elev), "Range" = max(elev) - min(elev))
  elevs <- rbind(elevs, summ)
}

# write.csv(elevs, "ElevationFireOccpoints.csv", row.names = F)

## Examine effect by elevation
setwd("~/Desktop/CocciWildfires/gsynth/")
elevs <- read.csv("~/Desktop/CocciWildfires/covariates/ElevationFireOccpoints.csv")
atts <- (read.csv("data/AttFullPeriod.csv") #%>% filter(period == "0-36")
         %>% filter(period %in% c("0-12","13-24","25-36","0-36"))
         %>% mutate(Significant = !((CI.lower < 0) & (CI.upper > 0)))
         %>% merge(., elevs %>% select(hexid, Elev), by = "hexid"))

(ggplot(atts) + geom_point(aes(x = Elev, y = ATT, color = Significant))
  + geom_errorbar(aes(x = Elev, ymin=CI.lower, ymax=CI.upper, color = Significant)) 
  + facet_wrap(~ period)
  + xlab("Elevation at fire source"))

covs <- read.csv("../Synth_analysis/MonthlyCovData/MonthlyCovariatesWPlacebos20km.csv")
mesmaElev <- (covs %>% filter(grepl("FR", id), month < 0) 
              %>% group_by(hexid, season) %>% summarise(MeanMesma = mean(mesma), Precip = mean(ppt))
              %>% merge(., elevs %>% dplyr::select(hexid, Elev), by = "hexid"))
ggplot(mesmaElev%>% filter(hexid %in% atts$hexid)) + geom_point(aes(x = Elev, y = MeanMesma)) + facet_wrap(~ season, scales = "free")

(mesmaElev %>% filter(MeanMesma > .7) %>% filter(hexid %in% atts$hexid)
  %>% merge(., thexes %>% data.frame() %>% select(HEXID, county), by.x = "hexid", by.y = "HEXID")
  %>% group_by(hexid) %>% summarise(county = head(county, 1), elev = head(Elev, 1), mesma = head(MeanMesma, 1))
)

sapply(c("f","sp","su","w"), function(s) cor(mesmaElev[mesmaElev$season==s, "MeanMesma"], mesmaElev[mesmaElev$season==s, "Elev"]))


