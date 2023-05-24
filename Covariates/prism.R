library(sf)
library(prism)
library(raster)
library(dplyr)
library(reshape2)
library(velox)

devtools::install_github("hunzikp/velox")

options(stringsAsFactors = F)
fires15 <- st_read("shapefiles/FireHexes15km")
fires25 <- st_read("shapefiles/FireHexes25km")
ct15 <- st_read("shapefiles/HexControls15km")
ct25 <- st_read("shapefiles/HexControls25km")

keep <- c("HEXID")
shapes15 <- rbind(fires15[,keep], ct15[,keep])
shapes25 <- rbind(fires25[,keep], ct25[,keep])

lst <- "prism/stacks/list_tmean_tmax_tmin_ppt"
prism15 <- data.frame()
prism25 <- data.frame()

for (year in 2000:2018) {
  prism_list <- readRDS(paste0(lst, year, ".rds"))
  names(prism_list) <- c("tmean", "tmax", "tmin", "ppt")
  annual15 <- data.frame()
  annual25 <- data.frame()
  for (var in c("tmean", "tmax", "ppt")) {
    start <- Sys.time()
    print(paste(year, var, start))
    r <- prism_list[[var]]
    vlr <- velox(r)
    vals15 <- vlr$extract(st_transform(shapes15, crs(r)), fun = base::mean)
    vals25 <- vlr$extract(st_transform(shapes25, crs(r)), fun = base::mean)
    formatData <- function(vals, shapes) {
      colnames(vals) <- names(r)
      vals <- (data.frame(vals) %>% mutate(HEXID = shapes$HEXID) %>% tidyr::gather(date, v, -HEXID)
               %>% mutate(date = readr::parse_number(gsub(".*4kmM3_", "", date)),
                          date = as.Date(paste0(date, "01"), "%Y%m%d"))
               )
      names(vals)[names(vals) == "v"] <- var
      return(vals)
    }
    vals15 <- formatData(vals15, shapes15)
    vals25 <- formatData(vals25, shapes25)
    if (var == "tmean") { # first loop
      annual15 <- vals15
      annual25 <- vals25
    } else {
      annual15 <- merge(annual15, vals15, by = c("HEXID", "date"))
      annual25 <- merge(annual25, vals25, by = c("HEXID", "date"))
    }
    print(Sys.time() - start)
  }
  prism15 <- rbind(prism15, annual15)
  prism25 <- rbind(prism25, annual25)
}
  
write.csv(prism15, "Prism15km.csv")
write.csv(prism25, "Prism25km.csv")

