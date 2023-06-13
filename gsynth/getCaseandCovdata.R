
library(rgdal)
library(rgeos)
library(reshape2)
library(zoo)
library(sf)
library(stringr)

##################################################################################
############ Prep data #######################################################
##############################################################################

# Determine eligible controls
# maxControls is an optional argument, include for the robustness check
getEligible <- function(fire_id, chexes, thexes, controlFires, n = 3) {
  # Over 10% of area overlaps with fire hexagon
  firehex <- thexes[thexes$fire_id==fire_id,]
  firehex <- as(firehex, "sf")
  chexesSf <- as(chexes, "sf")
  insex <- sf::st_intersection(firehex, chexesSf)
  intersecting_controls <- insex$HEXID.1
  percents <- sf::st_area(insex)/sf::st_area(chexesSf[chexesSf$HEXID %in% intersecting_controls,])
  overlaps <- intersecting_controls[as.numeric(percents) > .1]
  # Experienced a fire within n years
  igYear <- thexes@data[thexes$fire_id==fire_id, "year"]
  overlaps <- c(overlaps, 
                controlFires[abs(controlFires$year - igYear) <= n, "FID"])
  eligibles <- chexes[!(chexes$HEXID %in% overlaps),]
  eligibles$HEXID <- gsub("ID", "CT", eligibles$HEXID)
  return(eligibles)
}


# Convert onset month to date  (synth needs day month year) 
# OR Convert covariates apr_2012 format to R date
convertDate <- function(s, onset = T) {
  if (onset) {
    onset <- as.integer(regmatches(s, regexpr("\\d+", s)))
    month <- (c(12, 1:11))[(onset %% 12) + 1]
    month <- c(paste0(0, 1:9), 10:12)[month]
    dt <- paste("01", month, 2000+onset%/%12, sep = "/")
    return(as.Date(dt, format = "%d/%m/%Y"))
    return(dt)
  } else {
    months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                "sep", "oct","nov", "dec")
    names(months) <- c(paste0(0, 1:9), 10:12)
    year <- substring(s, 7, 8)
    mon <- sapply(s, function(s) {names(months)[which(months == substr(s, 1, 3))]})
    s <- paste("01", mon, year, sep = "_")
    s <- gsub("_", "/", s)
    return(as.Date(s, "%d/%m/%y"))
  }
}

# Get case data +/- n years after fire
# nbefore: num months before fire
# removes polygons with zero cases (necessary for placeboSynth function)
getCaseData <- function(fireid, fcases, ctcases, thexes, nbefore = 36, nafter = 36, warning = F) {
  endPre <- fcases[fcases$fire_id == fireid, "OnsetMonth"]
  startPost <- endPre
  
  # modified fireid handles case when fires were so close they should've been counted as one: 
  # instead you "cut out" the time between fires ie don't match or measure an effect in this time
  modifiedFireid <- fireidForDates(fireid)
  if (length(modifiedFireid) == 2) {
    endPre <- fcases[fcases$fire_id == modifiedFireid[1], "OnsetMonth"]
    startPost <- fcases[fcases$fire_id == modifiedFireid[2], "OnsetMonth"]
  } else if (fireid == "CA3571711855520140613") { # FR92
    startPost <- startPost + 3 # overlapped with a fire we missed on 8/18/14
  }
  startPre <- max(1, endPre - nbefore) 
  endPost <- min(216, startPost + nafter) 
  # if there were fires in 210/214, we want to end matching 209 and start effect 214
  keepCols <- paste0("OnsetMonth", c(startPre:(endPre - 1), startPost:endPost))

  caseDf <- rbind(fcases[fcases$fire_id == fireid,keepCols], ctcases[,keepCols])
  rownames(caseDf) <- c(fcases[fcases$fire_id == fireid, "HEXID"], ctcases$HEXID)
  caseDf <- reshape2::melt(t(caseDf))
  colnames(caseDf) <- c("date", "HEXID", "CaseCount")
  caseDf$date <- convertDate(caseDf$date,  onset = T)
  return(caseDf)
}

# This is for fires that need to be merged together due to spatial proximity
hexidToFireid <- function(k, thexes = thexes20) {
  return( sapply(k, function(hex) thexes@data[thexes$HEXID == hex, "fire_id"]) )
}

fireidToHexid <- function(k, thexes = thexes20) {
  return( sapply(k, function(id) thexes@data[thexes$fire_id == id, "HEXID"]) )
}

fireidForDates <- function(fireid) {
  problems <- c("FR7","FR12","FR40","FR44")
  hexid <- fireidToHexid(fireid, thexes20)
  # Formatted earlier, later
  if (any(hexid %in% problems[1:2])) {
    return( hexidToFireid(c("FR12","FR7"), thexes20) ) 
  } else if(any(hexid %in% problems[3:4])) {
    return( hexidToFireid(c("FR40","FR44"), thexes20) )
  } else {
    return(fireid)
  }
}


getDf <- function(fireid, thexes, chexes, wind, pop, prism, fcases, ctcases, mesma, 
                  controlFires, maxControls = data.frame()) {
  elControls <- getEligible(fireid, chexes, thexes, controlFires = controlFires, 
                            maxControls = maxControls)
  keepIds <- c(thexes@data[thexes$fire_id==fireid, "HEXID"], elControls$HEXID)
  combd <- merge(wind[wind$HEXID %in% keepIds,2:4], 
                 prism[prism$HEXID %in% keepIds,c("HEXID", "date","tmax", "tmean", "ppt")],
                 by = c("HEXID", "date"))
  combd <- merge(combd, mesma[mesma$HEXID %in% keepIds, c("HEXID", "date","mesma")],
                 by = c("HEXID", "date"),all.x=T)
  # I didn't extract all 19 years of cov data for later fires
  nhexmonths <- nrow(combd[combd$HEXID==thexes@data[thexes$fire_id==fireid, "HEXID"],])
  nmonths <-c(rep(228,(length(unique(combd$HEXID))-1)), nhexmonths)
  cases <- getCaseData(fireid, fcases, ctcases[ctcases$HEXID %in% keepIds,], thexes = thexes)
  cases$HEXID <- as.character(cases$HEXID)
  if (any(grepl("a", combd$date))) { # I used to encode date as apr_2011 etc
    combd$date <- convertDate(combd$date, onset = F)
  } else {
    combd$date <- as.Date(combd$date)
  }
  combd <- combd[order(combd$HEXID, combd$date), ]
  combd$year <- substring(combd$date,1,4)
  pop <- pop[order(pop$HEXID, pop$year),]
  pop <- pop[pop$HEXID %in% keepIds,]
  combd <- merge(combd, pop, by = c("HEXID","year"), all.x = T)
  combd <- combd[,!(names(combd)=="year")]
  # Merging with cases filters dates to nbefore, nafter and removes polygons with zero cases
  #return(list("a" = combd, "b" = cases))
  combd <- merge(combd, cases, by = c("HEXID", "date"))
  combd$date <- as.numeric(combd$date)
  # Add 000 to fire polygon id to distinguish from control ids
  combd$HEXID[grepl("FR", combd$HEXID)] <- paste0(combd$HEXID[grepl("FR", combd$HEXID)], "000")
  combd$NumId <- as.integer(regmatches(combd$HEXID, regexpr("\\d+", combd$HEXID)))
  fromLast <- is.na(combd[1, "mesma"])
  # some mesma values start with NA so there's no previous value to replace with,
  # in this case, replace with next value
  combd <- apply(combd, 2, zoo::na.locf, fromLast = fromLast)
  combd <- data.frame(combd)
  numerics <- lapply(combd[,-1], as.numeric)
  df <- data.frame(HEXID = combd$HEXID, numerics)
  #df$population <- log10(df$population) 
  return(df)
}

###########################################################################################
############ Driver code #################################################################
###########################################################################################
setwd("~/Desktop/CocciWildfires/ExtractGaps/")
options(stringsAsFactors = F)

# 20km
chexes20 <- readOGR("../hexagons/HexControls20km")
thexes20 <- readOGR("../hexagons/FireHexes20km/")
chexes20 <- spTransform(chexes20, CRS("+init=epsg:3488"))
thexes20 <- spTransform(thexes20, CRS("+init=epsg:3488"))

fix <- function(path){
  df <- read.csv(path)
  df <- df[, 2:ncol(df)]
  df$HEXID <- gsub("TR", "FR", df$HEXID)
  return(df)
}

wind20 <- read.csv("covariates/Wind5_20km.csv")
names(wind20)[4] <- "NDaysAbove5ms"
prism20 <- read.csv("covariates/Prism_20km.csv")
#pop20 <- read.csv("covariates/Population.csv")
pop20 <- read.csv("covariates/GpwPopulation_20km.csv")
mesma20 <- read.csv("covariates/MESMA_20km.csv")
impervious20 <- read.csv("covariates/Impervious_20km.csv")
impervious20 <- impervious20[,2:ncol(impervious20)]
ctcases20 <- read.csv("covariates/ControlCases_20km.csv")
fcases20 <- read.csv("covariates/FireCases_20km.csv")
soil20 <- read.csv("covariates/Soil_20km.csv")
controlFires20 <- read.csv("covariates/ControlHexesWFires_20km.csv")
maxControls20 <- read.csv("covariates/HighlyWeightedControls20km.csv")
fcases20$OnsetMonth <- fcases20$OnsetMonth + 1

# 15km
chexes15 <- readOGR("../hexagons/HexControls15km")
thexes15 <- readOGR("../hexagons/FireHexes15km")
wind15 <- read.csv("covariates/Wind_15km.csv")
prism15 <- read.csv("covariates/Prism_15km.csv")
pop15 <- read.csv("covariates/Population_15km.csv") 
mesma15 <- read.csv("covariates/MESMA_15km.csv")
soil15 <- read.csv("covariates/Soil_15km.csv")
impervious15 <- read.csv("covariates/Impervious_15km.csv")
ctcases15 <- read.csv("covariates/ControlCases_15km.csv")
fcases15 <- read.csv("covariates/FireCases_15km.csv")
fcases15$OnsetMonth <- fcases15$OnsetMonth + 1
controlFires15 <- read.csv("covariates/ControlHexesWFires_15km.csv")
## There are too many control regions when using 15km hexes so I choose the 
## 60 most populous regions
meanp <- aggregate(population~HEXID, pop15[grepl("CT", pop15$HEXID),], base::mean)
hexes60 <- meanp[order(meanp$population),][(nrow(meanp)-59):nrow(meanp), "HEXID" ]
chexes15 <- chexes15[chexes15$HEXID%in%hexes60,]

# 25km
chexes25 <- readOGR("../hexagons/HexControls25km")
thexes25 <- readOGR("../hexagons/FireHexes25km")
wind25 <- read.csv("covariates/Wind_25km.csv")
prism25 <- read.csv("covariates/Prism_25km.csv")
pop25 <- read.csv("covariates/Population_25km.csv")
mesma25 <- read.csv("covariates/MESMA_25km.csv")
impervious25 <- read.csv("covariates/Impervious_25km.csv")
ctcases25 <- read.csv("covariates/ControlCases_25km.csv")
fcases25 <- read.csv("covariates/FireCases_25km.csv")
fcases25$OnsetMonth <- fcases25$OnsetMonth + 1
soil25 <- read.csv("covariates/Soil_25km.csv")
controlFires25 <- read.csv("covariates/ControlHexesWFires_25km.csv")

# To handle different firedates from mergings
getFireDate <- function(fireid) {
  modifiedId <- fireidForDates(fireid)
  if (length(modifiedId) == 2) { 
    fireid <- modifiedId[2] # Define onset at date of second fire (don't start matching until all fires have started)
  }
  onsetmonth <- fireOnsetMonth(fireid)
  date <- convertDate(onsetmonth, onset = T)
  #date <- str_sub(fireid, start = -8, end = -2)
  #date <- as.Date(str_sub(fireid, start = -8), "%Y%m%d")
  return(date)
}


# Used to determine when fires overlapped in the post period
fireOnsetMonth <- function(fireid) {
  onsetmonth <- fcases20[fcases20$fire_id == fireid, "OnsetMonth"]
  return(as.integer(onsetmonth))
}

NameFire <- function(hexids) {
  names <- c()
  for (hexid in hexids) {
    row <- thexes20 %>%as.data.frame() %>%  filter(HEXID == hexid)
    months <-  c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec")
    m <- lubridate::month(row$igntn_d)
    county <- gsub(" County", "", row$county)
    name <- paste0(county, ", ", months[m], ".", row$year)
    if( hexid == "FR41") {
      name <- paste(name, "(2)")
    } else if (hexid == "FR7") {
      name <- 'San Joaquin, Jun.2009'
    }
    names <- c(names, name)
  }
  return(names)
}

hexidToHexid20 <- function(hexids, size) {
  thexes <- thexes20
  hexid20 <- c()
  if (size == "15km") {
    thexes <- thexes15
  } else if (size == "25km") {
    thexes <- thexes25
  }
  for (hexid in hexids) {
    fireid <- hexidToFireid(hexid, thexes)
    hexid20 <- c(hexid20, fireidToHexid(fireid, thexes20))
  }
  return(hexid20)
}

