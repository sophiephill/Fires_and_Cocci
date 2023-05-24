
library(sf)
library(rgdal)
library(sp)
library(rgeos)
library(hexbin)

setwd("~/Desktop/CocciWildfires/")
thexes <- read.csv("thexes_data.csv",stringsAsFactors = F)
fires <- readOGR("shapes/92FiresBoundaries", "fires92",stringsAsFactors = F)
fires <- spTransform(fires, CRS("+init=epsg:3488 +units=m"))

#### Get all fires 2003-2018 ###################################################################
occpoints <- readOGR("data/wildfire_data/S_USA.MTBS_BURN_AREA_BOUNDARY.shp")
occpoints <- occpoints[occpoints$FIRE_TYPE=="Wildfire",]
occpoints0315 <- occpoints[occpoints$YEAR %in% as.character(2003:2015),]
counties <- readOGR("~/Desktop/CocciWildfires/data/CA_counties","CA_Counties_TIGER2016")
counties <- spTransform(counties, proj4string(occpoints0315))
subset <-c("Kern", "Kings", "Fresno", "San Joaquin","San Luis Obispo", "Tulare")
counties6 <- counties[counties$NAME %in% subset,]
insex <- gIntersects(counties6, occpoints0315, byid = T)
rownames(insex) <- occpoints0315$FIRE_ID
colnames(insex) <- counties6$NAME
insex <- insex[rowSums(insex)>0, ]
ids <- rownames(insex)
tids <- thexes$fire_id
omits <- ids[!(ids %in% tids)]


#### Grouping Sizes ###################################################################
setwd("~/Desktop/CocciWildfires")
thexes <- thexes[order(thexes$acres, decreasing = T),]
quantiles <- quantile(thexes$acres)
jenks <- classInt::classIntervals(thexes$acres, 3, style = "jenks")$brks
thexes$JenksSizeCat <- cut(thexes$acres, breaks = jenks, 
                       labels = c("S","M","L"),
                       include.lowest = T)
(ggplot(thexes,aes(x=reorder(FID, -acres), y=acres)) 
  +geom_bar(stat = "identity", colour='black', fill='gray') 
  + theme_bw() + xlab("fire id") 
  + geom_hline(yintercept = jenks[2:4], color = "red",size=.8)
  #+ geom_hline(yintercept = 2000, color = "red",size=.8)
  + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  + ggtitle("Jenks Breaks"))
# I think jenks looks good

### Choosing sizes #####################################################################
getHexSideLengthfromArea <- function(area) {
  s2 <- (2/3)*area / sqrt(3)
  return(sqrt(s2))
}

getInscribedCircleRadiusfromArea <- function(area) {
  side <- getHexSideLengthfromArea(area)
  return(side*sqrt(3) / 2)
}
getSideLengths <- function(poly) {
  vertices <- poly@polygons[[1]]@Polygons[[1]]@coords
  for(i in 1:(nrow(vertices)-1)) {
    diff <- vertices[i,]-vertices[i+1,]
    diff <- sqrt(sum(diff**2))
    print(diff)
  }
}
acres2m2 <- function(a) 4047*a
area <- gArea(h22[1,]) / 1000000
getHexSideLengthfromArea(area)

g1 <- thexes[thexes$JenksSizeCat=="S", "acres"]
g2 <- thexes[thexes$JenksSizeCat=="M", "acres"]
g3 <- thexes[thexes$JenksSizeCat=="L", "acres"]
g4 <- thexes[thexes$JenksSizeCat=="XL", "acres"]

######################################################
#### Trying different radii around the fires 
######################################################
biggestfire <- fires[fires$FIRE_NAME=="PIUTE",]
# measures IF fire boundary were hexagonal, don't take this too seriously
side <- getHexSideLengthfromArea(acres2m2(biggestfire$ACRES))
r <- getInscribedCircleRadiusfromArea(acres2m2(biggestfire$ACRES)) 

# fire buffer
centroid <-  gCentroid(biggestfire, byid=TRUE)
coords <- centroid@coords
#with this ratio of dx,dy, you get a regular hexagon = all sides same length
#this is preferable because then you have a formula for the inscribed radius
dytodx <- 3.4641/2
dy <- 15000/dytodx
hexagons <- hexpolygon(coords[,"x"], coords[,"y"], dx=dytodx*dy, dy=dy)
hexmat <- matrix(c(hexagons$x, hexagons$y), ncol=2)
hexmat <- matrix(c(hexmat[,1], hexmat[1, 1], hexmat[, 2], hexmat[1, 2]), ncol=2)
hexpoly <- SpatialPolygons(list(Polygons(list(Polygon(hexmat)), 1)),
                           proj4string= CRS("+init=epsg:3488"))
hexradius <- getInscribedCircleRadiusfromArea(gArea(hexpoly))
plot(hexpoly)
plot(biggestfire, add=T)

fires <- fires[order(fires$ACRES),]
percents <- data.frame(acres =fires$ACRES, 
                       perc15 = 100*sapply(slot(fires, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area")[1])/gArea(hexpoly))
percents
#verify equal side lengths
getSideLengths(hexpoly)

#### What if I used circle buffers
sfcentroid <- st_sfc(st_point(centroid@coords[1,]))
circle <- st_buffer(sfcentroid, hexradius*1.05)
hexarea <- gArea(hexpoly)
circarea <- st_area(circle)
plot(circle)
plot(biggestfire, add = T)
plot(hexpoly, add = T)

#grid
hexside <- getHexSideLengthfromArea(gArea(hexpoly))
r <- getInscribedCircleRadiusfromArea(gArea(hexpoly))
HexPts <-spsample(border2, type="hexagonal",cellsize=2*r)
hexes <- HexPoints2SpatialPolygons(HexPts)
hexes <- spTransform(hexes, proj4string(biggestfire))
plot(hexes)
plot(hexpoly,add=T, col = "red")

###########################################################################################
##### Making Fire hexagons #################################################################
###########################################################################################

df <- fires@data[,1:7]
pre0 <- ifelse(nchar(df$STARTMONTH)==1, "0", "")
df$STARTMONTH <- paste0(pre0, df$STARTMONTH)
df$igntn_d <- paste(df$YEAR, df$STARTMONTH, df$STARTDAY,sep="/")
names(df) <- tolower(names(df))
df <- plyr::join(df,thexes[,c("fire_id", "mn_svrt","md_svrt","JenksSizeCat")], by = "fire_id")
df <- df[,c("fire_id", "year", "igntn_d","acres","fire_type","JenksSizeCat")]

polylist <- list()
hexsize <- 25000
for(i in 1:nrow(fires)) {
  fire <- fires[i,]
  side <- getHexSideLengthfromArea(acres2m2(fire$ACRES))
  r <- getInscribedCircleRadiusfromArea(acres2m2(fire$ACRES)) 
  coords <- gCentroid(fire, byid=TRUE)@coords
  dytodx <- 3.4641/2
  #dy <- side*1.5 ## makes hexagons vary with fire size
  dy <- hexsize/dytodx 
  hexagons <- hexpolygon(coords[,"x"], coords[,"y"], dx=dytodx*dy, dy=dy)
  hexmat <- matrix(c(hexagons$x, hexagons$y), ncol=2)
  hexmat <- matrix(c(hexmat[,1], hexmat[1, 1], hexmat[, 2], hexmat[1, 2]), ncol=2)
  hexpoly <- Polygons(list(Polygon(hexmat)), i)
  polylist <- c(polylist, hexpoly)
}
hexdf <- SpatialPolygonsDataFrame(SpatialPolygons(polylist, proj4string= CRS("+init=epsg:3488")), df)
hexdf$HexArea<-sapply(slot(hexdf, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))
writeOGR(hexdf, "hexagons/FireHexes25km", layer = "FireHexes25km", driver = "ESRI Shapefile",
         overwrite_layer = T)

###########################################################################################
### Forming Control grid ############################################################################
###########################################################################################
counties <- readOGR("~/Desktop/research/wildfire\ map/CA_counties",
                    "CA_Counties_TIGER2016")
subset <-c("Kern", "Kings", "Fresno", "San Joaquin","San Luis Obispo", "Tulare")
counties6 <- counties[counties$NAME %in% subset,]
border <- rgeos::gUnaryUnion(counties6)
border2 <- spTransform(border, CRS("+init=EPSG:3488"))

#### Identically sized controls
area <- hexdf@data[1, "HexArea"]
r <- getInscribedCircleRadiusfromArea(area)
HexPts <-spsample(border2, type="hexagonal",cellsize=2*r)
hexes <- HexPoints2SpatialPolygons(HexPts)
ids <- sapply(slot(hexes, "polygons"), function(x) slot(x, "ID"))
hexes <- SpatialPolygonsDataFrame(hexes, data.frame(FID=ids, row.names=ids))
writeOGR(hexes, "hexagons/HexControls15km", layer = "HexControls15km", 
         driver = "ESRI Shapefile",overwrite_layer = T)
