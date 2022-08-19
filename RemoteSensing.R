pck <- (c("tidyr","rgdal","ggplot2","raster","leaflet","rasterVis","gridExtra", "maptools","RColorBrewer","plotly", "rgeos"))
new_pck <- pck[!pck %in% installed.packages()[,"Package"]]
if(length(new_pck)){install.packages(new_pck)}

sapply(pck, require, character.only=TRUE)

setwd("/Users/mostafa/Desktop/Github/RS/ProS2L2A")

Study_Area<-readOGR("/Users/mostafa/Desktop/Github/New Folder/saa.shp")
point_SA<-readOGR("/Users/mostafa/Desktop/pointtt/points.shp")


Study_Area2 <- spTransform(Study_Area, CRS("+proj=longlat +datum=WGS84"))

m <- leaflet(sizingPolicy = leafletSizingPolicy(defaultHeight = 200, viewer.suppress = TRUE, knitr.figure=FALSE)) %>%
  addProviderTiles(providers$Esri.WorldStreetMap) %>%
  addPolygons(data = Study_Area2,
              stroke = FALSE,
              smoothFactor = 0.5
  )
m

#Load data
S2 <- "/Users/mostafa/Desktop/Github/RS/S2L2A"
S2 <- list.files(S2, recursive = TRUE, full.names = TRUE, pattern = "B0[2348]_10m.jp2$")
S2 <- lapply(1:length(S2), function (x) {raster(S2[x])})
S2[1]

#set layout
options(repr.plot.width=41, repr.plot.height=20)
m <- rbind(c(1,2))
layout(m)

#Stack
S2_stack <- stack(S2)

#plot True/False Color
plotRGB(S2_stack, r=3, g=2, b=1, scale=maxValue(S2[[2]]), stretch='hist') #stretch='lin'
plot(Study_Area, add=TRUE, border='yellow', lwd=5)
plotRGB(S2_stack, r=4, g=3, b=2, scale=maxValue(S2[[2]]), stretch='hist') #stretch='lin'
plot(Study_Area, add=TRUE, border='yellow', lwd=5)

#Set layout
options(repr.plot.width = 35, repr.plot.height = 10)
m <- rbind(c(1, 2))
layout(m)
#crop and plot
S2_stack_crop <- crop(S2_stack, Study_Area)
plotRGB(S2_stack_crop, r=3, g=2, b=1, scale=maxValue(S2[[2]]), stretch='hist')
plotRGB(S2_stack_crop, r=4, g=3, b=2, scale=maxValue(S2[[2]]), stretch='hist')

#Derive NDVI
NDVI<-list()

for (i in 1:(length(S2)/4)) {
  NDVI[[i]] <- overlay(S2_stack_crop[[((i-1)*4+3)]], S2_stack_crop[[((i-1)*4+4)]], fun=function(x,y) (y-x) / (y+x))
  names(NDVI[[i]]) <- paste0("NDVI_", strsplit(strsplit(names(S2_stack_crop[[(i-1)*4+4]]), "_")[[1]][2], "T")[[1]][1])
}
NDVI

#Set Layout
options(repr.plot.width = 100, repr.plot.height = 90)

#plot
breaks <- c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)
pal <- brewer.pal(11,"RdYlGn")
mapTheme <- rasterTheme(region = pal)
mapTheme$fontsize$text = 9

#Plot NDVI and Points of reference
levelplot(stack(NDVI), scales=list(draw=FALSE), colorkey=FALSE, par.settings=mapTheme)

NDVI_points <- lapply(NDVI, FUN=function (NDVI) {raster::extract(NDVI, point_SA, method = 'bilinear', df=TRUE)})
NDVI_points

#Combine df
NDVI_points_df <- do.call("cbind", NDVI_points)

#Clean df - remove duplicated columns
NDVI_points_df <- NDVI_points_df[, !duplicated(colnames(NDVI_points_df))]
NDVI_points_df

#Arrange df
NDVI_points_df <-gather(NDVI_points_df, key = Date, value = value, -ID)
NDVI_points_df

#plot NDVI Temporal series
ndvi_plot <- ggplot(data=NDVI_points_df, aes(x=Date, y=value, group=ID, color=ID)) + geom_line() + geom_point()
ndvi_plot

#(#f <- S2_stack_crop
#bandass <- f
#class(bandass)
#slotNames(bandass)
#res(bandass)

#plot(bandass)
#projection(bandass)
#S2_stack_crop <- projectRaster(bandass, crs=CRS("+init=epsg:4326"), method = "bilinear")
#plot(S2_stack_crop))
