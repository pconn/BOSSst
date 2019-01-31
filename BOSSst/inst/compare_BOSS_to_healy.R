# compare BOSS estimates to Healy estimates

library(sp)
library(rgdal)
library(rgeos)

Healy_area = readOGR(dsn="./Healy",layer="AMSR_4x4_StudyArea")
#reproject
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0","+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Healy_area <- spTransform(Healy_area, CRS(laea_180_proj))

load('./AlaskaBeringData2012_2013_14Dec2015.Rdat')  #boss grid w covariates by day

Grid = Data$Grid$y2012[[1]]

plot(Grid)
plot(Healy_area,col='blue',add=TRUE)
Healy_comb = gUnionCascaded(Healy_area)

Intersects = gIntersects(Grid,Healy_comb,byid=TRUE)
IDs = which(Intersects)

load("tweedie_ObsVar_2012.RData")  #


t_steps = 29
D=rep(0,29)
for(i in 1:29){
  n_cells = ncol(Out$Report$Z_s)/t_steps
  pred_day = i
  Which_cells =  (n_cells)*(pred_day-1)+IDs
  N = rowSums(Out$Report$Z_s[,Which_cells])
  D[i] = (N/625/length(Which_cells))[3]
}

load("tweedie_ObsVar_Day_2013.RData")  #
t_steps = 28
D=rep(0,t_steps)
for(i in 1:t_steps){
  n_cells = ncol(Out$Report$Z_s)/t_steps
  pred_day = i
  Which_cells =  (n_cells)*(pred_day-1)+IDs
  N = rowSums(Out$Report$Z_s[,Which_cells])
  D[i] = (N/625/length(Which_cells))[3]
}
