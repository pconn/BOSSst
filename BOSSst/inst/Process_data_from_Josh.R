#Calculate effective area surveyed based on field of view calculations;
# the approach works by joining the footprints associated with each photograph
library(sp)
library(rgeos)


load('./data_from_JML/boss_hotspots_sp.Rda') #hotspots
load('./data_from_JML/boss_geo_sp.Rda') #on effort points
load('./data_from_JML/boss_grid_env.Rda') # effort summary (without area surveyed)
load('./data_from_JML/grid_spdf.rda') # BOSS Grid

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
#reproject 
grid_spdf=spTransform(grid_spdf, CRS(laea_180_proj))
boss_geo=spTransform(boss_geo, CRS(laea_180_proj))
boss_hotspots_sp=spTransform(boss_hotspots_sp, CRS(laea_180_proj))


year_func<-function(x)strsplit(x,'_')[[1]][1]

aircraft_func<-function(x){
  str1 = strsplit(x,'_')[[1]][2]
  strsplit(str1,"F")[[1]][1]
}

#get a BOSS grid ID for each on effort point by overlaying points on grid
Study.area=gUnionCascaded(grid_spdf)
int=gIntersects(boss_geo,Study.area,byid=TRUE)
boss_geo=boss_geo[apply(int,2,'sum')==1,]  #only apply to points that intersect the grid
save(boss_geo,file="./boss_geo_sp_tmp1.Rda")

load("./boss_geo_sp_tmp1.Rda")
my.which<-function(x)which(x==1)  #run gIntersects in blocks of 1000 (slightly faster than 10000, MUCH faster than 1)
Which_cell = rep(0,length(boss_geo))
for(i in 1:1754){
  Tmp_range = c(((i-1)*1000+1):(i*1000))
  Tmp_mat = gIntersects(boss_geo[Tmp_range,],grid_spdf,byid=TRUE)
  Which_cell[Tmp_range]=apply(Tmp_mat,2,'my.which')
}
Tmp_range=c(1754001:1754038)
Tmp_mat = gIntersects(boss_geo[Tmp_range,],grid_spdf,byid=TRUE)
Which_cell[Tmp_range] = apply(Tmp_mat,2,'my.which')
#translate from which() to cell ID
ID_grid = getSpPPolygonsIDSlots(grid_spdf) #ID name for each grid cell
ID_photo_grid = ID_grid[Which_cell]
boss_geo[["GridID"]]=ID_photo_grid
save(boss_geo,file="./boss_geo_sp_tmp2.Rda")

load(boss_geo,file="./boss_geo_sp_tmp2.Rda")
for(ifl in 1:length(boss_grid_env)){
  Cur_pts <- boss_geo[which(boss_geo$flightid==boss_geo$flightid[1]),]
  Aircraft=sapply(as.character(points_fl1),aircraft_func)


# x - longitude
# y - latitude
# aircraft - 
get_vertices <- function(Pts){
  Aircraft <- substr(Pts,)
  
  
}

for(i in 1:length(points_fl1)){
  #get vertices
}