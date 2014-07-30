
#script to process 2012 ice data and attach spatio-temporal covariates as
#well as several additional structures for space-time modeling; in particular
#the output Data$Which.distances holds the entries for which the kernel
#redistribution matrix is positive, and Data$Dist.entries holds the actual distance
#values for those cells (in meters).  

library(devtools)
install_github('evaluate','hadley')

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(nPacMaps)  #from Josh London
library(maptools)
library(gpclib)
library(automap)
library(ggplot2)
library(Matrix)
#library(bass)

source("c:/users/paul.conn/git/bass/bass/R/bass.R")
MEAN_ADJUST=TRUE  #if TRUE, standardizes "distance to" 

r<-raster("//nmfs/akc-nmml/Polar/Data/Environ/SeaIce/SSMI_SIC/2012/nt_20120401_f17_nrt_n.bin.reproj.tif")

transform_ice <- function(x) {
  ifelse(x < 251, x/2.5, NA)
}
sic_raster <- calc(r, transform_ice)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
sic_raster <- projectRaster(sic_raster, res = c(25067.53, 25067.53),
                            crs = laea_180_proj)

data(alaska_dcw)
data(russia_dcw)
alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))

x_min <- -1600000
x_max <- 1500000
y_min <- -4000000
#y_max <- -1600000
y_max <- -2500000

pep_ext <- extent(x_min, x_max, y_min, y_max)
sic_raster <- crop(sic_raster, pep_ext, snap = "near")
dim_raster=dim(sic_raster)
pep_ext <- extent(sic_raster)

plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
plot(pep_ext, col = "red", add = TRUE)

Grid_poly<-rasterToPolygons(sic_raster,na.rm=FALSE) #convert to SpPolyDF for compatibility with rgeos
Land=list(alaska=alaska_dcw,russia=russia_dcw)

#the following takes awhile.  Instead, consider loading cur.Grid.Rdat
Grid_poly=add.prop.land(Grid=Grid_poly,Land=Land)
#remove initial layer
Grid_poly=Grid_poly[,-1]
save.image("Bering_st_Grid_2012.Rdat")


Area_alaska=gArea(alaska_dcw,byid=TRUE)
Area_russia=gArea(russia_dcw,byid=TRUE)
Alaska_mainland=alaska_dcw[which(Area_alaska==max(Area_alaska)),]
Russia_mainland=russia_dcw[which(Area_russia==max(Area_russia)),]

#Attach distance to mainland for each cell (distance from grid cell centroid to landscape polygon)
Grid_points=gCentroid(Grid_poly,byid=TRUE)
Dist_AK=gDistance(Grid_points,Alaska_mainland,byid=TRUE)
Dist_Rus=gDistance(Grid_points,Russia_mainland,byid=TRUE)
Dist_mainland=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
if(MEAN_ADJUST)Dist_mainland=Dist_mainland/mean(Dist_mainland)
Grid_poly[["dist_mainland"]]=Dist_mainland

#Attach distance to land (including islands) for each cell
Dist_AK=apply(gDistance(Grid_points,alaska_dcw,byid=TRUE),2,'min')
Dist_Rus=apply(gDistance(Grid_points,russia_dcw,byid=TRUE),2,'min')
Dist_land=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
if(MEAN_ADJUST)Dist_land=Dist_land/mean(Dist_land)
Grid_poly[["dist_land"]]=Dist_land

#Attach distance to shelf break (1000m depth contour) for each cell
  Shelf_break=readOGR("c:/users/paul.conn/git/bass/shapefiles/1000m_depth_contour.shp",layer="1000m_depth_contour")
Shelf_break=spTransform(Shelf_break, CRS(laea_180_proj))
Shelf_break1=gBoundary(Shelf_break[1,]) #US side of EEZ
Shelf_break2=gBoundary(Shelf_break[2,]) #russia side
#remove part of Shelf break that intersect with land or the EEZ line
Buffered1=gBuffer(alaska_dcw,width=10000,byid=TRUE)
for(i in 1:length(Buffered1))Shelf_break1=gDifference(Shelf_break1,Buffered1[i,])
Buffered2=gBuffer(russia_dcw,width=10000,byid=TRUE)
for(i in 1:length(Buffered2))Shelf_break1=gDifference(Shelf_break1,Buffered2[i,])
for(i in 1:length(Buffered1))Shelf_break2=gDifference(Shelf_break2,Buffered1[i,])
for(i in 1:length(Buffered2))Shelf_break2=gDifference(Shelf_break2,Buffered2[i,])

EEZs=readOGR("c:/users/paul.conn/git/bass/shapefiles/EEZ_Albers83_BS.shp",layer="EEZ_Albers83_BS")
EEZs=spTransform(EEZs, CRS(laea_180_proj))
EEZ_Alaska=EEZs[1,]  #limit to Alaska EEZ
Shelf_break1=gDifference(Shelf_break1,gBuffer(gBoundary(EEZ_Alaska),width=30000))
Shelf_break2=gDifference(Shelf_break2,gBuffer(gBoundary(EEZ_Alaska),width=30000))

Dist_shelf1=apply(gDistance(Grid_points,Shelf_break1,byid=TRUE),2,'min')
Dist_shelf2=apply(gDistance(Grid_points,Shelf_break2,byid=TRUE),2,'min')
Dist_shelf=apply(cbind(Dist_shelf1,Dist_shelf2),1,'min')

if(MEAN_ADJUST)Dist_shelf=Dist_shelf/mean(Dist_shelf)
Grid_poly[["dist_shelf"]]=Dist_shelf

#not doing depth right now - would need to update from bathymetry (see J. London's Grid_Bathymetry.pdf)
#load('c:\users\paul.conn\git\bass\bass\inst\Depth_Bering_static.Rdat') 
#Grid_poly$depth <-Depth_Bering_poly[["depth"]]
#if(MEAN_ADJUST)Grid_poly$depth=Grid_poly$depth/mean(Grid_poly$depth) #standradize by dividing by its mean

save(Grid_poly,file="BOSSst_Grid_Bering_static.Rdat")
#load("BOSSst_Grid_Bering_static.Rdat")
Adj=rect_adj(x=dim_raster[2],y=dim_raster[1],byrow=TRUE)
Adj2=rect_adj_RW2(x=dim_raster[2],y=dim_raster[1],byrow=TRUE)


#define indicator vector for which cells are to be included in analysis 
#Include = larger pool of cells we might be interested in
#Include2 = those that will actually be modeled
# we'll use the larger pool to get at the 3 cell buffer
Include=(coordinates(Grid_points)[,2]<(-2500000)) 
Include[which(coordinates(Grid_points)[,1]<(-300000))]=0
Include2 = Include
Include2[which(Grid_poly[["land_cover"]]>0.99)]=0  #dont include points with > 99% land
Include2[which(coordinates(Grid_points)[,2]>(-2630000))]=0 
#take out some cells in the bottom right quadrant that are over open water but near land so they have NA ice values
Include2[which(coordinates(Grid_points)[,1]>800000 & coordinates(Grid_points)[,2]<(-3700000))]=0 
Include2[which(coordinates(Grid_points)[,1]>1000000 & coordinates(Grid_points)[,2]<(-3650000))]=0 
Include2[which(coordinates(Grid_points)[,1]>1100000 & coordinates(Grid_points)[,2]<(-3600000))]=0 
Include2[which(coordinates(Grid_points)[,1]>1380000)]=0


#add back in cells over nunavak island
Include2[which(coordinates(Grid_points)[,1]>700000 & coordinates(Grid_points)[,1]<850000 & coordinates(Grid_points)[,2]>(-3300000) & coordinates(Grid_points)[,2]<(-3150000))]=1
#scroll through sea ice data to determine which cells always occur in open water
Date=c(paste("040",c(7:9),sep=''),paste("04",c(10:30),sep=''),paste("050",c(1:9),sep=''),paste("05",c(10:12),sep='')) #limit dates to 4/4-5/22 (first-last BOSS flights in 2012) NOTE: ice data 5/15-5/19 missing
#Date=c(paste("04",c(20:27),sep='')) #limit dates to 4/20-4/27 (first-last BOSS flights used in power analysis)
str1="//nmfs/akc-nmml/Polar/Data/Environ/SeaIce/SSMI_SIC/2012/nt_2012"
str2="_f17_nrt_n.bin.reproj.tif"
Ice1=rep(0,length(Include))
for(idate in 1:length(Date)){
  filename=paste(str1,Date[idate],str2,sep='')
  r<-raster(filename)
  tmp_raster <- calc(r,transform_ice)
  tmp_raster <- projectRaster(tmp_raster, res = c(25067.53, 25067.53),
                              crs = laea_180_proj)
  tmp_raster <- crop(tmp_raster, pep_ext, snap = "near")
  tmp_raster[is.na(tmp_raster)]=1
  Ice1=Ice1+values(tmp_raster>0)
}
Include2[which(Ice1==0)]=0
I.Alaska=as.vector(gIntersects(Grid_poly,EEZ_Alaska,byid=TRUE))
Include2[which(I.Alaska==0)]=0

Grid_reduced=Grid_poly[Include==1,]
Grid_reduced2=Grid_poly[Include2==1,]
plot(Grid_reduced)
plot(Grid_reduced2,col='blue',add=TRUE)

#Grid_reduced2 now has all cells we really want to analyze, but we need to
#select a 3 cell buffer around it (using cells from Grid_reduced) to help
#in computing with a equi-sized redistribution kernel NOT REALLY NEEDED

#Buffered=gBuffer(gCentroid(Grid_reduced2,byid=TRUE),width=75000) 
#I.inter=as.vector(gIntersects(Grid_reduced,Buffered,byid=TRUE))
#Grid_reduced3=Grid_reduced[I.inter==TRUE,]
#plot(Grid_reduced3)
#plot(Grid_reduced2,col="blue",add=TRUE)

#kernel distance - effective distance from one cell to another
Distances=gDistance(gCentroid(Grid_reduced2,byid=TRUE),gCentroid(Grid_reduced2,byid=TRUE),byid=TRUE)
Distances=Matrix(Distances)
#the following neighborhood uses all cells that intersect with a 75km radius circle surrounding a given grid cells centroid
Distances[which(Distances<0.01)]=9
Distances[which(Distances>25067.5 & Distances<25067.6)]=8
Distances[which(Distances>35450 & Distances<35451)]=7
Distances[which(Distances>50135 & Distances<50136)]=6
Distances[which(Distances>56052 & Distances<56053)]=5
Distances[which(Distances>70901 & Distances<70902)]=4
Distances[which(Distances>75202 & Distances<75203)]=3
Distances[which(Distances>79270 & Distances<79271)]=2
Distances[which(Distances>90382 & Distances<90383)]=1
Distances[which(Distances>90383)]=NA  #replace all distances >90383 with NA (not part of my kernel)
Which.distances=which(is.na(Distances)==FALSE)
Dist.entries=Distances[Which.distances]


Adj_reduced=Adj[Include2==1,Include2==1]
Adj2_reduced=Adj2[Include2==1,Include2==1]
plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
plot(pep_ext, col = "red", add = TRUE)
plot(Grid_reduced2,add=TRUE)
plot(EEZ_Alaska,col="yellow",add=TRUE)

Data=list(Adj=Adj_reduced,Adj2=Adj2_reduced,Which.distances=Which.distances,Dist.entries=Dist.entries)
Data$Grid=vector("list",length(Date))  #allocate space for one SpatialPolygonsDataframe for each date in study
rownames(Grid_reduced2@data)=c(1:nrow(Grid_reduced2))
for(irow in 1:nrow(Grid_reduced2))Grid_reduced2@polygons[[irow]]@ID=as.character(irow)
for(idate in 1:length(Date))Data$Grid[[idate]]=Grid_reduced2

filename=paste(str1,Date[1],str2,sep='')
r<-raster(filename)
tmp_raster <- calc(r, transform_ice)
tmp_raster <- projectRaster(tmp_raster, res = c(25067.53, 25067.53),
                            crs = laea_180_proj)
tmp_raster <- crop(tmp_raster, pep_ext, snap = "near")
ice_conc=values(tmp_raster)[Include2==1]
sum(is.na(ice_conc))
sic_poly<-rasterToPolygons(tmp_raster,na.rm=FALSE) #convert to SpPolygonsDF (conversion to Points removes NAs)
ice.poly=sic_poly[Include2==1,]
ice.poly@data$id=rownames(ice.poly@data)
library(ggplot2)
library(plyr)
tmp1<-fortify(ice.poly,region="id")
tmp2<-join(tmp1,ice.poly@data,by="id")
ggplot(tmp2)+aes(long,lat,fill=value)+geom_raster()

#krige missing ice covariates
ice.conc=sic_poly@data[,1]
sic_points=gCentroid(sic_poly,byid=TRUE)
Include=Include2
Include2[is.na(ice.conc)]=0
Include3=Include
Include3[is.na(ice.conc)==FALSE]=0  
n.na=sum(Include3)
na.dm=data.frame(matrix(rep(0,n.na),ncol=1))
colnames(na.dm)="ice.conc"
ice.conc=data.frame(matrix(ice.conc[Include2==1],ncol=1))
colnames(ice.conc)="ice.conc"
coords=coordinates(sic_points)[Include2==1,]
coords.pred=coordinates(sic_points)[Include3==1,]
rownames(coords)=c(1:nrow(coords))
rownames(coords.pred)=c(1:nrow(coords.pred))
sic_points=SpatialPointsDataFrame(coords,ice.conc,proj4string=CRS(laea_180_proj))
pred_loc=SpatialPointsDataFrame(coords.pred,na.dm,proj4string=CRS(laea_180_proj))
krige_out=autoKrige(ice.conc~1,input_data=sic_points,new_data=pred_loc)$krige_output[["var1.pred"]] 
ice.poly=sic_poly[Include==1,]
ice.poly[[1]][which(is.na(ice.poly[[1]]==TRUE))]=krige_out
ice.poly@data$id=rownames(ice.poly@data)
tmp.df<-fortify(ice.poly,region="id")
tmp.df<-join(tmp.df,ice.poly@data,by="id")
ggplot(tmp.df)+aes(long,lat,fill=value)+geom_raster()    

for(idate in 1:length(Date)){
  filename=paste(str1,Date[idate],str2,sep='')
  r<-raster(filename)
  tmp_raster <- calc(r, transform_ice)
  tmp_raster <- projectRaster(tmp_raster, res = c(25067.53, 25067.53),
                              crs = laea_180_proj)
  tmp_raster <- crop(tmp_raster, pep_ext, snap = "near")
  sic_poly<-rasterToPolygons(tmp_raster,na.rm=FALSE) 
  sic_points=gCentroid(sic_poly,byid=TRUE)
  ice.conc=sic_poly@data[,1]
  Include2=Include
  Include2[is.na(ice.conc)]=0
  Include3=Include
  Include3[is.na(ice.conc)==FALSE]=0  
  n.na=sum(Include3)
  na.dm=data.frame(matrix(rep(0,n.na),ncol=1))
  colnames(na.dm)="ice.conc"
  ice.conc=data.frame(matrix(ice.conc[Include2==1],ncol=1))
  colnames(ice.conc)="ice.conc"
  coords=coordinates(sic_points)[Include2==1,]
  coords.pred=coordinates(sic_points)[Include3==1,]
  rownames(coords)=c(1:nrow(coords))
  rownames(coords.pred)=c(1:nrow(coords.pred))
  sic_points=SpatialPointsDataFrame(coords,ice.conc,proj4string=CRS(laea_180_proj))
  pred_loc=SpatialPointsDataFrame(coords.pred,na.dm,proj4string=CRS(laea_180_proj))
  krige_out=autoKrige(ice.conc~1,input_data=sic_points,new_data=pred_loc)$krige_output[["var1.pred"]] 
  ice.poly=sic_poly[Include==1,]
  ice.poly[[1]][which(is.na(ice.poly[[1]]==TRUE))]=krige_out
  ice_conc=ice.poly[[1]]/100
  if(length(which(ice_conc<0))>0)ice_conc[which(ice_conc<0)]=0
  if(length(which(ice_conc>1))>0)ice_conc[which(ice_conc>1)]=1
  Data$Grid[[idate]]=spCbind(Data$Grid[[idate]],ice_conc)
}

#distance to countour
tmp.raster<-raster(nrow=dim_raster[1],ncol=dim_raster[2])
extent(tmp.raster)=extent(Data$Grid[[1]])
for(idate in 1:length(Date)){
  tmp.raster=rasterize(Data$Grid[[idate]],tmp.raster,'ice_conc')
  sic_contour=rasterToContour(tmp.raster,levels=c(0.1))
  dist_contour=as.vector(gDistance(Grid_points[Include==1,],sic_contour,byid=TRUE))
  if(MEAN_ADJUST)dist_contour=dist_contour/mean(dist_contour)
  Data$Grid[[idate]]=spCbind(Data$Grid[[idate]],dist_contour)
}  

#distance to southern ice edge
for(idate in 1:length(Date)){
  edge.line=get.edge.Bering(SpDF=Data$Grid[[idate]],proj=laea_180_proj)        
  dist_edge=as.vector(gDistance(Grid_points[Include==1,],edge.line,byid=TRUE))
  if(MEAN_ADJUST)dist_edge=dist_edge/mean(dist_edge)
  Data$Grid[[idate]]=spCbind(Data$Grid[[idate]],dist_edge)
}  

#ecoregion
Ecoreg=readOGR("c:/users/paul.conn/git/bass/shapefiles/Marine_Ecoregions_AK.shp",layer="Marine_Ecoregions_AK")
Ecoreg=spTransform(Ecoreg, CRS(laea_180_proj))
Ecoregion=over(Data$Grid[[1]],Ecoreg)
Ecoregion=Ecoreg[["ECOREGION"]][Ecoregion]
Ecoregion[which(Ecoregion==12)]=24  #move one cell that was given 
for(idate in 1:length(Date))Data$Grid[[idate]]=spCbind(Data$Grid[[idate]],Ecoregion)

#metadata
Data$Meta=list(date.start="7Apr2012",date.end="12May2012")
Data$Meta$info="2012 PEP BOSS effective area survey grid data for spatio-temporal modeling"

#save
save(Data,file="BOSSst_2012data.Rdata")