#plot_survey_area.R


library(sp)
library(rgdal)
library(nPacMaps)  #from Josh London
library(rgeos)
library(doBy)
library(ggplot2)
library(plyr)
library(RPostgreSQL)
library(sf)


#load('./data_from_JML/boss_grid_env.Rda') # effort summary (without area surveyed)
load('c:/users/paul.conn/git/STabundance/Data_for_ST_plot.Rdata')
load('boss_geo_sp_tmp5.Rda')  #produced in BOSShotspots package where turn filter etc. is used
load('boss_flt_dates.Rda')

Cur.grid=Data$Grid[[1]]

#sample 1 out every 10 photographs to reduce plotting burden
Photo_pts = boss_geo[c(1:171352)*10,]

#remove flts before or after study window, fast ice flights
Photo_pts = Photo_pts[-which(Photo_pts$flightid %in% as.character(Flt_table$Flt[c(19,34:39,57,70,71,73)])),]

#remove photos that don't intersect with grid
I.inter=gIntersects(Photo_pts,Cur.grid,byid=TRUE)
Points=Photo_pts[which(apply(I.inter,2,'sum')>0),]  #only include effort points that are 'on grid'

#data(alaska_dcw)
#data(russia_dcw)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))


ak<-nPacMaps::fortifi(alaska_dcw,tol=2,minarea=1000)
rus<-nPacMaps::fortifi(russia_dcw,tol=2,minarea=1000)

boss_grid <- fortify(Cur.grid)
#boss_grid <- join(boss_grid,Cur.grid@data,by="id")



# use the nPackMaps ggExpansion function to 
# expand the plot extent by 0.1 in x and y
#lims <- nPacMaps::ggExpansion(Cur.grid,x="x","y",x_fac=0.1,y_fac=0.1)
#xlim <- lims$xlim
#ylim <- lims$ylim


Points.df=data.frame(Points)
Yr.tmp = matrix(unlist(strsplit(Points.df$flightid,'_')),2,nrow(Points.df))
Yr.tmp = Yr.tmp[1,]
Yr.tmp[Yr.tmp=="13"]="2013"
Yr.tmp[Yr.tmp=="12"]="2012"
Points.df$Year = as.factor(Yr.tmp)



#ggplot()+geom_point(data=Points.df,aes(x=GPSLONG_INT,y=GPSLAT_INT))

shelf <- SpatialLinesDataFrame(gBoundary(Shelf_break[1,]),data.frame(id=1))
eez <- SpatialLinesDataFrame(gBoundary(EEZ_Alaska),data.frame(id=1))

# fortify shelf and eez
shelf <- fortify(shelf)
eez <- fortify(eez)

# use some colorbrewer colors
red <- "#e41a1c"
blue <- "#377eb8"
orange <- "#ff7f00"
brown <- "#a65628"

# make our plot
p <- ggplot() + geom_polygon(data=boss_grid,
                             aes(x=long,y=lat,group=group),
                             fill="white",color="gray")  
p <- p + geom_path(data=shelf,
                   aes(x=long,y=lat,group=group),
                   color="orange",size=0.75)
p <- p + geom_path(data=eez,
                   aes(x=long,y=lat,group=group),
                   color="brown",size=1)
p <- p + geom_polygon(data=ak,
                      aes(x=long,y=lat,group=group),
                      fill="grey60",color="grey60",size=1.2)
p <- p + geom_polygon(data=rus,
                      aes(x=long,y=lat,group=group),
                      fill="grey60",color="grey60",size=1.2)
p <- p + geom_point(data=Points.df,
                    aes(x=coords.x1,y=coords.x2),color="darkblue",
                    size=0.005)

xlim=c(-230000,1400000)
ylim=c(-3750000,-2500000)
p <- p + coord_equal(xlim=xlim,ylim=ylim)
#p <- p + scale_x_continuous(expand=c(0,0),labels=nPacMaps::to_km())
#p <- p + scale_y_continuous(expand=c(0,0),labels=nPacMaps::to_km())
p <- p + labs(x="Easting",y="Northing")
p <- p + theme(text=element_text(size=16),title=element_text(size=16),axis.ticks = element_blank(),axis.text=element_blank())
p <- p + facet_grid(Year ~ .)

#p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) #remove slash in legend
p

# save a pdf
#ggsave(file="BOSS_survey_ribbon.tiff",dpi=300)

pdf(file="BOSS_effort_map.pdf")
p
dev.off()



