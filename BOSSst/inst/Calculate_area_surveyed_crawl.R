#Calculate effective area surveyed based on field of view calculations;
# the approach works by joining the footprints associated with each photograph
library(sp)
library(rgeos)
library(nPacMaps)
library(crawl)

#load('./data_from_JML/boss_hotspots_sp.Rda') #hotspots
load('./data_from_JML/boss_geo_sp.Rda') #on effort points
#load('./data_from_JML/boss_grid_env.Rda') # effort summary (without area surveyed)
load('./data_from_JML/grid_spdf.rda') # BOSS Grid

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

data(alaska_dcw)
#data(russia_dcw)
alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
#russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))

#reproject 
grid_spdf=spTransform(grid_spdf, CRS(laea_180_proj))
boss_geo=spTransform(boss_geo, CRS(laea_180_proj))
#boss_hotspots_sp=spTransform(boss_hotspots_sp, CRS(laea_180_proj))

year_func<-function(x)strsplit(x,'_')[[1]][1]

aircraft_func<-function(x){
  str1 = strsplit(x,'_')[[1]][2]
  strsplit(str1,"F")[[1]][1]
}

#remove locations that don't overlap with grid or are over land
Study.area=gUnionCascaded(grid_spdf)
int=gIntersects(boss_geo,Study.area,byid=TRUE)
boss_geo=boss_geo[apply(int,2,'sum')==1,]  #only apply to images that intersect the grid
cat(paste("eliminated",sum(apply(int,2,'sum')==0),"images that don't intersect grid"))
int1=gIntersects(boss_geo,alaska_dcw[1,],byid=TRUE) #images that intersect the AK mainland
boss_geo=boss_geo[apply(int1,2,'sum')==0,]  
cat(paste("eliminated",sum(apply(int1,2,'sum')==1),"images that intersect mainland"))
int2=gIntersects(boss_geo,alaska_dcw[108,],byid=TRUE) #images that intersect St. Lawrence Island
boss_geo=boss_geo[apply(int2,2,'sum')==0,]  
cat(paste("eliminated",sum(apply(int2,2,'sum')==1),"images that intersect SL"))
int3=gIntersects(boss_geo,alaska_dcw[151,],byid=TRUE) #images that intersect Nunavak
boss_geo=boss_geo[apply(int3,2,'sum')==0,]  
cat(paste("eliminated",sum(apply(int3,2,'sum')==1),"images that intersect Nunavak"))
#remove AeroFl14(2013) starboard due to lens messup (24mm)
boss_geo=boss_geo[-which(boss_geo$flightid=="13_AeroFl14" & boss_geo$camera=="Starboard"),]
#remove AeroFl05 (port camera) - "unrecoverable exposure issues"
boss_geo=boss_geo[-which(boss_geo$flightid=="13_AeroFl05" & boss_geo$side=="Port"),]
save(boss_geo,file="./boss_geo_sp_tmp1.Rda")




#get a BOSS grid ID ("GridID") for each on effort photograph by overlaying points on grid
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

load("./boss_geo_sp_tmp2.Rda")

#reorder points within flightID by date / time
Order_all = rep(0,length(boss_geo))
FlightSide_ID = paste(boss_geo$flightid,boss_geo$side)
FlightSide_unique_IDs = unique(FlightSide_ID)
#boss_geo2 = boss_geo
Coords=matrix(0,length(boss_geo),2)
Data_ordered = boss_geo@data
cur_pl=1
for(ifl in 1:length(FlightSide_unique_IDs)){
  Which = which(FlightSide_ID == FlightSide_unique_IDs[ifl])
  Tmp=boss_geo[Which,]
  Order = order(boss_geo[Which,]$img_dt)
  Coords[cur_pl:(cur_pl+length(Which)-1),]=coordinates(Tmp[Order,])
  Data_ordered[cur_pl:(cur_pl+length(Which)-1),]=Tmp@data[Order,]
  #Order_all[Which]=cur_pl+Order-1
  cur_pl=cur_pl+length(Which)
}

#reform boss_geo SpPointsObject with ordered entries
#New_x = matrix(as.numeric(coordinates(boss_geo)[,1])[Order_all],ncol=1)
#New_y = matrix(as.numeric(coordinates(boss_geo)[,2])[Order_all],ncol=1)
#New_coords = cbind(New_x,New_y)
#Data=boss_geo@data[Order_all,]
#rownames(Data)=c(1:length(boss_geo))
boss_geo=SpatialPointsDataFrame(coords=Coords,data=Data_ordered,proj4str=CRS(laea_180_proj))
save(boss_geo,file='boss_geo_sp_tmp3.Rda')

load('boss_geo_sp_tmp3.Rda')



#produce some unique identifiers for points obtained for each flightID, camera, and pass through a grid cell
FlightCellSide_ID = paste(boss_geo$flightid,boss_geo$GridID,boss_geo$side)
FlightCell_ID = paste(boss_geo$flightid,boss_geo$GridID)
FlightCellSide_unique_IDs = unique(FlightCellSide_ID)

#for each identifier, break flights into "legs" based on number of segments >5 seconds apart [so that different aircraft interpolations can be done for each]
for(iid in 1:length(FlightCellSide_unique_IDs)){
  Which_photos = which(FlightCellSide_ID==FlightCellSide_unique_IDs[iid])
  n_photos=length(Which_photos)
  DT<-boss_geo@data[Which_photos,"img_dt"]
  DT_diff=difftime(DT[2:length(DT)],DT[1:(length(DT)-1)],units="secs")
  Breaks = which(abs(DT_diff)>5)
  if(length(Breaks)==1){
    FlightCellSide_ID[Which_photos[1:Breaks]]=paste0(FlightCellSide_ID[Which_photos[1:Breaks]],'leg1')
    FlightCellSide_ID[Which_photos[(Breaks+1):length(Which_photos)]]=paste0(FlightCellSide_ID[Which_photos[(Breaks+1):length(Which_photos)]],'leg2')
  }
  if(length(Breaks)>1){
    FlightCellSide_ID[Which_photos[1:Breaks[1]]]=paste0(FlightCellSide_ID[Which_photos[1:Breaks[1]]],'leg1')
    for(ibreak in 2:(length(Breaks))){
      FlightCellSide_ID[Which_photos[(Breaks[ibreak-1]+1):Breaks[ibreak]]]=paste0(FlightCellSide_ID[Which_photos[(Breaks[ibreak-1]+1):Breaks[ibreak]]],paste0('leg',ibreak))
    } 
    FlightCellSide_ID[Which_photos[(Breaks[length(Breaks)]+1):length(Which_photos)]]=paste0(FlightCellSide_ID[Which_photos[(Breaks[length(Breaks)]+1):length(Which_photos)]],paste0('leg',length(Breaks)+1))
  }
}

save(boss_geo,FlightCellSide_ID,file='boss_geo_sp_tmp3b.Rda')

#for each "leg", base plane x, y position for duplicated x,y coords based on flight track and frequency of camera firings
FlightCellSide_unique_IDs = unique(FlightCellSide_ID) 
Interp_loc = coordinates(boss_geo)
for(iid in 1:length(FlightCellSide_unique_IDs)){
  Which_photos = which(FlightCellSide_ID==FlightCellSide_unique_IDs[iid])
  n_photos=length(Which_photos)
  if((n_photos>1 & n_photos<10)|iid==8237){
    Coords = coordinates(boss_geo[Which_photos,])
    I_dup = (Coords[2:n_photos,1]==Coords[1:(n_photos-1),1] & Coords[2:n_photos,2]==Coords[1:(n_photos-1),2])
    Duplicated = as.numeric(which(I_dup==1))+1
    if(length(Duplicated)>0){
      #first.pl = 1
      #last.pl = max(which(I_dup==0))  #first and last non-duplicated
      #number of duplicated sequences
      I.str = rep(0,length(Duplicated))
      if(length(Duplicated)>1){
        I.str[2:length(Duplicated)]=(Duplicated[2:length(Duplicated)]==(Duplicated[1:(length(Duplicated)-1)]+1))
      }
      n.strings = sum(I.str==0)
      Begin=End = which(I.str==0)  #beginning and end of each missing (well, erroneous) gps stream
      for(istr in 1:(n.strings-1))End[istr]=Begin[istr+1]-1
      End[n.strings]=length(Duplicated)
      
      #first, treat case where there is only one record with usable GPS info.
      #in this case, can't interpolate or extrapolate; simply add 150m to each lat and long so that will be adding up the total photo area
      if((n_photos-length(Duplicated))==1){
        Coords[2:n_photos,1]=Coords[1,1]+150*(1:length(Duplicated))
        Coords[2:n_photos,2]=Coords[1,2]+150*(1:length(Duplicated))
      }
      else{  #>1 record with usable GPS info: interpolation and/or extrapolation possible
        #fill in last observation if it is a duplicate
        if(Duplicated[End[n.strings]]==n_photos){
          Which.not = c(1,which(I_dup==0)+1)
          Which.not = Which.not[c(length(Which.not)-1,length(Which.not))]
          Coords[n_photos,1] = Coords[Which.not[2],1]+(n_photos-Which.not[2])*(Coords[Which.not[2],1]-Coords[Which.not[1],1])/(Which.not[2]-Which.not[1])
          Coords[n_photos,2] = Coords[Which.not[2],2]+(n_photos-Which.not[2])*(Coords[Which.not[2],2]-Coords[Which.not[1],2])/(Which.not[2]-Which.not[1])
        
          if(End[n.strings]==Begin[n.strings]){
            End=End[1:(n.strings-1)]
            Begin=Begin[1:(n.strings-1)]
            n.strings=n.strings-1
          }
          else End[n.strings]=End[n.strings]-1
        }
        
        #fill in internal missing records via interpolation if there are any
        if(n.strings>0){
          for(istr in 1:n.strings){
            begin.pl = Duplicated[Begin[istr]]-1
            end.pl = Duplicated[End[istr]]+1
            Last.coords = Coords[begin.pl,]
            Next.coords = Coords[end.pl,]
            slope.incr = end.pl - begin.pl
            Coords[Duplicated[Begin[istr]]:Duplicated[End[istr]],1]=Last.coords[1]+(Next.coords[1]-Last.coords[1])/slope.incr * c(1:(end.pl-begin.pl-1))
            Coords[Duplicated[Begin[istr]]:Duplicated[End[istr]],2]=Last.coords[2]+(Next.coords[2]-Last.coords[2])/slope.incr * c(1:(end.pl-begin.pl-1))
          }
        }
      }
    }
    Interp_loc[Which_photos,]=Coords
  }
  if(n_photos>=10 & iid!=8237){
    Coords = coordinates(boss_geo[Which_photos,])
    I_dup = (Coords[2:n_photos,1]==Coords[1:(n_photos-1),1] & Coords[2:n_photos,2]==Coords[1:(n_photos-1),2])
    Duplicated = as.numeric(which(I_dup==1))+1
    Cur_df=data.frame(img_dt=as.numeric(difftime(boss_geo[Which_photos,]$img_dt,boss_geo[Which_photos[1],]$img_dt,units="secs")))
    Cur_df$x = Coords[,"coords.x1"]
    Cur_df$y = Coords[,"coords.x2"]
    if(length(Duplicated)>0)Cur_df=Cur_df[-Duplicated,]
    Duplicated = which(duplicated(Cur_df$img_dt)==TRUE)
    if(length(Duplicated)>0)Cur_df=Cur_df[-Duplicated,]
    Inits = list(a1.x = c(Coords[1,1],0), a1.y =c(Coords[1,2],0),P1.x=diag(c(10000,10000)),P1.y=diag(c(10000,10000)))
    #now run through crawl
    Fix=c(log(100),log(100),NA,NA)
    #if(iid %in% c(2133,2176,2196,2364,2573,2595,2876,5160,5898,6180,6241,6329,6634,6651,6678,6878,6903,7054,7179))Fix=c(log(150),log(150),NA,NA)
    my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1,y=~1),data=Cur_df,coord=c("x","y"),Time.name="img_dt",need.hess=FALSE,initial.state=Inits,polar.coord=FALSE,fixPar=Fix)
    if(is.character(my_fit)==TRUE){
      Fix=c(log(150),log(150),NA,NA)  
      my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1,y=~1),data=Cur_df,coord=c("x","y"),Time.name="img_dt",need.hess=FALSE,initial.state=Inits,polar.coord=FALSE,fixPar=Fix)
    }
    elapsed = Cur_df[nrow(Cur_df),"img_dt"]-Cur_df[1,"img_dt"]
    t_diff=elapsed/(n_photos-1)
    predTime=as.numeric(Cur_df[1,"img_dt"]+c(0:(n_photos-1))*t_diff)
    Preds=crwPredict(my_fit,predTime=predTime)
    Which_pred = which(Preds$locType=="p")
    Interp_loc[Which_photos,]=cbind(Preds$mu.x[Which_pred],Preds$mu.y[Which_pred])
  }
}
#formulate new sp points object with interpolated locations
boss_geo=SpatialPointsDataFrame(coords=Interp_loc,data=boss_geo@data,proj4str=CRS(laea_180_proj))
#fill in missing altitudes
Which.missing = which(is.na(boss_geo$interp_alt))  
for(i in 1:length(Which.missing))boss_geo$interp_alt[Which.missing[i]]=boss_geo$interp_alt[Which.missing[i]-1]
save(boss_geo,FlightCellSide_ID,FlightCell_ID,FlightCellSide_unique_IDs,file='boss_geo_sp_tmp4.Rda')

load('boss_geo_sp_tmp4.Rda')


#formulate lookup table that gives camera angles, offset by aircraft, lens, and camera position
Angles = expand.grid(aircraft=c("Aero","Otter"),camera=c("Starboard","Port","Center"),lens=c(85,100))
Angles$vert = Angles$horiz = Angles$offset = rep(0,nrow(Angles)) 
Angles[which(Angles$lens==85),"vert"]=16.07  #angles based on Mike's "field of view calcs" spreadshets
Angles[which(Angles$lens==85),"horiz"]=23.85
Angles[which(Angles$lens==100),"vert"]=13.69
Angles[which(Angles$lens==100 & Angles$aircraft=="Otter"),"horiz"]=20.41
Angles[which(Angles$lens==100 & Angles$aircraft=="Aero"),"horiz"]=20.35
Angles[which(Angles$aircraft=="Aero" & Angles$camera=="Port"),"offset"]=13  #offsets based on data from Erin
Angles[which(Angles$aircraft=="Otter" & Angles$camera=="Port"),"offset"]=26
Angles[which(Angles$aircraft=="Aero" & Angles$camera=="Starboard"),"offset"]=13
Angles[which(Angles$aircraft=="Otter" & Angles$camera=="Starboard"),"offset"]=26  
Angles[,c("vert","horiz","offset")]=Angles[,c("vert","horiz","offset")]/360*2*pi  #covert to radians



get_footprint_corners_port <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_far = alt/cos(Fl_angles[,"offset"]+0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_near = alt/cos(Fl_angles[,"offset"]-0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  #horiz_width = alt*tan(0.5*Fl_angles[,"horiz"])
  horiz_dist_far = alt*tan(0.5*Fl_angles[,"horiz"]+Fl_angles[,"offset"]+roll)
  horiz_dist_near = alt*tan(Fl_angles[,"offset"]+roll-0.5*Fl_angles[,"offset"])
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(cos(1.5*pi-bearing),-sin(1.5*pi-bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(cos(bearing-1.5*pi),sin(bearing-1.5*pi)) 
  }
  C_far = horiz_dist_far * Tmp + Tmp_coords
  C_near = horiz_dist_near * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_far/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_far - Tmp
  Poly_points[2,] = C_far + Tmp
  Tmp = rep(vert_width_near/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_near + Tmp
  Poly_points[4,] = C_near - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}

get_footprint_corners_starboard <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_far = alt/cos(Fl_angles[,"offset"]+0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_near = alt/cos(Fl_angles[,"offset"]-0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  horiz_dist_far = alt*tan(0.5*Fl_angles[,"horiz"]+Fl_angles[,"offset"]-roll)
  horiz_dist_near = alt*tan(Fl_angles[,"offset"]-roll-0.5*Fl_angles[,"offset"])
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(cos(0.5*pi-bearing),-sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(cos(bearing-0.5*pi),sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(-cos(1.5*pi - bearing),sin(1.5*pi - bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(-cos(bearing-1.5*pi),-sin(bearing-1.5*pi)) 
  }
  C_far = horiz_dist_far * Tmp + Tmp_coords
  C_near = horiz_dist_near * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_far/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_far - Tmp
  Poly_points[2,] = C_far + Tmp
  Tmp = rep(vert_width_near/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_near + Tmp
  Poly_points[4,] = C_near - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}

get_footprint_corners_center <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_port = alt/cos(0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_starboard = alt/cos(0.5*Fl_angles[,"horiz"]-roll)*sin(0.5*Fl_angles[,"vert"])
  horiz_dist_port = alt*tan(0.5*Fl_angles[,"horiz"]+roll)
  horiz_dist_starboard = alt*tan(0.5*Fl_angles[,"horiz"]-roll)
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(cos(1.5*pi-bearing),-sin(1.5*pi-bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(cos(bearing-1.5*pi),sin(bearing-1.5*pi)) 
  }
  C_port = horiz_dist_port * Tmp + Tmp_coords
  C_starboard = -horiz_dist_starboard * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_port/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_port - Tmp
  Poly_points[2,] = C_port + Tmp
  Tmp = rep(vert_width_starboard/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_starboard + Tmp
  Poly_points[4,] = C_starboard - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}


#plot(Tmp_coords[1],Tmp_coords[2],xlim=c(Tmp_coords[1]-500,Tmp_coords[1]+500),ylim=c(Tmp_coords[2]-500,Tmp_coords[2]+500))
#points(C_far[1],C_far[2],col='orange')
#points(C_near[1],C_near[2],col='orange')
#points(C_port[1],C_port[2],col='cyan')
#points(C_starboard[1],C_starboard[2],col='cyan')
#points(Poly_points[1,1],Poly_points[1,2],col="red")
#points(Poly_points[2,1],Poly_points[2,2],col="purple")
#points(Poly_points[3,1],Poly_points[3,2],col="blue")
#points(Poly_points[4,1],Poly_points[4,2],col="green")

#  Formulate a list of times within flights for which banking turns are happening (for filter)
# base banking inference on camera for unique flight & cell combo that has the most photographs
# step 1: calculate change in bearing
FlightCell_unique_IDs = unique(FlightCell_ID) 
Photos = matrix(0,length(FlightCell_unique_IDs),3)
Delta_bearing=0
Delta_flight="dummy"
Delta_dt=boss_geo[1,]$img_dt
for(iid in 1:length(FlightCell_unique_IDs)){
  Which_photos = which(FlightCell_ID==FlightCell_unique_IDs[iid])
  n_photos=length(Which_photos)
  aircraft=as.character(aircraft_func(boss_geo[Which_photos[1],]$flightid))
  n_cameras=2
  if(aircraft=="Otter")n_cameras=3
  N_photos=rep(0,3)
  N_photos[1]=length(which(boss_geo[Which_photos,]$side=="Port"))
  N_photos[2]=length(which(boss_geo[Which_photos,]$side=="Starboard"))
  if(aircraft=="Otter")N_photos[3]=length(which(boss_geo[Which_photos,]$side=="Center"))
  Photos[iid,]=N_photos
  max_n = max(N_photos)
  if(Photos[iid,3]==max_n)cam_name="Center"  
  if(Photos[iid,1]==max_n)cam_name="Port"
  if(Photos[iid,2]==max_n)cam_name="Starboard"
  if(max_n>3){
    Which_cam = which(boss_geo[Which_photos,]$side==cam_name)
    Coords = Interp_loc[Which_photos[Which_cam],]
    Bearing=rep(0,max_n)
    Bearing = as.numeric(atan2(Coords[2:max_n,2]-Coords[1:(max_n-1),2],Coords[2:max_n,1]-Coords[1:(max_n-1),1]))
    Bearing[max_n]=Bearing[max_n-1]
    Bearing[2:(max_n-1)] = 0.5*(Bearing[1:(max_n-2)]+Bearing[2:(max_n-1)])
    #if(sum(Bearing<0)>0){  #make sure bearing is in [0,2*pi]
    #  Bearing[which(Bearing<0)]=Bearing[which(Bearing<0)]+2*pi
    #}
    Delta_dt = c(Delta_dt,boss_geo[Which_photos[Which_cam],]$img_dt[1:(max_n-1)])
    Delta_flight=c(Delta_flight,rep(boss_geo[Which_photos[1],]$flightid,max_n-1))
    Delta_bearing=c(Delta_bearing,Bearing[2:max_n]-Bearing[1:(max_n-1)])
  }
}
Delta_bearing=Delta_bearing[-1]
Delta_flight=Delta_flight[-1]
Delta_dt = Delta_dt[-1]

#filter based on 10 values with abs(Delta_bearing)
#first make a pass through to document sequence
I_filter=Seq=rep(0,length(Delta_bearing))
Abs_bearing=abs(Delta_bearing)
for(i in 2:length(Delta_bearing)){
  if(Abs_bearing[i]>0.01){
    Seq[i]=1
    if(Seq[i-1]>0 & abs(as.numeric(difftime(Delta_dt[i],Delta_dt[i-1],units="secs")))<5)Seq[i]=Seq[i-1]+1
  }
}
Which.eq.10 = which(Seq==10)  #566 turns
n_turns=length(Which.eq.10)
Filter_df=data.frame(flight=rep("a",n_turns),start=rep(Delta_dt[1],n_turns),end=rep(Delta_dt[1],n_turns),stringsAsFactors=FALSE)
for(ipeak in 1:n_turns){
  cur.pl=Which.eq.10[ipeak]
  I_filter[(cur.pl-9):cur.pl]=1
  Filter_df[ipeak,"flight"]=Delta_flight[cur.pl]
  Filter_df[ipeak,"start"]=Delta_dt[cur.pl-9]
  Filter_df[ipeak,"end"]=Delta_dt[cur.pl]
  if(Seq[cur.pl+1]==11){
    flag=0
    cur.pl2=1
    while(flag==0){
      if(Seq[cur.pl+cur.pl2]>Seq[cur.pl+cur.pl2-1])I_filter[cur.pl+cur.pl2]=1
      else flag=1
      cur.pl2=cur.pl2+1
    }
    Filter_df[ipeak,"end"]=Delta_dt[cur.pl+cur.pl2-2]
  }
}
save(Filter_df,file="BOSS_turn_filter.Rda")
  

#filter out all photos that occur within turns
I_filter=rep(0,length(boss_geo))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(boss_geo$flightid==Filter_df[irow,"flight"] & boss_geo$img_dt>=Filter_df[irow,"start"] & boss_geo$img_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}

boss_geo=boss_geo[which(I_filter==0),]
save(boss_geo,file='boss_geo_sp_tmp5.Rda')

load('boss_geo_sp_tmp5.Rda')
#####  Now, produce footprint given flight track
Alt = .3048*boss_geo$interp_alt  #convert to meters
attr(boss_geo$img_dt, "tzone") <- "UTC"  #convert time to UTC
Sides=c("Port","Starboard","Center")
FlightCell_ID = paste(boss_geo$flightid,boss_geo$GridID)
FlightCell_unique_IDs=unique(FlightCell_ID)
Photos = matrix(0,length(FlightCell_unique_IDs),3)
Photo_area = Flight_ID=Cell_ID=rep(0,length(FlightCell_unique_IDs))
for(iid in 1:length(FlightCell_unique_IDs)){
  Which_photos = which(FlightCell_ID==FlightCell_unique_IDs[iid])
  n_photos=length(Which_photos)
  aircraft=as.character(aircraft_func(boss_geo[Which_photos[1],]$flightid))
  n_cameras=2
  if(aircraft=="Otter")n_cameras=3
  N_photos=rep(0,3)
  N_photos[1]=length(which(boss_geo[Which_photos,]$side=="Port"))
  N_photos[2]=length(which(boss_geo[Which_photos,]$side=="Starboard"))
  if(aircraft=="Otter")N_photos[3]=length(which(boss_geo[Which_photos,]$side=="Center"))
  Photos[iid,]=N_photos
  lens = 100
  flightid=boss_geo[Which_photos[1],]$flightid
  Flight_ID[iid]=flightid
  Cell_ID[iid]=boss_geo[Which_photos[1],]$GridID
  Polys=vector("list",n_photos)
  cur_pl = 0
  for(iside in 1:n_cameras){
    side=Sides[iside]
    if(aircraft=="Aero"){
      tmp_time = as.POSIXct("2012-04-13 20:01:00",tz="UTC")
      if(flightid=="12_AeroFl04" & boss_geo[Which_photos,]$img_dt[1]>tmp_time)lens=85
      if(flightid=="12_AeroFl09")lens=85
      if(flightid=="13_AeroFl06" & side=="Port")lens=85
      if(flightid=="13_AeroFl07" & side=="Port")lens=85
      if(flightid=="13_AeroFl08" & side=="Port")lens=85
      if(flightid=="13_AeroFl09" & side=="Port")lens=85
      if(flightid=="13_AeroFl10")lens=85
      if(flightid=="13_AeroFl11")lens=85
    }
    if(iside>1)cur_pl=cur_pl+N_photos[iside-1]
    Which_cam = which(boss_geo[Which_photos,]$side==Sides[iside])
    if(N_photos[iside]>0){
      Coords = boss_geo[Which_photos[Which_cam],]@coords
      Bearing=rep(0,N_photos[iside])
      if(N_photos[iside]>1){
        Bearing = as.numeric(atan2(Coords[2:N_photos[iside],2]-Coords[1:(N_photos[iside]-1),2],Coords[2:N_photos[iside],1]-Coords[1:(N_photos[iside]-1),1]))
        Bearing[N_photos[iside]]=Bearing[N_photos[iside]-1]
        Bearing[2:(N_photos[iside]-1)] = 0.5*(Bearing[1:(N_photos[iside]-2)]+Bearing[2:(N_photos[iside]-1)])
      }
      if(sum(Bearing<0)>0){  #make sure bearing is in [0,2*pi]
        Bearing[which(Bearing<0)]=Bearing[which(Bearing<0)]+2*pi
      }
      if(side=="Port")footprint_fn = get_footprint_corners_port
      if(side=="Center")footprint_fn = get_footprint_corners_center
      if(side=="Starboard")footprint_fn = get_footprint_corners_starboard
      Tmp_angles = as.matrix(Angles[which(Angles[,"aircraft"]==aircraft & Angles[,"camera"]==side & Angles[,"lens"]==lens),c("vert","horiz","offset")])
      #now determine footprints
      for(iphoto in 1:N_photos[iside]){
        Polys[[cur_pl+iphoto]]=Polygons(list(Polygon(footprint_fn(Fl_angles=Tmp_angles,roll=0,bearing=Bearing[iphoto],Tmp_coords=Coords[iphoto,],alt=Alt[Which_photos[Which_cam[iphoto]]]))),paste(cur_pl+iphoto))
      }
    }
  }
  SPDF = SpatialPolygons(Polys,proj4string=CRS(laea_180_proj))  
  if(iid==10)save(SPDF,file="12_AeroFl01_701_spdf.Rda")
  Union = gUnionCascaded(SPDF)
  Photo_area[iid]=gArea(Union)
}  

Area_table = data.frame(flightid=Flight_ID,Grid_ID=Cell_ID,Area_m2 = Photo_area)
save(Area_table,file="Area_photographed.Rda")
  
