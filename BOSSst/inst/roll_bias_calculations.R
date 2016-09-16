#Calculate effective area surveyed based on field of view calculations;
# the approach works by joining the footprints associated with each photograph
library(sp)
library(rgeos)
library(nPacMaps)
library(crawl)


laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")



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
  if(bearing<(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<(1.5*pi)){
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
  if(bearing<(0.5*pi)){
    Tmp = c(cos(0.5*pi-bearing),-sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<pi){
    Tmp = c(cos(bearing-0.5*pi),sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<(1.5*pi)){
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
  if(bearing<(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<(1.5*pi)){
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

#produce table of area covered by aircraft, lens, altitude, roll
Area_df = expand.grid(aircraft=c("Aero","Otter"),lens=c(85,100),altitude=c(800,900,1000,1100),roll=c(0,1,3,5,10))
Area_df=Area_df[-which(Area_df[,"aircraft"]=="Otter" & Area_df[,"lens"]==85),]
Area_df$Area=0
cur_bearing=0
for(irow in 1:nrow(Area_df)){
  cur_roll=Area_df[irow,"roll"]/360*2*pi
  Polys=vector("list",2)
  if(Area_df[irow,"aircraft"]=="Otter")Polys=vector("list",3)
  Pts_star=get_footprint_corners_starboard(Fl_angles=Angles[which(Angles[,"aircraft"]==Area_df[irow,"aircraft"] & Angles[,"lens"]==Area_df[irow,"lens"] & Angles[,"camera"]=="Starboard"),c("vert","horiz","offset")],roll=cur_roll,bearing=cur_bearing,Tmp_coords=c(1000000,-3300000),alt=.3048*Area_df[irow,"altitude"])
  Pts_port=get_footprint_corners_port(Fl_angles=Angles[which(Angles[,"aircraft"]==Area_df[irow,"aircraft"] & Angles[,"lens"]==Area_df[irow,"lens"] & Angles[,"camera"]=="Port"),c("vert","horiz","offset")],roll=cur_roll,bearing=cur_bearing,Tmp_coords=c(1000000,-3300000),alt=.3048*Area_df[irow,"altitude"])
  Polys[[1]]=Polygons(list(Polygon(Pts_star)),"1")
  Polys[[2]]=Polygons(list(Polygon(Pts_port)),"2")
  if(Area_df[irow,"aircraft"]=="Otter"){
    Pts_center=get_footprint_corners_center(Fl_angles=Angles[which(Angles[,"aircraft"]==Area_df[irow,"aircraft"] & Angles[,"lens"]==Area_df[irow,"lens"] & Angles[,"camera"]=="Center"),c("vert","horiz","offset")],roll=cur_roll,bearing=cur_bearing,Tmp_coords=c(1000000,-3300000),alt=.3048*Area_df[irow,"altitude"])
    Polys[[3]]=Polygons(list(Polygon(Pts_center)),"3")
  }
  SP=SpatialPolygons(Polys,proj4string=CRS(laea_180_proj))
  Area_df[irow,"Area"]=gArea(gUnionCascaded(SP))
}

plot_margin=1000
plot(Pts_star,xlim=c(1000000-plot_margin,1000000+plot_margin),ylim=c(-3300000-plot_margin,-3300000+plot_margin))
points(Pts_port,col='blue')
points(Pts_center,col="red")

