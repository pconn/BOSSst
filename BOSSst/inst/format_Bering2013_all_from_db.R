#format data for BOSS spatio-temporal analysis
load('./data_from_JML/boss_hotspots_sp.Rda') #hotspots
#load('./data_from_JML/boss_geo_sp.Rda') #on effort points
load('./data_from_JML/boss_grid_env.Rda') # effort summary (without area surveyed)
load('./AlaskaBeringData2012_2013_14Dec2015.Rdat')  #boss grid w covariates by day
load('./Area_photographed.rda')  #produced with 'Calculate_area_surveyed_crawl.R'
load('BOSS_Flt_dates.Rda')
load('Knot_cell_distances.Rdata') #load object giving K matrix, Q for knots
load('p13_obseffect.RData')  #read in confusion array
load('Haulout_samples.Rdat')  #read in haulout proportion MCMC samples

library(rgeos)
library(sp)
library(RPostgreSQL)
library(sf)

set.seed(12345)  #need because pseudo-zeroes are selected randomly

# Run code -------------------------------------------------------
# Extract data from DB ------------------------------------------------------------------
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0","+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
                       
#boss <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM surv_boss.hotspot")

boss.sf <- sf::st_read_db(con, 
                          query = "SELECT * FROM surv_boss.hotspot", 
                          geom_column = "geom")
boss.sf2 = boss.sf[which(boss.sf$hotspot_type=="seal"),]
boss.sf2 = sf::st_transform(boss.sf2,laea_180_proj)

#remove duplicated entries - these are mom-pup pairs  #####  Stacie fixed in DB
# Dup = which(duplicated(boss.sf2$hotspot_id))
# for(i in 1:length(Dup)){
#   Cur.which = which(boss.sf2$hotspot_id==boss.sf2$hotspot_id[Dup[i]])
#   boss.sf2$gross_age[Cur.which[1]]="mom_pup"
# }
# boss.sf2=boss.sf2[-Dup,]
                       
#restrict to 'on effort'
boss.sf2=boss.sf2[which(boss.sf2$effort=="On"),]

#limit to 2012
boss.sf2=boss.sf2[which(lubridate::year(boss.sf2$hotspot_dt)==2013),]

boss_grid_env = data.frame(boss_grid_env)
#reproject 
#grid_spdf=spTransform(grid_spdf, CRS(laea_180_proj))

#boss_hotspots_sp=spTransform(boss_hotspots_sp, CRS(laea_180_proj))

#cur_fl="12_OtterFl10"
#crap=which(Area_table[,"flightid"]==cur_fl)
#length(unique(Area_table[crap,"Grid_ID"]))
#crap=which(boss_grid_env[,"flightid"]==cur_fl)
#length(crap)


# 2013: Run analysis for 4/7 - 5/4
date.start=as.Date("2013-04-07")
date.end=as.Date("2013-05-04")
t.steps=as.numeric(date.end-date.start)+1
Data$Grid = Data$Grid$y2013
Day_PST = as.numeric(as.Date(boss_grid_env$grid_dt,tz="PST8PDT")-date.start)
Which_dates = which(Day_PST>=0 & Day_PST<t.steps)
Area_table[,"flightid"]=as.character(Area_table[,"flightid"])
#Area_yr = matrix(unlist(strsplit(Area_table[,"flightid"],'_')),nrow=2)[1,]
#Area_table=Area_table[which(Area_yr=="12"),] #limit to 2012
Flts_include=Flt_table[which(Flt_table[,"Dates"]>=date.start & Flt_table[,"Dates"]<=date.end),"Flt"]
Flts_include = Flts_include[-which(Flts_include %in% c("13_AeroFl19","13_OtterFl14","13_OtterFl15"))]
Area_table=Area_table[which(Area_table$flightid%in%Flts_include),]
boss_grid_env=boss_grid_env[which(boss_grid_env$flightid %in% Flts_include),]
Day_PST = as.numeric(as.Date(boss_grid_env$grid_dt,tz="PST8PDT")-date.start)


Cell_IDs=as.numeric(rownames(Data$Grid[[1]]@data))
S=nrow(Data$Grid[[1]])
for(it in 1:t.steps){ #rename cell IDs so they go 1:S
  for(ipoly in 1:S){
    Data$Grid[[it]]@polygons[[ipoly]]@ID=as.character(ipoly)
  }
}


#remove effort when altitude >1250ft
Which_remove = which(boss_grid_env[,"interp_alt"]>1250)
Flt_remove= boss_grid_env[Which_remove,"flightid"]
Cell_remove=boss_grid_env[Which_remove,"objectid"]
for(ifl in 1:length(Which_remove))Area_table=Area_table[-which(Area_table[,"Grid_ID"]==Cell_remove[ifl] & Area_table[,"flightid"]==Flt_remove[ifl]),]
if(length(Which_remove)>0){
  boss_grid_env=boss_grid_env[-Which_remove,]
  Day_PST = Day_PST[-Which_remove]
}

#remove cells not in the survey grid
Which_remove = which(!(boss_grid_env[,"objectid"]%in%Cell_IDs))
Flt_remove= boss_grid_env[Which_remove,"flightid"]
Cell_remove=boss_grid_env[Which_remove,"objectid"]
for(ifl in 1:length(Which_remove))Area_table=Area_table[-which(Area_table[,"Grid_ID"]==Cell_remove[ifl] & Area_table[,"flightid"]==Flt_remove[ifl]),]
if(length(Which_remove)>0){
  boss_grid_env=boss_grid_env[-Which_remove,]
  Day_PST = Day_PST[-Which_remove]
}


comb_row <- function(Rows){
  New_row=Rows[1,1:(length(Rows)-1)]
  New_row["grid_dt"]=mean(Rows[,"grid_dt"])
  New_row[c("narr_air2m","narr_prmsl","narr_uwnd","narr_vwnd")]=apply(Rows[,c("narr_air2m","narr_prmsl","narr_uwnd","narr_vwnd")],2,'mean')
  New_row[c("img_count_port","img_count_center","img_count_starboard")]=apply(Rows[,c("img_count_port","img_count_center","img_count_starboard")],2,'sum')
  Tmp=Rows$seal_id_array[[1]]
  for(irow in 2:nrow(Rows))Tmp=rbind(Tmp,Rows$seal_id_array[[irow]])
  New_row$seal_id_array=vector("list",1)
  New_row$seal_id_array[[1]]=Tmp
  New_row
}

#combine data for cells sampled >1 time in a given day
Day_Cell=paste(Day_PST,boss_grid_env[,"objectid"])
Duplicated=which(duplicated(Day_Cell))
cur_pl=1
Firsts=Firsts_area=Duplicated_area=rep(0,length(Duplicated))
Area_table$Grid_ID = as.numeric(as.character(Area_table$Grid_ID))
Replace=boss_grid_env[1:length(Duplicated),]
Replace_area=rep(0,length(Duplicated))
for(i in 1:length(Duplicated)){
  Firsts[i]=which(Day_Cell==Day_Cell[Duplicated[i]])[1]
  Replace[i,]=comb_row(rbind(boss_grid_env[Firsts[i],],boss_grid_env[Duplicated[i],]))
  Firsts_area[i] = which(Area_table$Grid_ID==boss_grid_env[Firsts[i],"objectid"] & Area_table$flightid==boss_grid_env[Firsts[i],"flightid"])
  Duplicated_area[i] = which(Area_table$Grid_ID==boss_grid_env[Duplicated[i],"objectid"] & Area_table$flightid==boss_grid_env[Duplicated[i],"flightid"])
}

boss_grid_env=boss_grid_env[-Duplicated,]
Row_index = Firsts - c(0:(length(Firsts)-1))
boss_grid_env[Row_index,]=Replace
#ONLY WORKS BECAUSE DUPLICATED ALL OCCUR AFTER FIRSTS in Area_table
Area_table[Firsts_area,"Area_m2"]=Area_table[Firsts_area,"Area_m2"]+Area_table[Duplicated_area,"Area_m2"]  
Area_table=Area_table[-Duplicated_area,]
Day_PST = Day_PST[-Duplicated]

n_surveyed=nrow(boss_grid_env)

#Produce Day, Hour in solar hour for haulout predictions
Longitude.grid = coordinates(spTransform(rgeos::gCentroid(Data$Grid[[1]],byid=TRUE), CRS("+proj=longlat +datum=WGS84")))[,1]
Grid_sampled = rep(0,n_surveyed)  #grid cell # for each row in boss_grid_env
for(i in 1:n_surveyed)Grid_sampled[i]=which(Cell_IDs==boss_grid_env$objectid[i])
Longitude.sampled = Longitude.grid[Grid_sampled]
SolarT = solaR::local2Solar(boss_grid_env$grid_dt,lon=Longitude.sampled)
DayHour=matrix(0,n_surveyed,2)
DayHour[,1] = as.numeric(as.Date(SolarT)-date.start+1)
DayHour[,2] = as.numeric(strftime(SolarT, format="%H"))
#this makes Hour = 0 ; need it to go 1:24  (no longer needed since not in UTC)
#for(i in 1:n_surveyed){
#  if(DayHour[i,2]==0)DayHour[i,]=c(DayHour[i,1]-1,24)
#}

#translate Area_table by new cell IDs
Grid_ID_new = rep(0,nrow(Area_table))
for(i in 1:nrow(Area_table))Grid_ID_new[i]=which(Cell_IDs==Area_table[i,"Grid_ID"])
Area_table[,2]=Grid_ID_new
#strcat <- function(x)paste(x[1],x[2])
#Flight_grid_id=apply(Area_table,1,"strcat")
#length(unique(Flight_grid_id))  #confirmed that there is one record per grid cell per flight

#assemble "seal_data_array" directly from database; note Josh originally did this but we need to replace values
# since database has been updated
#get a cell # for each hotspot
hotspots.sp = as(as(boss.sf2,"Spatial"),'SpatialPoints')
Grid.sp = as(Data$Grid[[1]],"SpatialPolygons")
Cell.hs = over(hotspots.sp,Grid.sp,fn=NULL)
#boss.sf2 = boss.sf2[-which(is.na(Cell.hs)),]
#Cell.hs = Cell.hs[-which(is.na(Cell.hs))]
SolarT = solaR::local2Solar(boss.sf2$gps_dt,lon=boss.sf2$longitude)
Day.hs = as.numeric(as.Date(SolarT)-date.start+1)
Sampled = rep(0,length(Cell.hs))
for(icell in 1:n_surveyed){
  Cells.surv = which(Cell_IDs==boss_grid_env$objectid[icell])
  Which_hotspots=which(Cell.hs %in% Cells.surv & Day.hs==DayHour[icell,1])
  if(length(Which_hotspots)>0){
    Sampled[Which_hotspots]=1
    #cat(paste(icell,'\n'))
    Cur.seals = boss.sf2[Which_hotspots,c("hotspot_id","hotspot_type","num_seals","species","species_user","species_conf")] 
    boss_grid_env$seal_id_array[[icell]]=Cur.seals
  }
  else boss_grid_env$seal_id_array[[icell]]=list()
}

#determine how many hot spot records there are for each cell surveyed, record grid cell surveyed and day, hour surveyed
N_records = rep(0,n_surveyed)
Mapping = matrix(0,n_surveyed,2)
for(i in 1:n_surveyed){
  Mapping[i,1]=which(Cell_IDs==boss_grid_env[i,"objectid"])
  if(length(boss_grid_env[i,"seal_id_array"][[1]])>0)N_records[i]=nrow(boss_grid_env[i,"seal_id_array"][[1]])
}
Mapping[,2]=Day_PST+1
n_records=sum(N_records)


#form matrix of observations
counter=1
Obs_mat = data.frame(matrix(0,n_records,7))
for(i in 1:n_surveyed){
  if(N_records[i]>0){
    Obs_mat[counter:(counter+N_records[i]-1),]=cbind(rep(i,N_records[i]),data.frame(boss_grid_env[i,"seal_id_array"][[1]])[,1:6])
    counter=counter+N_records[i]
  }
}
colnames(Obs_mat)=c("id_row","hotspotid","hotspot_type","numseals","species","species_user","species_conf")
if(length(which(Obs_mat[,"hotspot_type"]=="other_animal"))>0)Obs_mat=Obs_mat[-which(Obs_mat[,"hotspot_type"]=="other_animal"),]

#check to make sure no hotspot IDs are duplicated
if(length(unique(Obs_mat[,"hotspotid"])) != nrow(Obs_mat))cat("WARNING!  some hotspot IDs duplicated")

#construct count data set
#Obs_mat2=Obs_mat2[-2364,]  #no seal in image
Count_data_BOSS = data.frame(matrix(1,nrow(Obs_mat),4)) 
Count_data_BOSS[,1]=Obs_mat[,"id_row"]
Count_data_BOSS[,4]=Obs_mat[,"numseals"]
n_hotspots = nrow(Obs_mat)
for(irow in 1:n_hotspots){
  #if(irow %in% c(1484))Obs_mat[irow,3]=13 #hot spots missing species info that were unknown
  #else{
    if(Obs_mat[irow,"species"]=="sd" & Obs_mat[irow,"species_conf"]=="guess")Count_data_BOSS[irow,3]=1
    if(Obs_mat[irow,"species"]=="sd" & Obs_mat[irow,"species_conf"]=="likely")Count_data_BOSS[irow,3]=2
    if(Obs_mat[irow,"species"]=="sd" & Obs_mat[irow,"species_conf"]=="pos")Count_data_BOSS[irow,3]=3
    if(Obs_mat[irow,"species"]=="rn" & Obs_mat[irow,"species_conf"]=="guess")Count_data_BOSS[irow,3]=4
    if(Obs_mat[irow,"species"]=="rn" & Obs_mat[irow,"species_conf"]=="likely")Count_data_BOSS[irow,3]=5
    if(Obs_mat[irow,"species"]=="rn" & Obs_mat[irow,"species_conf"]=="pos")Count_data_BOSS[irow,3]=6
    if(Obs_mat[irow,"species"]=="bd" & Obs_mat[irow,"species_conf"]=="guess")Count_data_BOSS[irow,3]=7
    if(Obs_mat[irow,"species"]=="bd" & Obs_mat[irow,"species_conf"]=="likely")Count_data_BOSS[irow,3]=8
    if(Obs_mat[irow,"species"]=="bd" & Obs_mat[irow,"species_conf"]=="pos")Count_data_BOSS[irow,3]=9
    if(Obs_mat[irow,"species"]=="rd" & Obs_mat[irow,"species_conf"]=="guess")Count_data_BOSS[irow,3]=10
    if(Obs_mat[irow,"species"]=="rd" & Obs_mat[irow,"species_conf"]=="likely")Count_data_BOSS[irow,3]=11
    if(Obs_mat[irow,"species"]=="rd" & Obs_mat[irow,"species_conf"]=="pos")Count_data_BOSS[irow,3]=12
    if(Obs_mat[irow,"species"]=="unk")Count_data_BOSS[irow,3]=13
  #}
}
colnames(Count_data_BOSS)=c("Transect","Photo","Obs","Group")

#calculate area surveyed for each grid cell and time
Flt_day = as.numeric(Flt_table[,"Dates"]-date.start)+1
Area = rep(NA,nrow(Mapping))
for(imap in 1:nrow(Mapping)){
  Which_flts = Flt_table[which(Flt_day==Mapping[imap,2]),"Flt"]
  Cur_which=which(Area_table[,"Grid_ID"]==Mapping[imap,1] & Area_table[,"flightid"]%in%Which_flts)
  Area[imap]=sum(Area_table[Cur_which,"Area_m2"])
}
Area=Area/gArea(Data$Grid[[1]][1,])


#Knot calculations
Coords=coordinates(Data$Grid[[1]])
x.min=min(Coords[,1])-100000
x.max=max(Coords[,1])+100000
y.min=min(Coords[,2])-100000
y.max=max(Coords[,2])+100000

X=x.min+(x.max-x.min)/8*c(0:8)
Y=y.min+(y.max-y.min)/8*c(0:8)
XY=expand.grid(x=X,y=Y)

Knots=SpatialPoints(coords=XY,proj4string=CRS(proj4string(Data$Grid[[1]])))
#save(Knots,file="BOSS_Knots_SP.Rda")

Distances=gDistance(Knots,Data$Grid[[1]],byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=150000
Which.include=which(Distances<my.buffer)
Knot.cell.distances=gDistance(Knots[Which.include,],Data$Grid[[1]],byid=TRUE)
diff.x=(x.max-x.min)/6
diff.y=(y.max-y.min)/6
sigma=(diff.x+diff.y)/2

K=dnorm(Knot.cell.distances,0,sigma)
K=K/apply(K,1,'sum')
save(K,file="Knot_cell_distances_Apr2018.Rdata")  #note these implicitly include assumed sigma value



##more for BOSS 2013 analysis
Dat=Count_data_BOSS
Area_hab = 1-Data$Grid[[1]]$land_cover
Area_trans = Area

#fomulate sample sizes on each camera (used to produce beta mixture prior on p owing to different detection probabilities on different cameras)
#note that 2013 entirely used Skeyes 2.0 approach (68/70 detections)
Flt_table$Flt=as.character(Flt_table$Flt)
Flt_table$det_port=0  #if 1, use manual approach (66/70)
Flt_table$det_center=0
Flt_table$det_star=0
Flt_table$I_otter=0
Flt_table$yr = matrix(unlist(strsplit(as.character(Flt_table[,"Flt"]),'_')),nrow=2)[1,]
Flt_table$I_otter[grep("Otter",as.character(Flt_table[,"Flt"]))]=1
Flt_table[which(Flt_table$I_otter==0 & Flt_table$yr=="12"),"det_port"]=1
Flt_table[which(Flt_table$I_otter==1 & Flt_table$yr=="12"),"det_center"]=1
Flt_table[which(Flt_table$Flt=="12_OtterFl12"),c("det_port","det_center","det_star")]=1
Det_wgt = rep(0,n_surveyed)


### produce a spatial points object giving grid cell centroid, date-time, and detection algortihm indicator for each cell surveyed
Centroids_sf = st_as_sf(gCentroid(Data$Grid[[1]],byid=TRUE))
Sampled_cells_ST = Centroids_sf[Grid_sampled,1]
Sampled_cells_ST$dt = boss_grid_env$grid_dt
Sampled_cells_ST$lon = Longitude.sampled
Sampled_cells_ST$Ep = rep(0,n_surveyed)
Sampled_cells_ST$Varp = rep(0,n_surveyed)
Ep1 = 66/70
Ep2 = 68/70
varp1 = Ep1*(1-Ep1)/70
varp2 = Ep2*(1-Ep2)/70
for(isamp in 1:n_surveyed){
  I_method = Flt_table[which(Flt_table$Flt==boss_grid_env[isamp,"flightid"]),c("det_port","det_center","det_star")]
  Wgt = c(sum(as.matrix(I_method,nrow=1)%*%matrix(as.numeric(boss_grid_env[isamp,c("img_count_port","img_count_center","img_count_starboard")]),ncol=1)),
          sum(as.matrix(1-I_method,nrow=1)%*%matrix(as.numeric(boss_grid_env[isamp,c("img_count_port","img_count_center","img_count_starboard")]),ncol=1)))
  Wgt=Wgt/sum(Wgt)    
  Sampled_cells_ST[isamp,"Ep"]=Wgt[1]*Ep1+Wgt[2]*Ep2
  Sampled_cells_ST[isamp,"Varp"]=Wgt[1]*varp1+Wgt[2]*varp2
}
save(Sampled_cells_ST,file="SampledCells_BOSS2013.RData")


# 
# #produce thinning array for Bayesian models - reformulating this to be an n_species * Cells sampled * n_iter array 
# n_species=4
# Thin = array(1,dim=c(n_species,n_surveyed,1000))
# P=rbeta(1000,67,5)  #conjugate beta(1,1) for binomial detection data (66/70 successes)
# for(isamp in 1:n_surveyed){
#   I_method = Flt_table[which(Flt_table$Flt==boss_grid_env[isamp,"flightid"]),c("det_port","det_center","det_star")]
#   Wgt = c(sum(as.matrix(I_method,nrow=1)%*%matrix(as.numeric(boss_grid_env[isamp,c("img_count_port","img_count_center","img_count_starboard")]),ncol=1)),
#           sum(as.matrix(1-I_method,nrow=1)%*%matrix(as.numeric(boss_grid_env[isamp,c("img_count_port","img_count_center","img_count_starboard")]),ncol=1)))
#   Wgt=Wgt/sum(Wgt)    
#   P1=rbeta(1000,67,5) #manual method
#   P2=rbeta(1000,69,3)  #Skeyes 2.0 method
#   I_method=(runif(1000)<Wgt[1])*1
#   P=I_method*P1+(1-I_method)*P2  #two point mixture based on photo sample size
#   for(isp in 1:4){
#     Thin[isp,isamp,]=P
#     ###TEMPORARY FIX: Set day 30 = Day 29; in future, should get haulout data for 1 day after last survey start day
#     if(DayHour[isamp,1]==30)DayHour[isamp,1]=29
#     if(isp==1)Thin[isp,isamp,]=Thin[isp,isamp,]*Haulout.samples$spotted[DayHour[isamp,1],DayHour[isamp,2],]
#     if(isp==2)Thin[isp,isamp,]=Thin[isp,isamp,]*Haulout.samples$ribbon[DayHour[isamp,1],DayHour[isamp,2],]
#     if(isp==3)Thin[isp,isamp,]=Thin[isp,isamp,]*Haulout.samples$bearded[DayHour[isamp,1],DayHour[isamp,2],]
#   }
# }

Prop_photo=rep(1,n_surveyed)

n_hab_col=ncol(Data$Grid[[1]]@data)
Hab_cov = data.frame(matrix(0,S*t.steps,n_hab_col+2))
colnames(Hab_cov)=c(colnames(Data$Grid[[1]]@data),"ice2","depth2")
counter=1
for(it in 1:t.steps){
  Tmp_data = Data$Grid[[it+1]]@data  #first grid entry on Apr 5 but first survey on Apr 6
  Tmp_data$ice2 = Tmp_data$ice_conc^2
  Tmp_data$depth2 = Tmp_data$depth^2
  Hab_cov[counter:(counter+S-1),]=Tmp_data
  counter=counter+S
}
Hab_cov$Ecoregion=factor(Hab_cov$Ecoregion)
Hab_cov$depth=Hab_cov$depth/mean(Hab_cov$depth)  #note: this is the 2nd time this has been standardized (1st was wrt US + Russia grid)
Hab_cov$depth2=Hab_cov$depth^2
Hab_cov$sqrt_edge = sqrt(Hab_cov$dist_edge)


#create Psi matrix as average of observer relative contributions to the dataset
Observer_counts = summary(as.factor(Obs_mat[,"species_user"]))
Observer_prop=as.numeric(Observer_counts[c("KYM.YANO","GAVIN.BRADY","ERIN.MORELAND","ERIN.RICHMOND")])
Observer_prop = Observer_prop/sum(Observer_prop)
Which_psi = sample(c(1:(dim(p13)[4])),10000)
p13_sampled=p13[,,,Which_psi]
Psi=Observer_prop[1]*p13_sampled[,,1,]
for(iobs in 2:4){
  Psi=Psi+Observer_prop[iobs]*p13_sampled[,,iobs,]
}
#Mat_psi = matrix(Psi,4*13,10000)
Mu_psi = apply(Psi,c(1,2),'mean')
#convert to multinomial logit scale
Cur_beta = Cur_beta2 = matrix(0,4,13)
MisID_zero_cols = c(3,6,9,12)
Beta_chain = matrix(0,36,10000)
for(isp in 1:4){
  log_Real0 = log(Mu_psi[isp,MisID_zero_cols[isp]])
  Cur_beta[isp,] = log(Mu_psi[isp,])-log_Real0
}
Beta_psi=as.vector(t(Cur_beta[,c(1,2,4,5,7,8,10,11,13)]))

for(iiter in 1:10000){
  for(isp in 1:4){ 
    log_Real0 = log(Psi[isp,MisID_zero_cols[isp],iiter])
    Cur_beta2[isp,]=log(Psi[isp,,iiter])-log_Real0
  }
  Beta_chain[,iiter]=as.vector(t(Cur_beta2[,c(1,2,4,5,7,8,10,11,13)]))
}

VC_psi = matrix(0,4*9,4*9)
for(i in 1:(4*9)){
  for(j in 1:(4*9)){
    VC_psi[i,j] = cov(Beta_chain[i,],Beta_chain[j,])
  }
}

#add in some 'zero' data to anchor model in places where there's no ice
Surveyed=c(1:nrow(Mapping))  #number index of mapping values where surveys actually occur
n_zeros=100  #number of 'extra zeros' to put in (max is around 3000; to use max set n.zeros=NA)
Temp=which(Data$Grid[[1]]@data[,"ice_conc"]<.001)
Which_no=matrix(1,length(Temp),2)
Which_no[,1]=Temp
for(it in 2:t.steps){
  Temp=which(Data$Grid[[it+5]]@data[,"ice_conc"]<.001)
  Cur_mat=matrix(it,length(Temp),2)
  Cur_mat[,1]=Temp
  Which_no=rbind(Which_no,Cur_mat)
}
if(is.na(n_zeros)==FALSE)Which_no=Which_no[sample(c(1:nrow(Which_no)),n_zeros),]


#get rid of these if they were actually sampled [already in dataset] and readjust sizes of various quantities
Mapping_1d=(Mapping[,2]-1)*S+Mapping[,1]  
Which_no_1d=(Which_no[,2]-1)*S+Which_no[,1]
I_sampled=Which_no_1d%in%Mapping_1d
I_sampled2=Mapping_1d%in%Which_no_1d
#Area.trans[which(I.sampled2==TRUE)]=0.9  #change area sampled in 0 ice cells to 0.9 #no, this messes up initial density estimates
Which_no_1d=Which_no_1d[which(I_sampled==0)]
Which_no=Which_no[which(I_sampled==0),]
#now add in 0 ice cells to list of sampled cells and adjust related quantities
Mapping=rbind(Mapping,Which_no)
Area_trans=c(Area_trans,rep(.9,length(Which_no_1d)))
#n_transects=n_transects+length(Which_no_1d)
Tmp=matrix(1,nrow(Mapping),2)
Tmp[1:nrow(DayHour),]=as.matrix(DayHour)
Tmp[,1]=Mapping[,2]
DayHour=Tmp
#Thin2=array(1,dim=c(n_species,nrow(DayHour),dim(Thin)[3]))
#Thin2[,1:length(Surveyed),]=Thin
#Thin=Thin2  


#look for seals seen where not much ice - cells that have ice<0.01 where seals are seen are replaced with U(0.05,0.3) random variates
Which_seals=unique(Dat[,"Transect"])
Mapping_1d=(Mapping[,2]-1)*S+Mapping[,1]  
I_replace=I_seal=rep(0,length(Mapping_1d))
Which_less=which(Hab_cov[Mapping_1d,"ice_conc"]<0.01)
I_replace[Which_less]=1
I_seal[Which_seals]=1
I_replace=I_seal*I_replace
Which_replace=which(I_replace==1)
for(iobs in 1:length(Which_replace)){ #8 such replacements
  cur_s=Mapping[Which_less[iobs],1]
  cur_t=Mapping[Which_less[iobs],2]
  new_ice=Data$Adj[cur_s,]%*%Data$Grid[[cur_t]][["ice_conc"]]/sum(Data$Adj[cur_s,])
  #new_ice=runif(1,0.05,0.3)
  Hab_cov[S*(cur_t-1)+cur_s,"ice_conc"]=new_ice
  Hab_cov[S*(cur_t-1)+cur_s,"ice2"]=new_ice^2
  cat(paste('replaced cell',cur_s,'time',cur_t,'with original ice_conc',Data$Grid[[cur_t+5]]@data[cur_s,"ice_conc"],'with ice=',new_ice,'\n'))
}


# ###output bearded data set for preferential sampling analysis
# n_bearded=length(which(Obs_mat2[,"species"]=="bd"))
# Obs_bearded=Obs_mat2[which(Obs_mat2[,"species"]=="bd"),]
# Count_bd = c(tabulate(Obs_bearded[,"id_row"]),rep(0,length(Mapping[,1])-max(Obs_bearded$id_row)))
# Count_data_bearded=matrix(0,length(Count_bd),3)
# colnames(Count_data_bearded)=c('Cell','AreaSurveyed','Count')
# Count_data_bearded[,'Count']=Count_bd
# Count_data_bearded[,'Cell']=Mapping[,1]
# Count_data_bearded[,'AreaSurveyed']=c(Area,rep(0,100))
# Bearded_effort = list(Mapping=Mapping,Count.data=Count_data_bearded,Area.hab=1-Data$Grid[[1]]$land_cover,Area.trans=Area,DayHour=DayHour)
# save(Bearded_effort,file="Bearded_effort.Rda")
# 

#take a look at where these occur
#plot_N_map(1,as.matrix(Data$Grid[[20]][["ice_conc"]],ncol=1),highlight=581,Grid=Data$Grid)
#plot_N_map(1,as.matrix(Data$Grid[[15]][["ice_conc"]],ncol=1),highlight=636,Grid=Data$Grid)

Knot_loc = XY[Which.include,]
save(Mapping,Surveyed,Dat,K,Knot_loc,Coords,Area_hab,Area_trans,DayHour,Prop_photo,Hab_cov,VC_psi,Beta_psi,Psi,file="BOSS_data_TMB_2013.Rda")

