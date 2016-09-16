#EDA for habitat covariates and time of year

data(AlaskaBeringData2012_17April2014)  #boss grid, ice data
data(Knot_cell_distances) #load object giving K matrix, Q for knots
load("Effort2012_BOSSst_5Sep2014.Rdata")  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)



#summarize number of each species seen by cell-day
#columns are Cell, Time, Area surveyed, Spotted count, Ribbon ct, Bearded ct, Ringed ct

Counts=data.frame(Effort$Mapping[,1:2])
colnames(Counts)=c("Cell","Time")
Counts$Area=Effort$Area.trans*Effort$Area.hab[Counts$Cell]
Counts$Spotted=0
Counts$Ribbon=0
Counts$Bearded=0
Counts$Ringed=0
for(i in 1:nrow(Counts)){
  Which.entries=which(Effort$Count.data[,"Transect"]==i)
  if(length(Which.entries)>0){
    for(j in 1:length(Which.entries)){
      if(Effort$Count.data[Which.entries[j],"Obs"]%in%c(1:3))Counts[i,"Spotted"]=Counts[i,"Spotted"]+1
      if(Effort$Count.data[Which.entries[j],"Obs"]%in%c(4:6))Counts[i,"Ribbon"]=Counts[i,"Ribbon"]+1
      if(Effort$Count.data[Which.entries[j],"Obs"]%in%c(7:9))Counts[i,"Bearded"]=Counts[i,"Bearded"]+1
      if(Effort$Count.data[Which.entries[j],"Obs"]%in%c(10:12))Counts[i,"Ringed"]=Counts[i,"Ringed"]+1        
    }
  }
}


#attach habitat covariates (ice, etc.)
Temp=data.frame(matrix(0,nrow(Effort$Mapping),ncol(Data$Grid[[1]]@data)))
colnames(Temp)=colnames(Data$Grid[[1]]@data)
for(i in 1:nrow(Counts))Temp[i,]=Data$Grid[[Counts$Time[i]]]@data[Counts$Cell[i],]
Counts=cbind(Counts,Temp)


#fit GAM models for each species
library(mgcv)
Ringed.Time=gam(Ringed~s(Time),family=poisson,data=Counts,offset=Area)
plot(Ringed.Time)

Species=c("Spotted","Ribbon","Bearded","Ringed")
Covs=colnames(Data$Grid[[1]]@data)[2:8]
for(isp in 1:4){
  for(icov in 1:length(Covs)){
    eval(parse(text=paste(Species[isp],'.',Covs[icov],'=gam(',Species[isp],'~s(',Covs[icov],'),family=poisson,data=Counts,offset=Area)',sep='')))
  }
}

plot(Spotted.dist_mainland)  #increasing, with a wiggle
plot(Spotted.dist_shelf) #nonlinear - decreases rapidly then increases then flat (quadratic maybe okay)
plot(Spotted.depth) #nonlinear - concave quadratic okay
plot(Spotted.ice_conc) #quadratic ok
plot(Spotted.dist_contour) #decreasing wiht some wiggles
plot(Spotted.dist_edge) #decreasing overall with concave quadratic near ice edge and then some wiggles

plot(Ribbon.dist_mainland)  #nonlinear - slight increase followed by large increase - quadratic ok
plot(Ribbon.dist_shelf) #nonlinear - decreases rapidly then rel flat
plot(Ribbon.depth) #nonlinear - increases rapidly then rel flat
plot(Ribbon.ice_conc) #increase to ~0.2 then rel flat
plot(Ribbon.dist_edge) #concave at small values, then flat or wobbly
plot(Ribbon.dist_contour) #concave at small values, then slight decrease

plot(Bearded.dist_mainland)  #nonlinear - downward trend maybe w concave quadratic thrown in
plot(Bearded.dist_shelf) #nonlinear - slight downward w pronouncedconcave quadratic in middle
plot(Bearded.depth) #linearly decreasing
#Dat.tmp=Counts[-which(Counts$depth>0.4),]
#Bearded.depth=gam(Bearded~s(depth,bs="cr",k=5),family=poisson,offset=Area,data=Dat.tmp)
#plot(Bearded.depth)
plot(Bearded.ice_conc) #typical concave quadratic
plot(Bearded.dist_contour) #bimodal - larger mode farther away
plot(Bearded.dist_edge) #trimodal!

plot(Ringed.dist_mainland)  #bimodal - largest near land, decreasing, increasing again at middle values, then declining
plot(Ringed.dist_shelf) #nonlinear - bimodal
#plot(Ringed.depth) #nonlinear - increases rapidly then rel flat
Dat.tmp=Counts[-which(Counts$depth>0.4),]
Ringed.depth=gam(Ringed~s(depth),family=poisson,offset=Area,data=Dat.tmp)
plot(Ringed.depth) #decreasing with depth, then increasing again
plot(Ringed.ice_conc) #quadratic should do
plot(Ringed.dist_edge) #multimodal - low peak near edge, medium peak at mid-range, highest farthest away from edge



#plot(predict(Ringed.Time))

