### analyze and visualize 'generic' simulation data
library(grid)
library(ggplot2)
library(animation)
library(gridExtra)
t.steps=29

#input effort and grid info for plotting
#load('c:/users/paul.conn/git/STabundance/BOSS_2012Effort_22Apr14.Rdata')  
load('BOSS_data_2012.Rda')
load("AlaskaBeringData2012_2013_14Dec2015.Rdat")  #boss grid, ice data

#load MCMC results from one chain
load('c:/users/paul.conn/git/BOSSst/output/BOSS2012_noST_mcmc_workspace.Rdata')
load('c:/users/paul.conn/git/BOSSst/output/BOSS2012_noST_MCMC_output.Rdata')
load('c:/users/paul.conn/git/BOSSst/Mapping2012_ordered.Rda')  #need ordered Mapping from hierchical_boss_st because Dat from the workspace is ordered 

n.mcmc.iter=nrow(MCMC$MCMC)
n.pred.iter=dim(Post$N)[2]
n.obs=nrow(Dat)

#assemble count dataset
n.sampled = nrow(Mapping)
Counts = data.frame(Ribbon = rep(0,n.sampled),
                    Ringed = rep(0,n.sampled),
                    Spotted = rep(0,n.sampled),
                    Bearded = rep(0,n.sampled),
                    Cell = Mapping[,1],
                    Day = Mapping[,2])
for(irow in 1:n.obs){
  if(Dat[irow,"Species"]==1)Counts[Dat[irow,"Transect"],"Spotted"]=Counts[Dat[irow,"Transect"],"Spotted"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==2)Counts[Dat[irow,"Transect"],"Ribbon"]=Counts[Dat[irow,"Transect"],"Ribbon"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==3)Counts[Dat[irow,"Transect"],"Bearded"]=Counts[Dat[irow,"Transect"],"Bearded"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==4)Counts[Dat[irow,"Transect"],"Ringed"]=Counts[Dat[irow,"Transect"],"Ringed"]+Dat[irow,"Group"]
}



###########
# plot ice covariate, surveys, estimates as a function of 4 different time steps
#########
S = nrow(Data$Grid$y2012[[1]])
require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),legend.title=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
library(RColorBrewer)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))
library(gridExtra)


#create time lapsed video of covariate, transect counts, abundance
Dates=data.frame(matrix(0,29,3))
colnames(Dates)=c("x","y","date")
Dates[,"date"]=c(paste(c(10:30),"April 2012"),paste(c(1:8),"May 2012"))
#spotted animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Spotted"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[1,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,2000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Spotted_Bering_2012.html",outdir=getwd()
)

#ribbon animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ribbon"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[2,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ribbon seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,2000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Ribbon_Bering_2012.html",outdir=getwd()
)

#bearded animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Bearded"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[3,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Bearded seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,600),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Bearded_Bering_2012.html",outdir=getwd()
)

#ringed animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ringed"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[4,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ringed seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Ringed_Bering_2012.html",outdir=getwd()
)


#####   2013 Animations  ############
load('BOSS_data_2013.Rda')
load("AlaskaBeringData2012_2013_14Dec2015.Rdat")  #boss grid, ice data

#load MCMC results from one chain
load('c:/users/paul.conn/git/BOSSst/output/BOSS2013_mcmc_workspace.Rdata')
load('c:/users/paul.conn/git/BOSSst/Mapping2013_ordered.Rda')  #need ordered Mapping from hierchical_boss_st because Dat from the workspace is ordered 

n.mcmc.iter=nrow(MCMC$MCMC)
n.pred.iter=dim(Post$N)[2]
n.obs=nrow(Dat)

#assemble count dataset
n.sampled = nrow(Mapping)
Counts = data.frame(Ribbon = rep(0,n.sampled),
                    Ringed = rep(0,n.sampled),
                    Spotted = rep(0,n.sampled),
                    Bearded = rep(0,n.sampled),
                    Cell = Mapping[,1],
                    Day = Mapping[,2])
for(irow in 1:n.obs){
  if(Dat[irow,"Species"]==1)Counts[Dat[irow,"Transect"],"Spotted"]=Counts[Dat[irow,"Transect"],"Spotted"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==2)Counts[Dat[irow,"Transect"],"Ribbon"]=Counts[Dat[irow,"Transect"],"Ribbon"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==3)Counts[Dat[irow,"Transect"],"Bearded"]=Counts[Dat[irow,"Transect"],"Bearded"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==4)Counts[Dat[irow,"Transect"],"Ringed"]=Counts[Dat[irow,"Transect"],"Ringed"]+Dat[irow,"Group"]
}



###########
# plot ice covariate, surveys, estimates as a function of 4 different time steps
#########
S = nrow(Data$Grid$y2013[[1]])
require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid$y2013[[1]],byid=TRUE))

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),legend.title=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
library(RColorBrewer)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))
library(gridExtra)


#create time lapsed video of covariate, transect counts, abundance
#2013 survey counts go from 4/7-5/6 but environmental covariate data goes 4/5-5/9
Dates=data.frame(matrix(0,30,3))
colnames(Dates)=c("x","y","date")
Dates[,"date"]=c(paste(c(7:30),"April 2013"),paste(c(1:6),"May 2013"))
t.steps=30
#spotted animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/9 but analysis goes 4/7-5/6
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Spotted"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[1,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1600),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Spotted_Bering_2013.html",outdir=getwd()
)

#ribbon animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ribbon"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[2,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ribbon seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1500),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Ribbon_Bering_2013.html",outdir=getwd()
)

#bearded animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Bearded"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[3,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Bearded seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,800),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Bearded_Bering_2013.html",outdir=getwd()
)

#ringed animation
saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ringed"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[4,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ringed seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,500),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,htmlfile="Ringed_Bering_2013.html",outdir=getwd()
)


####redo all as MP4
library(grid)
library(ggplot2)
library(animation)
library(gridExtra)
t.steps=29

#input effort and grid info for plotting
#load('c:/users/paul.conn/git/STabundance/BOSS_2012Effort_22Apr14.Rdata')  
load('BOSS_data_2012.Rda')
load("AlaskaBeringData2012_2013_14Dec2015.Rdat")  #boss grid, ice data

#load MCMC results from one chain
load('c:/users/paul.conn/git/BOSSst/output/BOSS2012_noST_mcmc_workspace.Rdata')
load('c:/users/paul.conn/git/BOSSst/output/BOSS2012_noST_MCMC_output.Rdata')
load('c:/users/paul.conn/git/BOSSst/Mapping2012_ordered.Rda')  #need ordered Mapping from hierchical_boss_st because Dat from the workspace is ordered 

n.mcmc.iter=nrow(MCMC$MCMC)
n.pred.iter=dim(Post$N)[2]
n.obs=nrow(Dat)

#assemble count dataset
n.sampled = nrow(Mapping)
Counts = data.frame(Ribbon = rep(0,n.sampled),
                    Ringed = rep(0,n.sampled),
                    Spotted = rep(0,n.sampled),
                    Bearded = rep(0,n.sampled),
                    Cell = Mapping[,1],
                    Day = Mapping[,2])
for(irow in 1:n.obs){
  if(Dat[irow,"Species"]==1)Counts[Dat[irow,"Transect"],"Spotted"]=Counts[Dat[irow,"Transect"],"Spotted"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==2)Counts[Dat[irow,"Transect"],"Ribbon"]=Counts[Dat[irow,"Transect"],"Ribbon"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==3)Counts[Dat[irow,"Transect"],"Bearded"]=Counts[Dat[irow,"Transect"],"Bearded"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==4)Counts[Dat[irow,"Transect"],"Ringed"]=Counts[Dat[irow,"Transect"],"Ringed"]+Dat[irow,"Group"]
}



###########
# plot ice covariate, surveys, estimates as a function of 4 different time steps
#########
S = nrow(Data$Grid$y2012[[1]])
require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),legend.title=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
library(RColorBrewer)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))
library(gridExtra)


#create time lapsed video of covariate, transect counts, abundance
Dates=data.frame(matrix(0,29,3))
colnames(Dates)=c("x","y","date")
Dates[,"date"]=c(paste(c(10:30),"April 2012"),paste(c(1:8),"May 2012"))
#spotted animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Spotted"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[1,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,2000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Spotted_Bering_2012.mp4"
)

#ribbon animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ribbon"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[2,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ribbon seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,2000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Ribbon_Bering_2012.mp4"
)

#bearded animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Bearded"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[3,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Bearded seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,600),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Bearded_Bering_2012.mp4"
)

#ringed animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2012[[it+5]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/22 but analysis goes 4/10-5/9
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ringed"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[4,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ringed seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1000),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Ringed_Bering_2012.mp4"
)


#####   2013 Animations  ############
load('BOSS_data_2013.Rda')
load("AlaskaBeringData2012_2013_14Dec2015.Rdat")  #boss grid, ice data

#load MCMC results from one chain
load('c:/users/paul.conn/git/BOSSst/output/BOSS2013_mcmc_workspace.Rdata')
load('c:/users/paul.conn/git/BOSSst/Mapping2013_ordered.Rda')  #need ordered Mapping from hierchical_boss_st because Dat from the workspace is ordered 

n.mcmc.iter=nrow(MCMC$MCMC)
n.pred.iter=dim(Post$N)[2]
n.obs=nrow(Dat)

#assemble count dataset
n.sampled = nrow(Mapping)
Counts = data.frame(Ribbon = rep(0,n.sampled),
                    Ringed = rep(0,n.sampled),
                    Spotted = rep(0,n.sampled),
                    Bearded = rep(0,n.sampled),
                    Cell = Mapping[,1],
                    Day = Mapping[,2])
for(irow in 1:n.obs){
  if(Dat[irow,"Species"]==1)Counts[Dat[irow,"Transect"],"Spotted"]=Counts[Dat[irow,"Transect"],"Spotted"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==2)Counts[Dat[irow,"Transect"],"Ribbon"]=Counts[Dat[irow,"Transect"],"Ribbon"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==3)Counts[Dat[irow,"Transect"],"Bearded"]=Counts[Dat[irow,"Transect"],"Bearded"]+Dat[irow,"Group"]
  if(Dat[irow,"Species"]==4)Counts[Dat[irow,"Transect"],"Ringed"]=Counts[Dat[irow,"Transect"],"Ringed"]+Dat[irow,"Group"]
}



###########
# plot ice covariate, surveys, estimates as a function of 4 different time steps
#########
S = nrow(Data$Grid$y2013[[1]])
require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid$y2013[[1]],byid=TRUE))

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),legend.title=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
library(RColorBrewer)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))
library(gridExtra)


#create time lapsed video of covariate, transect counts, abundance
#2013 survey counts go from 4/7-5/6 but environmental covariate data goes 4/5-5/9
Dates=data.frame(matrix(0,30,3))
colnames(Dates)=c("x","y","date")
Dates[,"date"]=c(paste(c(7:30),"April 2013"),paste(c(1:6),"May 2013"))
t.steps=30
#spotted animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  #covariate grid goes 4/5-5/9 but analysis goes 4/7-5/6
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Spotted"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[1,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1600),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Spotted_Bering_2013.mp4"
)

#ribbon animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ribbon"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[2,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ribbon seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1500),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Ribbon_Bering_2013.mp4"
)

#bearded animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Bearded"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[3,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Bearded seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,800),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Bearded_Bering_2013.mp4"
)

#ringed animation
saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid$y2013[[it+2]]@data[,"ice_conc"]  
    # Count - unobserved cells remain NAs
    Cur.count.data=Counts[which(Counts[,"Day"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Day"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Ringed"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=apply(Post$N[4,,((it-1)*S+1):(it*S)],2,'median')
    
    Centroids=data.frame(gCentroid(Data$Grid$y2012[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=bluePalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Ringed seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=orangePalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,500),colours=orangePalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Ringed_Bering_2013.mp4"
)


