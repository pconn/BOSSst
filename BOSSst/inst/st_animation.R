#make ice + tracks animation

library(animation)
library(ggplot2)
library(plyr)
load("BOSSst_2012data.Rdata")  #boss grid, ice data

load('c:/users/paul.conn/git/BOSS/BOSS/data/Transect_Data.Rdat') #read in effort data (held in "flight_segs" SpLinesDF)

#assemble a list of unique flight names and first/last on effort dates for each
Flt.ids=unique(flight_segs[["flightid"]])
Flt.table=data.frame(id=Flt.ids,min.time=rep(flight_segs[["min_dt"]][1],length(Flt.ids)),max.time=rep(flight_segs[["max_dt"]][1],length(Flt.ids)))
for(i in 1:nrow(Flt.table)){
  Cur.which=which(flight_segs[["flightid"]]==Flt.table[i,"id"])
  Flt.table[i,"min.time"]=min(flight_segs[["min_dt"]][Cur.which])
  Flt.table[i,"max.time"]=max(flight_segs[["max_dt"]][Cur.which])
}
Day=as.numeric(ceiling((Flt.table[,"min.time"]-min(Flt.table[,"min.time"])+.01)/(3600*24)))
Flt.table=cbind(Flt.table,Day)
Airplane=substr(as.character(flight_segs@data[,"flightid"]),1,2)
flight_segs@data=cbind(flight_segs@data,Airplane)

#fix some issues with OtterFl06 having effort off grid
#which.rows=which(flight_segs@data[,"flightid"]=="OtterFl06")
#plot(flight_segs[which.rows[c(1:2,3,4:7)],])
flight_segs=flight_segs[-which.rows[3],]

#Ice - spdf of ice data, Flts = spatial lines data frame giving all flights on a given date
make_plot<-function(Ice,Flts=NULL){
  Ice@data$id=rownames(Ice@data)
  tmp1<-fortify(Ice,region="id")
  tmp2<-join(tmp1,Ice@data,by="id")
  colnames(tmp2)[1:2]=c("Easting","Northing")
  tmp.plot=ggplot(data=tmp2,aes(x=Easting,y=Northing))+geom_raster(aes(fill=ice_conc))
  if(is.null(Flts)==FALSE){
    n.fl=nrow(Flts@data)
    coords=coordinates(Flts[1,])[[1]][[1]]
    tmp1=data.frame(matrix(0,nrow(coords)-1,5))
    colnames(tmp1)=c("x","y","xend","yend","Airplane")
    for(iseg in 1:(nrow(coords)-1)){
      tmp1[iseg,]=c(coords[iseg,],coords[iseg+1,],Flts@data[1,"Airplane"])
    }
    if(n.fl>1){
      for(ifl in 2:n.fl){
        coords=coordinates(Flts[ifl,])[[1]][[1]]
        tmp3=data.frame(matrix(0,nrow(coords)-1,5))
        colnames(tmp3)=c("x","y","xend","yend","Airplane")
        for(iseg in 1:(nrow(coords)-1)){
          tmp3[iseg,]=c(coords[iseg,],coords[iseg+1,],Flts@data[ifl,"Airplane"])
        }
        tmp1=rbind(tmp1,tmp3)
      }
    }
    tmp1[,"Airplane"]=factor(as.character(tmp1[,"Airplane"]),levels=c("1","2"))
    tmp.plot=tmp.plot+geom_segment(data=tmp1,mapping=aes(x=x, y=y, xend=xend, yend=yend,colour=Airplane),size=0.5,show_guide=FALSE)+scale_colour_manual(labels=c("Twin otter","Aerocommander"),values = c("1"="red","2"="orange"))
  }
  tmp.plot
}

Flight.days=sort(unique(Flt.table[,"Day"]))
setwd('c:/users/paul.conn/git/BOSSst/ice_flt_plots')
for(iday in 1:36){ #limit to 4/7-5/12 (missing ice after 5/13)
  if(iday%in%Flight.days){
    cur.id=as.character(Flt.table[which(Flt.table[,"Day"]==iday),"id"])
    cur.segs=which(as.character(flight_segs@data[,"flightid"])%in%cur.id)
    cur.plot=make_plot(Ice=Data$Grid[[3+iday]],Flts=flight_segs[cur.segs,])  
  }
  else cur.plot=make_plot(Ice=Data$Grid[[3+iday]]) 
  pdf(file=paste("ice_flts_2013_",iday,".pdf",sep=''))
  print(cur.plot,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  dev.off()
} 




