## produce flight table (date each flight started)

library(sp)
library(rgeos)
library(nPacMaps)

#load('./data_from_JML/boss_hotspots_sp.Rda') #hotspots
load('./data_from_JML/boss_geo_sp.Rda') #on effort points

Flts = sort(unique(boss_geo@data[,"flightid"]))
Times = rep(boss_geo$img_dt[1],length(Flts))

for(i in 1:length(Flts)){
  crap=boss_geo[which(boss_geo$flightid==Flts[i]),]
  Times[i]=min(crap$img_dt)
}

Times=as.Date(Times,tz="PST8PDT")

Flt_table = data.frame(Flt = Flts, Dates = Times)
save(Flt_table,file="BOSS_Flt_dates.Rda")