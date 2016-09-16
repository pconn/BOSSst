require(STabundance)

source('./BOSSst/R/util_funcs.R')
data(AlaskaBeringData2012_17April2014)  #boss grid, ice data

Control=list(fname='./output/Bering2012_1camera_PostPred_noST_1week_fixTau_Outliers_oldCovs',n.files=10) 

PostN=cat_preds(fname=Control$fname,n.files=Control$n.files)
load("./output/Bering2012_1camera_PostPred_noST_1week_fixTau_Outliers_oldCovs.Rdata")

#rows of PostN: spotted, ribbon, bearded, ringed
isp=4
S=1299
t.steps=7
#N=apply(PostN,c(1,2),'sum')
N=rowSums(PostN[isp,,])/t.steps
summary(N)
Cur.G=matrix(apply(PostN[isp,1:dim(PostN)[2],],2,'median'),S,t.steps)
#Cur.G[,1]=Data$Grid[[1]]@data[,"ice_conc"]
#Cnt=rep(0,S)
#for(i in 1:nrow(Dat))Cnt[Mapping[Dat[i,"Transect"],1]]=Cnt[Mapping[Dat[i,"Transect"],1]]+1
#Cur.G[,1]=Cnt
#for(i in 1:t.steps){
plot_N_map(1,Cur.G,Grid=Data$Grid)
plot_N_map(5,Cur.G,Grid=Data$Grid)
plot_N_map(10,Cur.G,Grid=Data$Grid)
plot_N_map(15,Cur.G,Grid=Data$Grid)
plot_N_map(22,Cur.G,Grid=Data$Grid)
plot_N_map(27,Cur.G,Grid=Data$Grid)
#}
#calculate posterior for N

#look at covariates day 1 = 4, day 29 = 32
iday=32
plot_N_map(1,matrix(Data$Grid[[iday]]@data[,"dist_shelf"],S,1),Grid=Data$Grid)


#makeshift estimate of N just using estimates from t=1
N=rowSums(PostN[isp,,1:S])
plot(N)

summary(N)

#MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
#N.est=apply(MCMC$MCMC$N,1,'mean')
#plot(N.est)

#  look at relative time effect on counts (controling for density) using mgcv - requires Dat, Mapping from e.g. run_Bering2012
if(isp==1)Which.spotted=which(Dat[,"Obs"]<4)
if(isp==2)Which.spotted=which(Dat[,"Obs"]>3 & Dat[,"Obs"]<7)
if(isp==3)Which.spotted=which(Dat[,"Obs"]>6 & Dat[,"Obs"]<10)
if(isp==4)Which.spotted=which(Dat[,"Obs"]>9 & Dat[,"Obs"]<13)
Dat.gam=Dat[Which.spotted,]
N=rep(0,nrow(Mapping))
for(i in 1:nrow(Mapping))N[i]=Cur.G[Mapping[i,1],Mapping[i,2]]
Count.df=data.frame(Count=0,Cell=Mapping[,1],Day=Mapping[,2],Offset=Area.trans,N=N)
for(i in 1:nrow(Dat.gam)){
  Count.df[Dat.gam[i,"Transect"],"Count"]=Count.df[Dat.gam[i,"Transect"],"Count"]+1
}
library(mgcv)
spotted.gam <- gam(Count~offset(log(Offset))+N+s(Day,k=5),data=Count.df,family="tw")

plot(spotted.gam)

par(mfrow=c(2,1))
hist(Count.df[,"Day"],breaks=28,main="Cells surveyed",xlab="Day")
Counts.by.day=rep(0,29)
for(i in 1:29)Counts.by.day[i]=sum(Count.df[which(Count.df[,"Day"]==i),"Count"])
plot(Counts.by.day,xlab="Day",ylab="Count")



#crap=acf(MCMC$MCMC$N[18,200:799],lag.max=6)
#(1+2*sum(crap$acf)*var(MCMC$MCMC$N[18,200:799])*10000)/(mean(MCMC$MCMC$N[18,200:799]))^2

#plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

