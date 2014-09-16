# run_BOSS_ST_sims.R
# script to run generic spatio-temporal count data simulations
require(STabundance)
source('./BOSSst/R/mcmc_boss_st.R')
source('./BOSSst/R/hierarchical_boss_st.R')  
source('./BOSSst/R/util_funcs.R')
source('./BOSSst/R/sim_funcs.R')


data(AlaskaBeringData2012_17April2014)  #boss grid, ice data
data(Knot_cell_distances) #load object giving K matrix, Q for knots
load("Effort2012_BOSSst_5Sep2014.Rdata")  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)
load('c:/users/paul.conn/git/BOSS/BOSS/data/p13.Rdata')  #read in confusion array
load('c:/users/paul.conn/git/BOSSst/Haulout_samples.Rdat')  #read in haulout proportion MCMC samples

set.seed(123454)
t.steps=29
n.species=5  #last is `other'
n.transects=length(Effort$Area.trans)
n.zeros=1000  #number of 'extra zeros' to put in (max is around 3000; to use max set n.zeros=NA)
Old.Grid=Data$Grid
Data$Grid=vector('list',t.steps)
for(it in 1:t.steps){
  Data$Grid[[it]]=Old.Grid[[it+3]]
  Data$Grid[[it]]@data$ice2=Data$Grid[[it]]@data$ice_conc^2
  Data$Grid[[it]]@data$sqrt_edge=Data$Grid[[it]]@data$dist_edge^(0.5)
}
#Data$Count.data=Effort$Count.data
rm(Old.Grid)

S=nrow(Data$Grid[[1]])

Effort$Area.hab[which(Effort$Area.hab==0)]=0.00001
n.knots=length(Data$Knot.locations)
Data$K=K.data$K
Area.trans=Effort$Area.trans
Area.hab=Effort$Area.hab
DayHour=data.frame(day=Effort$Mapping[,2],hour=Effort$Mapping[,3])
Mapping=Effort$Mapping[,1:2]
Dat=Effort$Count.data
rm(Effort)
#change observations from NA to 'unknown' for observations with photo =1 and obs=NA to other
Temp=which(Dat[,"Photo"] & is.na(Dat$Obs))
if(length(Temp>0))Dat[Temp,"Obs"]=14

#add in some 'zero' data to anchor model in places where there's no ice
Temp=which(Data$Grid[[1]]@data[,"ice_conc"]<.001)
Which.no=matrix(1,length(Temp),2)
Which.no[,1]=Temp
for(it in 2:length(Data$Grid)){
  Temp=which(Data$Grid[[it]]@data[,"ice_conc"]<.001)
  Cur.mat=matrix(it,length(Temp),2)
  Cur.mat[,1]=Temp
  Which.no=rbind(Which.no,Cur.mat)
}
if(is.na(n.zeros)==FALSE)Which.no=Which.no[sample(c(1:nrow(Which.no)),n.zeros),]


#get rid of these if they were actually sampled [already in dataset]
Mapping.1d=(Mapping[,2]-1)*S+Mapping[,1]  
Which.no.1d=(Which.no[,2]-1)*S+Which.no[,1]
I.sampled=Which.no.1d%in%Mapping.1d
I.sampled2=Mapping.1d%in%Which.no.1d
Area.trans[which(I.sampled2==TRUE)]=0.9  #change area sampled in 0 ice cells to 0.9
Which.no.1d=Which.no.1d[which(I.sampled==0)]
Which.no=Which.no[which(I.sampled==0),]
#now add in 0 ice cells to list of sampled cells and adjust related quantities
Mapping=rbind(Mapping,Which.no)
Area.trans=c(Area.trans,rep(.9,length(Which.no.1d)))
n.transects=n.transects+length(Which.no.1d)
Tmp=matrix(1,nrow(Mapping),2)
Tmp[1:nrow(DayHour),]=as.matrix(DayHour)
Tmp[,1]=Mapping[,2]
DayHour=Tmp


#Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
#Distances=gDistance(Data$Knot.locations,Cell.centroids,byid=TRUE)
#K=matrix(dnorm(Distances,0,5),S,n.knots)  #knot sd=5 
#K=K/rowSums(K)          
#Data$K=K

#generate thinning priors using haulout and det prob data
Thin=array(1,dim=c(n.species,29,24,1000))
P=rbeta(1000,67,5)  #conjugate beta(1,1) for binomial detection data (66/70 successes)
for(isp in 1:4){
  for(iday in 1:29){
    for(ihr in 1:24){
      Thin[isp,iday,ihr,]=P
      if(isp==1)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$spotted[iday,ihr,]
      if(isp==2)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$ribbon[iday,ihr,]
      if(isp==3)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$bearded[iday,ihr,]
    }
  }
}
rm(Haulout.samples) 



#Define misclassification matrix  [species observations follow p13 order - spotted (guess, likely pos), ribbon, bearded, ringed]
Psi=array(0,dim=c(5,14,dim(p13)[3]))
Psi[1:4,1:13,]=p13[1:4,1:13,]
Psi[5,14,]=1  #'other' always gets a 14 
rm(p13)




Hab.cov=data.frame(matrix(0,S*t.steps,ncol(Data$Grid[[1]]@data)))
colnames(Hab.cov)=colnames(Data$Grid[[1]]@data)
for(it in 1:t.steps){
  Hab.cov[((it-1)*S+1):((it-1)*S+S),]=Data$Grid[[it]]@data
}
#TEMPORARILY INCREASE ICE in 2 cells to decrease sea ice effect at 0 sea ice concentration
Hab.cov[c(18822,21975),"ice_conc"]=0.5
Hab.cov[c(18822,21975),"ice2"]=0.25

hab.formula=vector("list",n.species)

for(isp in 1:n.species)hab.formula[[isp]]=~0+ice_conc+ice2+sqrt_edge+dist_shelf+dist_mainland+dist_contour
Cov.prior.parms=array(0,dim=c(n.species,2,1))
Cov.prior.parms[,1,1]=1  #gamma(1,1) on E(group size)-1
Cov.prior.parms[,2,1]=1
#Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=matrix(0,n.species,1)  
Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
Cov.prior.n=matrix(2,n.species,1)
spat.ind=FALSE


Prior.pars=list(beta.tau=0.0001,
                a.eps=1,
                b.eps=0.01,
                a.eta=1,
                b.eta=0.01,
                beta0.tau.rw2=1,
                beta1.tau.rw2=10)

Control=list(iter=125000,burnin=2000,thin=50,n.adapt=2000,predict=TRUE,MH.N=rep(0.2,n.species),adapt=TRUE,fix.tau.epsilon=FALSE,species.optim=TRUE)        


Dat[,4]=as.numeric(as.character(Dat[,4]))
Prop.photo=rep(0.5,n.transects)
Obs.cov=NULL
Inits=NULL
grps=TRUE
K=Data$K
Hab.formula=hab.formula
n.obs.cov=0
Psi=Psi
post.loss=FALSE
#True.sp=Sim.data$True.sp #if not null, sets all observations to have true species values (for debugging)
True.sp=NULL
True.species=True.sp
#Omega.true=Sim.data$Omega.true
Omega.true=NULL
#Eta.true=Sim.data$Eta.true
Eta.true=NULL
#Alpha.true=Sim.data$Alpha.true
Alpha.true=NULL
Surveyed=c(1:1197)
DayHour[which(DayHour[,2]==0),2]=24

set.seed(12345)
MCMC=hierarchical_boss_st(Dat=Dat,K=Data$K,Area.hab=Area.hab,Area.trans=Area.trans,Mapping=Mapping,DayHour=DayHour,Thin=Thin,Prop.photo=rep(0.5,n.transects),Hab.cov=Hab.cov,Obs.cov=NULL,Hab.formula=hab.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,n.species=n.species,n.obs.cov=0,spat.ind=spat.ind,Psi=Psi,Inits=NULL,grps=TRUE,Control=Control,Prior.pars=Prior.pars,post.loss=post.loss,Surveyed=Surveyed,True.species=True.sp,Omega.true=Omega.true,Eta.true=Eta.true,Alpha.true=Alpha.true,DEBUG=TRUE)


out.file="./Bering2012_1camera_MCMC.Rdata"
save(MCMC,file=out.file)



Cur.G=matrix(apply(MCMC$Post$G[isp,,],2,'median'),S,t.steps)
for(i in 1:t.steps){
  plot_N_map(i,Cur.G,Grid=Data$Grid)
}
#calculate posterior for N
#MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
#N.est=apply(MCMC$MCMC$N,1,'mean')
#plot(N.est)

#crap=acf(MCMC$MCMC$N[18,200:799],lag.max=6)
#(1+2*sum(crap$acf)*var(MCMC$MCMC$N[18,200:799])*10000)/(mean(MCMC$MCMC$N[18,200:799]))^2

#plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

#plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
#plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(17,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(29,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")

