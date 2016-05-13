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
load('p13.RData')  #read in confusion array
load('Haulout_samples.Rdat')  #read in haulout proportion MCMC samples

set.seed(123454)
t.steps=29
n.species=5  #last is `other'
n.transects=length(Effort$Area.trans)
n.zeros=100  #number of 'extra zeros' to put in (max is around 3000; to use max set n.zeros=NA)
Old.Grid=Data$Grid
S=nrow(Data$Grid[[1]])
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

#look for density outliers
C.tot=rep(0,n.transects)
for(i in 1:n.transects)C.tot[i]=length(which(Dat[,"Transect"]==i))
C.tot=C.tot/Area.trans
summary(C.tot)

#delete effort for cells for which Area.trans<0.001
Which.yes=which(Area.trans>0.001)
New.pl=rep(0,n.transects)
Area.trans=Area.trans[Which.yes]
Mapping=Mapping[Which.yes,]
New.pl[Which.yes]=c(1:length(Which.yes))
Dat[,"Transect"]=New.pl[Dat[,"Transect"]]
Dat=Dat[-which(Dat[,"Transect"]==0),]  #remove observations where transect area <0.001
n.transects=length(Area.trans)

#restrict to seals w photos only
Dat=Dat[which(Dat[,"Photo"]==1),]
Dat=Dat[-which(Dat[,"Obs"]==14),]

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
#Area.trans[which(I.sampled2==TRUE)]=0.9  #change area sampled in 0 ice cells to 0.9 #no, this messes up initial density estimates
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
Which.sample=sample(c(1:dim(p13)[3]),1000)  #reduce posterior to 1000 samples to reduce memory requirements
Psi=array(0,dim=c(5,14,1000))
Psi[1:4,1:13,]=p13[1:4,1:13,Which.sample]  
Psi[5,14,]=1  #'other' always gets a 14 
rm(p13)




Hab.cov=data.frame(matrix(0,S*t.steps,ncol(Data$Grid[[1]]@data)))
colnames(Hab.cov)=colnames(Data$Grid[[1]]@data)
for(it in 1:t.steps){
  Hab.cov[((it-1)*S+1):((it-1)*S+S),]=Data$Grid[[it]]@data
}

#look for seals seen where not much ice - cells that have ice<0.01 where seals are seen are replaced with average ice values for neighboring cells
Which.less=which(Hab.cov[Mapping.1d,"ice_conc"]<0.01)
for(iobs in 1:length(Which.less)){
  cur.s=Mapping[Which.less[iobs],1]
  cur.t=Mapping[Which.less[iobs],2]
  new.ice=Data$Adj[cur.s,]%*%Data$Grid[[cur.t]][["ice_conc"]]/sum(Data$Adj[cur.s,])
  Hab.cov[S*(cur.t-1)+cur.s,"ice_conc"]=new.ice
  Hab.cov[S*(cur.t-1)+cur.s,"ice2"]=new.ice^2
  cat(paste0('replaced cell ',cur.s,' time ',cur.t,' with original ice_conc=',Data$Grid[[cur.t]]@data[cur.s,"ice_conc"]," with ice=",new.ice,'\n'))
}
#take a look at where these occur
#plot_N_map(1,as.matrix(Data$Grid[[17]][["ice_conc"]],ncol=1),highlight=1191,Grid=Data$Grid)
#plot_N_map(1,as.matrix(Data$Grid[[15]][["ice_conc"]],ncol=1),highlight=636,Grid=Data$Grid)


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
                beta0.tau.rw2=2,
                beta1.tau.rw2=50)

Control=list(iter=550000,burnin=50000,thin=250,n.adapt=2500,predict=TRUE,MH.N=rep(0.2,n.species),adapt=TRUE,fix.tau.epsilon=FALSE,species.optim=TRUE,update.sp=TRUE,est.alpha=FALSE,n.files=10,fname="./output/Bering2012_1camera_PostPred")        


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

#set initial true species to help with algorithm convergence
True.sp=NULL
True.sp=rep(0,nrow(Dat))
for(i in 1:nrow(Dat)){
  if(is.na(Dat[i,"Obs"])==FALSE){  
    if(Dat[i,"Obs"]%in%c(1,2,3))True.sp[i]=1
    if(Dat[i,"Obs"]%in%c(4,5,6))True.sp[i]=2
    if(Dat[i,"Obs"]%in%c(7,8,9))True.sp[i]=3
    if(Dat[i,"Obs"]%in%c(10,11,12))True.sp[i]=4
    if(Dat[i,"Obs"]==13)True.sp[i]=sample(c(1:4),1)
    if(Dat[i,"Obs"]==14)True.sp[i]=5   
  }
}
#now tabulate counts by species and ecoregion to help set species for non-photographed
n.eco.bins=length(tabulate(Hab.cov[,"Ecoregion"]))
Tab.species=matrix(0,5,n.eco.bins)
for(isp in 1:5){
  Cur.dat=Dat[which(is.na(Dat[,"Obs"])==FALSE & True.sp==isp),]
  Cur.map=Mapping[Cur.dat$Transect,]
  Cur.eco=Hab.cov[(Cur.map[,2]-1)*S+Cur.map[,1],"Ecoregion"]
  Tab.species[isp,]=tabulate(Cur.eco,nbins=n.eco.bins)
}
for(i in 1:nrow(Dat)){
  if(is.na(Dat[i,"Obs"])){  
    cur.map=Mapping[Dat[i,"Transect"],]
    cur.eco=Hab.cov[(cur.map[2]-1)*S+cur.map[1],"Ecoregion"]
    True.sp[i]=which(rmultinom(1,1,prob=Tab.species[,cur.eco])==1)
  }
}

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


out.file="./output/Bering2012_1camera_MCMC.Rdata"
save(MCMC,file=out.file)


#PostN=cat_preds(fname=Control$fname,n.files=Control$n.files)
  
#Cur.G=matrix(apply(PostN[4,1000:2000,],2,'median'),S,t.steps)
#for(i in 1:t.steps){
#plot_N_map(1,Cur.G,Grid=Data$Grid)
#plot_N_map(5,Cur.G,Grid=Data$Grid)
#plot_N_map(10,Cur.G,Grid=Data$Grid)
#plot_N_map(15,Cur.G,Grid=Data$Grid)
#plot_N_map(22,Cur.G,Grid=Data$Grid)
#plot_N_map(27,Cur.G,Grid=Data$Grid)
#}
#calculate posterior for N
#MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
#N.est=apply(MCMC$MCMC$N,1,'mean')
#plot(N.est)

#crap=acf(MCMC$MCMC$N[18,200:799],lag.max=6)
#(1+2*sum(crap$acf)*var(MCMC$MCMC$N[18,200:799])*10000)/(mean(MCMC$MCMC$N[18,200:799]))^2

#plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

#Plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
#plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(1,apply(MCMC$Post$$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(17,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(29,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")

