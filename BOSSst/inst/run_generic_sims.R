# run_generic_sims.R
# script to run generic spatio-temporal count data simulations
require(STabundance)

set.seed(12345)
n.sims=1 #number of simulations at each design point
#Delta=c(-0.02,-0.02+0.04*c(1:(n.sims-1))/(n.sims-1)) #equally spaced from 2% decrease to 2% increase in habitat
Delta=c(-0.01)
S=400
t.steps=20
N.transects=c(2)
line.width=0.05
n.species=4
GENERATE=TRUE
if(GENERATE==TRUE){
  for(igen in 1:1){  #loop over generating model to generate data sets
    for(itrans in 1:1){ #loop over number of transects in each cell
      for(isim in 1:n.sims){
        Sim.data=sim_data_generic(n.species=4,S=S,t.steps=t.steps,n.transects=N.transects[itrans],line.width=line.width,buffer=0,delta=Delta[isim])
        #fname=paste("simdata_gen",Model.list[igen],"_trans",N.transects[itrans],"_sim",isim,sep='')
        #save(Sim.data,file=paste("./sim_generic_data/",fname,sep=''))
      }
    }  
  }
}

n.transects=nrow(Sim.data$Effort)


#for(igen in 1:3){  #loop over generating model to generate data sets
#  for(itrans in 1:1){ #loop over number of transects in each cell
#    for(isim in 1:(n.sims/2)){  

#      if(igen==3){
#        fname=paste("simdata_gen",Model.list[igen-1],"_trans",N.transects[itrans],"_sim",isim,sep='')
#        load(paste("./sim_generic_data/",fname,sep=''))
#        Which.distances=Sim.data$Data$Which.distances  
#      }

#      fname=paste("simdata_gen",Model.list[igen],"_trans",N.transects[itrans],"_sim",isim,sep='')
#      load(paste("./sim_generic_data/",fname,sep=''))
Data=Sim.data$Data
#      if(igen==3)Data$Which.distances=Which.distances

n.knots=length(Data$Knot.locations)
#calculate kernel densities at grid cell centroids 
Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
Distances=gDistance(Data$Knot.locations,Cell.centroids,byid=TRUE)
K=matrix(dnorm(Distances,0,5),S,n.knots)  #knot sd=5 
K=K/rowSums(K)          
Data$K=K

#Data$Count.data=Data$Count.data[-which(Data$Count.data[,"Time"]%in%c(3,4,10)),]
DayHour=data.frame(day=rep(1,n.transects),hour=rep(1,n.transects))
Thin=array(1,dim=c(n.species,2,2,1000))	#need >1 day, hour for matrix subscripting to work right in mcmc function
Psi=0.8*diag(4)
Psi[1,2]=0.1
Psi[1,3]=Psi[1,4]=0.05
Psi[2,1]=0.1
Psi[2,3]=Psi[2,4]=0.05
Psi[3,1]=Psi[3,2]=0.05
Psi[3,4]=0.1
Psi[4,1]=Psi[4,2]=0.05
Psi[4,3]=0.1
Psi=array(Psi,dim=c(4,4,1))

Hab.cov=data.frame(matrix(0,S*t.steps,ncol(Data$Grid[[1]]@data)))
colnames(Hab.cov)=colnames(Data$Grid[[1]]@data)
for(it in 1:t.steps){
  Hab.cov[((it-1)*S+1):((it-1)*S+S),]=Data$Grid[[it]]@data
}
hab.formula=vector("list",n.species)
for(isp in 1:n.species)hab.formula[[isp]]=~0+matern+matern2
Cov.prior.parms=array(0,dim=c(n.species,2,1))
Cov.prior.parms[,1,1]=2  #we'll put priors a little off from their true values; expected group sizes are 4 and 2 for each species
Cov.prior.parms[,2,1]=1
#Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=matrix(0,n.species,1)  
Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
Cov.prior.n=matrix(2,n.species,1)
spat.ind=FALSE
srr.tol=0.5

Prior.pars=list(beta.tau=0.01,
                a.eps=1,
                b.eps=0.01,
                a.eta=1,
                b.eta=0.01,
                beta0.tau.rw2=1,
                beta1.tau.rw2=10)

Control=list(iter=5000,burnin=10,thin=10,predict=TRUE,MH.N=rep(0.2,n.species),MH.omega=matrix(0.05,n.species,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE,species.optim=TRUE)        

Dat=Sim.data$Obs
Area.hab=rep(1,S*t.steps)
Mapping=Sim.data$Effort[,c("Cell","Time")]
Area.trans=Sim.data$Effort$AreaSurveyed
Prop.photo=rep(0.5,n.transects)
Obs.cov=NULL
Hab.formula=hab.formula
n.obs.cov=0
Psi=Psi
Inits=NULL
grps=TRUE
post.loss=TRUE
Area.trans=Sim.data$Effort[,"AreaSurveyed"]

set.seed(12345)
MCMC=hierarchical_boss_st(Dat=Sim.data$Obs,K=Data$K,Area.hab=rep(1,S*t.steps),Area.trans=Sim.data$Effort[,"AreaSurveyed"],Mapping=Sim.data$Effort[,c("Cell","Time")],DayHour=DayHour,Thin=Thin,Prop.photo=rep(0.5,n.transects),Hab.cov=Hab.cov,Obs.cov=NULL,Hab.formula=hab.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,n.species=n.species,n.obs.cov=0,spat.ind=spat.ind,Psi=Psi,Inits=NULL,grps=TRUE,Control=Control,Prior.pars=Prior.pars,post.loss=TRUE)
    

#Obs.data
#Obs.data=Data$Count.data[which(Data$Count.data[,"Count"]>0),]
#Cov=rep(0,nrow(Obs.data))
#for(i in 1:nrow(Obs.data))Cov[i]=Data$Grid[[Obs.data[i,"Time"]]]@data[Obs.data[i,"Cell"],"matern"]

#calculate posterior for N
MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
N.true=apply(Sim.data$N,2,'sum')
N.est=apply(MCMC$MCMC$N,1,'mean')
plot(N.est)
lines(N.true)

#plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

#plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
#plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(15,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(30,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")


