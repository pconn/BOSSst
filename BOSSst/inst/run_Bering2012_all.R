#format data for BOSS spatio-temporal analysis
load('./BOSS_data_2012.Rda')
library(hierarchicalDS) #a couple of utility functions are used from here
source('./BOSSst/R/util_funcs.R')
source('./BOSSst/R/hierarchical_boss_st.R')
source('./BOSSst/R/mcmc_boss_st.R')

t_steps=29
n_species=4  
n_transects=length(Area_trans)
S=length(Area_hab)
set.seed(12345)


#take a look at where these occur
#plot_N_map(1,as.matrix(Data$Grid[[17]][["ice_conc"]],ncol=1),highlight=1191,Grid=Data$Grid)
#plot_N_map(1,as.matrix(Data$Grid[[15]][["ice_conc"]],ncol=1),highlight=636,Grid=Data$Grid)


hab_formula=vector("list",n_species)

for(isp in 1:n_species)hab_formula[[isp]]=~0+ice_conc+ice2+sqrt_edge+dist_shelf+dist_mainland+dist_contour
Cov_prior_parms=array(0,dim=c(n_species,2,1))
Cov_prior_parms[,1,1]=1  #gamma(1,1) on E(group size)-1
Cov_prior_parms[,2,1]=1
#Cov_prior_parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
Cov_prior_fixed=matrix(0,n_species,dim(Cov_prior_parms)[3])
Cov_prior_pdf=matrix(0,n_species,1)  
Cov_prior_pdf[,1]="pois1"  #model group size as a zero truncated poisson
Cov_prior_n=matrix(2,n_species,1)
spat_ind=FALSE


Prior_pars=list(beta.tau=0.0001,
                a.eps=1,
                b.eps=0.01,
                a.eta=100,
                b.eta=0.01,
                beta0.tau.rw2=2,
                beta1.tau.rw2=5000)

Control=list(iter=1200000,burnin=700000,thin=5,n.adapt=600000,predict=TRUE,MH.N=rep(0.2,n_species),adapt=TRUE,fix.tau.epsilon=rep(TRUE,4),species.optim=TRUE,update.sp=TRUE,est.alpha=F,n.files=10,fname="./output/Bering2012all_test")        
Control$fix.tau.linear=TRUE
Control$Tau.eps.end=c(20,20,20,20)
Control$save.file="BOSS2012_mcmc_workspace.RData"
Control$save.every=10000
Control$post.loss=TRUE
Control$GOF=TRUE
Control$gIVH=0
#Control$start.file="BOSS2012_mcmc_workspace2.RData"
Control$start.file=NULL

Dat[,4]=as.numeric(Dat[,4])
Prop_photo=rep(0.5,n_transects)
Obs_cov=NULL
grps=TRUE
Hab_formula=hab_formula
n_obs_cov=0
#Psi=Psi

#True.sp=Sim.data$True.sp #if not null, sets all observations to have true species values (for debugging)

#set initial true species to help with algorithm convergence
#True_sp=NULL
True_sp=rep(0,nrow(Dat))
for(i in 1:nrow(Dat)){
  if(is.na(Dat[i,"Obs"])==FALSE){  
    if(Dat[i,"Obs"]%in%c(1,2,3))True_sp[i]=1
    if(Dat[i,"Obs"]%in%c(4,5,6))True_sp[i]=2
    if(Dat[i,"Obs"]%in%c(7,8,9))True_sp[i]=3
    if(Dat[i,"Obs"]%in%c(10,11,12))True_sp[i]=4
    if(Dat[i,"Obs"]==13)True_sp[i]=sample(c(1:4),1)
    #if(Dat[i,"Obs"]==14)True.sp[i]=5   
  }
}


True_species=True_sp
#Omega.true=Sim.data$Omega.true
Omega_true=NULL
#Eta.true=Sim.data$Eta.true
Eta_true=NULL
#Alpha.true=Sim.data$Alpha.true
Alpha_true=NULL
#DayHour[which(DayHour[,2]==0),2]=24



Area.hab=Area_hab
Area.trans=Area_trans
Prop.photo=rep(1.0,n_transects)
Hab.cov=Hab_cov
Obs.cov=NULL
Hab.formula=hab_formula
Cov.prior.pdf=Cov_prior_pdf
Cov.prior.parms=Cov_prior_parms
Cov.prior.fixed=Cov_prior_fixed
Cov.prior.n=Cov_prior_n
n.species=n_species
n.obs.cov=0
spat.ind=spat_ind
Inits=list(tau.eps=c(2,2,2,2))
grps=TRUE
Prior.pars=Prior_pars
True.species=True_sp
Omega.true=Omega_true
Eta.true=Eta_true
Alpha.true=Alpha_true
DEBUG=TRUE


set.seed(12346)
#Control$start.file=NULL
MCMC=hierarchical_boss_st(Dat=Dat,K=K,Area.hab=Area_hab,Area.trans=Area_trans,Mapping=Mapping,DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Hab.cov=Hab_cov,Obs.cov=NULL,Hab.formula=hab_formula,Cov.prior.pdf=Cov_prior_pdf,Cov.prior.parms=Cov_prior_parms,Cov.prior.fixed=Cov_prior_fixed,Cov.prior.n=Cov_prior_n,n.species=n_species,n.obs.cov=0,spat.ind=spat_ind,Psi=Psi,Inits=Inits,grps=TRUE,Control=Control,Prior.pars=Prior_pars,Surveyed=Surveyed,True.species=True_sp,Omega.true=Omega_true,Eta.true=Eta_true,Alpha.true=Alpha_true,DEBUG=TRUE)
save(MCMC,file="BOSS2012_MCMC_output2.RData")

load('./AlaskaBeringData2012_2013_14Dec2015.Rdat')
PostN=cat_preds(fname=Control$fname,n.files=Control$n.files)
Cur.G=matrix(apply(PostN[2,1:80,],2,'median'),S,t_steps)
#for(i in 1:t.steps){
plot_N_map(20,Cur.G,Grid=Data$Grid$y2012)

#plot from workspace
load('BOSS2012_mcmc_workspace.rdata')
Cur.G=matrix(Par$G[4,],Meta$S,Meta$t.steps)
plot_N_map(1,Cur.G,Grid=Data$Grid$y2012)

