# Abundance estimation for multiple species from count data 
# using spatial regression with prior distributions on detection probability at
# each location sampled.  In this version counts are disassociated from species; a
# prior distribution on confusion matrix parameters provides the link to species-specific 
# counts.
# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# Eta -- spatial random effects
# Beta -- spatial regression parameters
# logtau -- precision of spatial random effects

#library( RandomFields )
library( TMB )
library( INLA )
library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
#library( mvtnorm)
#library(TMBdebug)

load('BOSS_data_TMB_2013.Rda')
load('detection_priors2013.RData')

# Compile
setwd("C:/Users/paul.conn/git/BOSSst/")

TmbFile = "C:/Users/paul.conn/git/BOSSst/BOSSst/src/BossTMB_noST_fixObs_tweedie"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

#TmbFile = "C:/Users/paul.conn/git/BOSSst/BOSSst/src/BossTMB_noST"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
#compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

n_cells = length(Area_hab)
t_steps = 30
#require(STabundance)
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
#source('./OkhotskSeal/R/sim_funcs.R')
#source('./OkhotskSeal/R/sim_data_generic.R')

#a few more data adjustments; don't model counts in cells where < 0.0002 of grid cell is sampled
area_thresh = 0.0002
Which_counts_model = which(Area_trans>area_thresh)

# Settings
n_species = 4
n_samp = nrow(DayHour)
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category
n_zeroes = sum(Area_trans==0.9)

#misID_fix = c(10,10,10)  #which column of confusion matrix to fix to 1.0 on the multinomial logit scale for each species
Thin_zero = rep(0.9,n_zeroes)
Thin = c(Det_priors$p_sd,Thin_zero,Det_priors$p_rn,Thin_zero,Det_priors$p_bd,Thin_zero,Det_priors$p_rd,Thin_zero)

Sigma_thin = vector("list",n_species)
Sigma_thin[[1]]=Det_priors$Var_sd
Sigma_thin[[2]]=diag(n_zeroes)*0.0001
Sigma_thin[[3]]=Det_priors$Var_rn
Sigma_thin[[4]]=diag(n_zeroes)*0.0001
Sigma_thin[[5]]=Det_priors$Var_bd
Sigma_thin[[6]]=diag(n_zeroes)*0.0001
Sigma_thin[[7]]=Det_priors$Var_sd  #set variance for ringed = the variance for spotted
Sigma_thin[[8]]=diag(n_zeroes)*0.0001
#convert into block diagonal 
Sigma_thin = as.matrix(bdiag(Sigma_thin))
Sigma_thin = as(Sigma_thin,"dgTMatrix")

#compute logit scale distribution using delta method
diff_logit <- function(x) -1/(x*(x-1))
Thin_logit = log(Thin/(1-Thin))
Diff_logit = diag(diff_logit(Thin))
Sigma_logit_thin = Diff_logit %*% Sigma_thin %*% t(Diff_logit)



#set up misID matrix using parameters on multinomial-logit scale
misIDcols = c(1:n_obs_types)
MisID_pos_rows = c(rep(1,9),rep(2,9),rep(3,9),rep(4,9))
MisID_zero_cols = c(3,6,9,12)
n_misID_par=length(MisID_pos_rows)
MisID_pos_cols = rep(c(1,2,4,5,7,8,10,11,13),4) 
MisID_pars = Beta_psi
# MisID_real = matrix(-20,n_species,n_obs_types)
# for(ipar in 1:n_misID_par)MisID_real[MisID_pos_rows[ipar],MisID_pos_cols[ipar]]=MisID_pars[ipar]
# for(isp in 1:n_species){
#   MisID_real[isp,MisID_zero_cols[isp]]=0
#   MisID_real[isp,] = exp(MisID_real[isp,])/sum(exp(MisID_real[isp,]))
# }
MisID_Sigma = VC_psi

#assemble hot spot data into count data by cell surveyed
C_i = matrix(0,n_samp,n_obs_types)  #one entry per sample
Dat[,"Group"]=as.numeric(Dat[,"Group"])
for(i in 1:n_samp){
  Cur_dat = Dat[which(Dat$Transect==i),]
  if(nrow(Cur_dat)>0){
    for(j in 1:nrow(Cur_dat))C_i[i,Cur_dat[j,"Obs"]]=C_i[i,Cur_dat[j,"Obs"]]+Cur_dat[j,"Group"]
  }
}

# Options
Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?

#for extra day effect on availability...
X_day = matrix(0,n_samp,2)
X_day[,1] = DayHour[,1]
X_day[,2] = (X_day[,1])^2 
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
X_day[(n_samp-n_zeroes+1):n_samp,]=0  #no day effects on pseudo-zeroes
X_day = kronecker(diag(n_species),X_day)

S_i = (DayHour[,1]-1)*n_cells+Mapping[,1]

X_s = Hab_cov[,c(2,4,5,6,7,8,10,11,12)]
#Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species)

#take cubed root of depth
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2


Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_logit_i"=Thin_logit,"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_pars"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species,"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars,"Which_counts_model"=Which_counts_model-1,"Which_obs_sp"=c(0,0,0,1,1,1,2,2,2,3,3,3,1),"n_i_real"=n_samp-n_zeroes,"h_mean"=c(mean(Det_priors$p_sd),mean(Det_priors$p_rn),mean(Det_priors$p_bd),mean(Det_priors$p_rd)))
save(Data,file="Data2013.RData")

Beta_init = matrix(0,n_species,ncol(Data$X_s))
# Parameters / Initial values - set close to truth for faster convergence
#load('tweedie_2013_noObsVar_day.RData')
#Params = list("log_N"=log(Report$N),"Beta"=Report$Beta,"thin_beta_day"=0*Report$thin_beta_day,"phi_log"=log(Report$phi),"p_logit"=log((Report$power-1)/(1-(Report$power-1))),"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars)

Params = list("log_N"=log(c(150000,50000,200000,150000)),"Beta"=Beta_init,"thin_beta_day"=rep(0,2*n_species),"phi_log"=rep(0.5,4),"p_logit"=rep(-2,4))

# Random
Random = NULL
#Random = c("thin_logit_i","MisID_pars")

# Fix parameters
Map = list()
Map[["thin_beta_day"]]=factor(rep(NA,length(Params$thin_beta_day)))

#Map[["p_logit"]]=factor(rep(NA,length(Params$p_logit)))
#Map[["phi_log"]]=factor(rep(NA,length(Params$phi_log)))


# Make object
#compile( paste0(Version,".cpp") )
dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
#setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
#image(Obj$env$spHess(random=TRUE))  #look at covariance structure

# Run
#Lower = -Inf
#Upper = Inf
start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
#Opt = Optimize(obj=Obj,newtonsteps=2,lower=-50,upper=50,loopnum=1,bias.correct=TRUE)
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
#Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report()

Sys.time() - start.time

SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )

Out = list(Report=Report,SD=SD)
save(Out,file="tweedie_noObsVar_noDay_2013.RData")


# Potentially fix random fields with zero sample or population variance
#if( any(Report$MargSD_z<0.001) ){
#Which = which(Report$MargSD_z<0.001)
# Map[["logtau_z"]] = factor( ifelse(1:2==Which,NA,1:2) )
#  if(length(Which)==2){
#    Map[["logkappa_z"]] = factor( c(NA,NA) )
# }
# Map[["beta_b"]] = factor(NA)
#  Params[["beta_b"]] = 0
#  Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
#  Params[["etainput_s"]][] = 0
#  # Re-run
#  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
#  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
#  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
#}

Coords2=Coords
Coords2[,1] = round(Coords2[,1])
Coords2[,2] = round(Coords2[,2])
plot_N_map_xy(N=Report$Z_s[1,1:n_cells],XY=Coords2,leg.title="Abundance")
#cur_day=10
#plot_N_map_xy(N=Report$Z_s[1,((cur_day-1)*n_cells+1):(cur_day*n_cells)],XY=loc_s,leg.title="Abundance")

Converge=Opt$convergence


#plot count data by apparent species
Count = rowSums(Data$C_i[,1:3])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,4:6])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,7:9])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,10:12])
crap = Count/Data$P_i
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=crap,XY=Coords3,leg.title="Count")
# 

Count = rowSums(Report$E_count_obs[,4:6])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = Report$E_count_sp[2,]
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

# 
#plot thinning by day
par(mfrow=c(2,2))
sp=1
plot(DayHour[1:Data$n_i_real,1],Report$Thin_trans[sp,1:Data$n_i_real],xlab='Survey day',ylab='Detection probability')
points(DayHour[1:Data$n_i_real,1],plogis(Data$thin_mu_logit[(sp-1)*n_samp+c(1:Data$n_i_real)]),col='red')
sp=2
plot(DayHour[1:Data$n_i_real,1],Report$Thin_trans[sp,1:Data$n_i_real],xlab='Survey day',ylab='Detection probability')
points(DayHour[1:Data$n_i_real,1],plogis(Data$thin_mu_logit[(sp-1)*n_samp+c(1:Data$n_i_real)]),col='red')
sp=3
plot(DayHour[1:Data$n_i_real,1],Report$Thin_trans[sp,1:Data$n_i_real],xlab='Survey day',ylab='Detection probability')
points(DayHour[1:Data$n_i_real,1],plogis(Data$thin_mu_logit[(sp-1)*n_samp+c(1:Data$n_i_real)]),col='red')
sp=4
plot(DayHour[1:Data$n_i_real,1],Report$Thin_trans[sp,1:Data$n_i_real],xlab='Survey day',ylab='Detection probability')
points(DayHour[1:Data$n_i_real,1],plogis(Data$thin_mu_logit[(sp-1)*n_samp+c(1:Data$n_i_real)]),col='red')

#look at proportion of zeroes
library(tweedie)
get_pred_zeros_tweedie <- function(E_count,phi,power,sim=TRUE){
  if(sim){
    Rand = matrix(0,nrow(E_count),ncol(E_count))
    for(i in 1:4)Rand[,i]=rtweedie(nrow(E_count),mu=E_count[,i],phi=rep(phi[i],nrow(E_count)),power=rep(power[i]))
  }
  if(!sim)Rand = E_count
  Prop_zero=rep(0,ncol(E_count))
  for(i in 1:ncol(E_count))Prop_zero[i]=sum(Rand[,i]==0)
  Prop_zero
}

E_count = cbind(rowSums(Report$E_count_obs[,1:3]),rowSums(Report$E_count_obs[,4:6]),rowSums(Report$E_count_obs[,7:9]),rowSums(Report$E_count_obs[,10:12]))
Data_sp = cbind(rowSums(Data$C_i[,1:3]),rowSums(Data$C_i[,4:6]),rowSums(Data$C_i[,7:9]),rowSums(Data$C_i[,10:12]))
get_pred_zeros_tweedie(E_count=E_count,phi=Report$phi,power=Report$power)
get_pred_zeros_tweedie(Data_sp,sim=FALSE)

par(mfrow=c(2,2))
plot(E_count[,1],Data_sp[,1],xlab="Expected count",ylab="Observed count")
abline(0,1)
plot(E_count[,2],Data_sp[,2],xlab="Expected count",ylab="Observed count")
abline(0,1)
plot(E_count[,3],Data_sp[,3],xlab="Expected count",ylab="Observed count")
abline(0,1)
plot(E_count[,4],Data_sp[,4],xlab="Expected count",ylab="Observed count")
abline(0,1)

#randomized quantile residuals
library(tweedie)
Resids = 0*Report$E_count_obs
for(icol in 1:13){
  Resids[,icol] = ptweedie(Data$C_i[,icol]-1,mu=Report$E_count_obs[,icol],phi=Report$phi[Data$Which_obs_sp[icol]+1],power=Report$power[Data$Which_obs_sp[icol]+1])+runif(nrow(Data$C_i))*dtweedie(Data$C_i[,icol],mu=Report$E_count_obs[,icol],phi=Report$phi[Data$Which_obs_sp[icol]+1],power=Report$power[Data$Which_obs_sp[icol]+1])
}

Resids[Resids>1]=0.999

Resid_binned = matrix(0,20,13)
for(irow in 1:nrow(Resids)){
  for(icol in 1:13){
    Resid_binned[ceiling(Resids[irow,icol]/0.05),icol]=Resid_binned[ceiling(Resids[irow,icol]/0.05),icol]+1
  }
}
Xsq = rep(0,13)
for(i in 1:13){
  Xsq[i]=20/nrow(Resids)*sum((Resid_binned[,i]-nrow(Resids)/20)^2)
}
Pval = 1-pchisq(Xsq,df=19)   #chi square p-value for each bin

Resids.df = data.frame(Residual = as.vector(Resids))
Labels1 = rep('',13)
for(i in 1:13){
  Labels1[i] = paste0('Obs = ',i,', p=',format(Pval[i],digits=2))
}
Resids.df$Labels=rep(Labels1,each=nrow(Resids))

ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:20)/20)+facet_wrap(~Labels)



# plots of spatial residuals
#plot count data by apparent species
Map.list = vector("list",13)
Coords3 = Coords2[Mapping[,1],]
for(i in 1:13){
  Map.list[[i]]=plot_N_map_xy(N=Resids[,i],XY=Coords3,leg.title="Randomized residual")
}

