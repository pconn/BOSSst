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

load('BOSS_data_TMB_2012.Rda')
load('detection_priors.RData')

# Compile
setwd("C:/Users/paul.conn/git/BOSSst/")

TmbFile = "C:/Users/paul.conn/git/BOSSst/BOSSst/src/BossTMB_spatial_fixObs_tweedie"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

#TmbFile = "C:/Users/paul.conn/git/BOSSst/BOSSst/src/BossTMB_noST"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
#compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

n_cells = length(Area_hab)
t_steps = 29
#require(STabundance)
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
#source('./OkhotskSeal/R/sim_funcs.R')
#source('./OkhotskSeal/R/sim_data_generic.R')

#a few more data adjustments; don't model counts in cells where < 0.0001 of grid cell is sampled
area_thresh = 0.0002
Which_counts_model = which(Area_trans>area_thresh)

# Settings
n_species = 4
n_samp = nrow(DayHour)
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category


#misID_fix = c(10,10,10)  #which column of confusion matrix to fix to 1.0 on the multinomial logit scale for each species
Thin_zero = rep(0.9,100)
Thin = c(Det_priors$p_sd,Thin_zero,Det_priors$p_rn,Thin_zero,Det_priors$p_bd,Thin_zero,Det_priors$p_rd,Thin_zero)

Sigma_thin = vector("list",n_species)
Sigma_thin[[1]]=Det_priors$Var_sd
Sigma_thin[[2]]=diag(100)*0.0001
Sigma_thin[[3]]=Det_priors$Var_rn
Sigma_thin[[4]]=diag(100)*0.0001
Sigma_thin[[5]]=Det_priors$Var_bd
Sigma_thin[[6]]=diag(100)*0.0001
Sigma_thin[[7]]=Det_priors$Var_sd  #set variance for ringed = the variance for spotted
Sigma_thin[[8]]=diag(100)*0.0001
#convert into block diagonal 
Sigma_thin = as.matrix(bdiag(Sigma_thin))
Sigma_thin = as(Sigma_thin,"dgTMatrix")

#compute logit scale distribution using delta method
diff_logit <- function(x) exp(x)/(1+exp(x))^2
Thin_logit = log(Thin/(1-Thin))
Diff_logit = diag(diff_logit(Thin))
Sigma_logit_thin = Diff_logit %*% Sigma_thin %*% t(Diff_logit)

#set up spatial model
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:

x.min=min(Coords[,1])-100000
x.max=max(Coords[,1])+100000
y.min=min(Coords[,2])-100000
y.max=max(Coords[,2])+100000
X=x.min+(x.max-x.min)/7*c(0:7)
Y=y.min+(y.max-y.min)/7*c(0:7)
XY=expand.grid(x=X,y=Y)
library(sp)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0","+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Knots=SpatialPoints(coords=XY,proj4string=CRS(laea_180_proj))
Distances=matrix(0,n_cells,length(Knots))
Coord_knots=coordinates(Knots)
for(icell in 1:n_cells){
  for(iknot in 1:length(Knots)){
    Distances[icell,iknot]=sqrt(sum((Coord_knots[iknot,]-Coords[icell,])^2))
  }
}
Min_dist = apply(Distances,2,'min')
my.buffer=150000
Which.include=which(Min_dist<my.buffer)
Knot_loc=Coord_knots[Which.include,]
Knot.cell.distances=Distances[,Which.include]
mesh = inla.mesh.create( Knot_loc )
n_knots = mesh$n


Knot.cell.distances = matrix(0,n_cells,n_knots)  #need to recalculate because delauney set of knots larger
Kmap = 0*Knot.cell.distances
for(i in 1:n_cells){
  for(j in 1:n_knots)Knot.cell.distances[i,j]=sqrt(sum((Coords[i,]-mesh$loc[j,1:2])^2))
  second_largest = sort(Knot.cell.distances[i,])[2]  #using inverse distance weighting of closest three knots hampers sparsity
  Which_largest = which(Knot.cell.distances[i,]<=second_largest)
  Kmap[i,Which_largest]=1/Knot.cell.distances[i,Which_largest]
  Kmap[i,]=Kmap[i,]/sum(Kmap[i,])
}
#now need to redo mesh with standardized distances; 1 = distance between cells
cell_dist = Coords[2,1]-Coords[1,1]
Knot_loc2 = Knot_loc/cell_dist
mesh = inla.mesh.create( Knot_loc2 )
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]


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
X_day = kronecker(diag(n_species),X_day)

S_i = (DayHour[,1]-1)*n_cells+Mapping[,1]

X_s = Hab_cov[,c(2,4,5,6,7,8,10,11,12)]
#Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species)

#take cubed root of depth
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2

Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"spde"=spde,"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"Kmap"=Kmap,"n_s" = n_cells, "n_sp" = n_species,"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars,"Which_counts_model"=Which_counts_model-1,"Which_obs_sp"=c(0,0,0,1,1,1,2,2,2,3,3,3,1))


# Parameters / Initial values - set close to truth for faster convergence
Beta_init = matrix(0,n_species,ncol(Data$X_s))
Etainput_s = matrix(0,n_species,mesh$n)  #spatio-temporal random effects
#Params = list("log_N"=log(c(200000,80000,100000,100000)),"Beta"=Beta_init,"thin_logit_i"=Thin_logit,"thin_beta_day"=rep(0,2*n_species),"MisID_pars"=MisID_pars)
Params = list("log_N"=log(c(200000,80000,100000,100000)),"Beta"=Beta_init,"logtau_eta"=rep(0,n_species),"logit_kappa_eta"=rep(-1,n_species),"Etainput_s"=Etainput_s,"thin_beta_day"=rep(0,2*n_species),"phi_log"=rep(-1,4),"p_logit"=rep(0,4))

# Random
Random = c( "Etainput_s") 

# Fix parameters
Map = list()
Map[["thin_beta_day"]]=factor(rep(NA,length(Params$thin_beta_day)))


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
#Opt[["diagnostics"]] = data.frame( "Paramrm(list="=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report()

cat(Sys.time() - start.time)
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
plot_N_map_xy(N=Report$Z_s[3,1:n_cells],XY=Coords2,leg.title="Abundance")
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
plot(DayHour[,1],Report$Thin_trans[2,])

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

