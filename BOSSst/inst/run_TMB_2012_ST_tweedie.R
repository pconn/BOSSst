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

TmbFile = "C:/Users/paul.conn/git/BOSSst/BOSSst/src/BossTMB_fixObs_eps"
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

# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
mesh = inla.mesh.create( Knot_loc )

n_knots = mesh$n

Knot.cell.distances = matrix(0,n_cells,n_knots)
Kmap = 0*Knot.cell.distances
for(i in 1:n_cells){
  for(j in 1:n_knots)Knot.cell.distances[i,j]=sqrt(sum((Coords[i,]-mesh$loc[j,1:2])^2))
  #cur.which = which(Knot.cell.distances[i,]==min(Knot.cell.distances[i,]))
  #Kmap[i,cur.which]=1
  #third_largest = sort(Knot.cell.distances[i,])[3]  #using inverse distance weighting of closest three knots hampers sparsity
  #Which_largest = which(Knot.cell.distances[i,]<=third_largest)
  #Kmap[i,Which_largest]=1/Knot.cell.distances[i,Which_largest]
  ##Kmap[i,(nrow(loc_k)+1):n_knots]=0       #take out influence of knots on edges
  #Kmap[i,]=Kmap[i,]/sum(Kmap[i,])
  second_largest = sort(Knot.cell.distances[i,])[2]  #using inverse distance weighting of closest three knots hampers sparsity
  Which_largest = which(Knot.cell.distances[i,]<=second_largest)
  Kmap[i,Which_largest]=1/Knot.cell.distances[i,Which_largest]
  ##Kmap[i,(nrow(loc_k)+1):n_knots]=0       #take out influence of knots on edges
  Kmap[i,]=Kmap[i,]/sum(Kmap[i,])
}

#now need to redo mesh with standardized distances; 1 = distance between cells
cell_dist = Coords[2,1]-Coords[1,1]
Knot_loc2 = Knot_loc/cell_dist
mesh = inla.mesh.create( Knot_loc2 )


spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]



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

#X_s = Hab_cov[,c(5,11)]
#take cubed root of depth
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2

Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"spde"=spde,"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"Kmap"=Kmap,"n_s" = n_cells, "n_sp" = n_species,"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars,"Which_counts_model"=Which_counts_model-1)


# Parameters / Initial values - set close to truth for faster convergence
Beta_init = matrix(0,n_species,ncol(Data$X_s))
Etainput_st = array(0,dim=c(n_species,mesh$n,t_steps))  #spatio-temporal random effects
RE_tot = matrix(0,n_species,n_cells*t_steps)
#Params = list("log_N"=log(c(200000,80000,100000,100000)),"Beta"=Beta_init,"thin_logit_i"=Thin_logit,"thin_beta_day"=rep(0,2*n_species),"MisID_pars"=MisID_pars)
Params = list("log_N"=log(c(200000,80000,100000,100000)),"Beta"=Beta_init,"logtau_eta"=rep(0,n_species), "logkappa_eta"=rep(-2,n_species),"logit_rho"=rep(2,n_species),"Etainput_st"=Etainput_st,"thin_beta_day"=rep(0,2*n_species),"RE_tot"=RE_tot,"log_sd"=rep(-2,n_species))

# Random
Random = c( "Etainput_st") #"thin_logit_i")
#if(Use_REML==TRUE) Random = c(Random,"Beta")

# Fix parameters
Map = list(thin_beta_day=factor(rep(NA,length(Params$thin_beta_day))),logkappa_eta = factor(rep(NA,4)))


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
it=15
plot_N_map_xy(N=Report$Z_s[4,n_cells*(it-1)+c(1:n_cells)],XY=Coords2,leg.title="Abundance")
plot_N_map_xy(N=Kmap%*%Report$Eta_st[3,,1],XY=Coords2,leg.title="Abundance")
#cur_day=10
#plot_N_map_xy(N=Report$Z_s[1,((cur_day-1)*n_cells+1):(cur_day*n_cells)],XY=loc_s,leg.title="Abundance")

Converge=Opt$convergence


#plot count data by apparent species
Count = rowSums(Data$C_i[,1:3])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")
plot_N_map_xy(N=Report$E_count_sp[1,],XY=Coords3,leg.title="Count")

#look at thinning parameters
plot_N_map_xy(N=Thin[1:1285],XY=Coords3,leg.title="Count")
plot_N_map_xy(N=Report$Thin_trans[1,],XY=Coords3,leg.title="Count")

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
get_pred_zeros <- function(E_count,sim=TRUE){
  if(sim)Rand = matrix(rpois(length(E_count),E_count),nrow(E_count),ncol(E_count))
  if(!sim)Rand = E_count
  Prop_zero=rep(0,ncol(E_count))
  for(i in 1:ncol(E_count))Prop_zero[i]=sum(Rand[,i]==0)
  Prop_zero
}

E_count = cbind(rowSums(Report$E_count_obs[,1:3]),rowSums(Report$E_count_obs[,4:6]),rowSums(Report$E_count_obs[,7:9]),rowSums(Report$E_count_obs[,10:12]))
Data_sp = cbind(rowSums(Data$C_i[,1:3]),rowSums(Data$C_i[,4:6]),rowSums(Data$C_i[,7:9]),rowSums(Data$C_i[,10:12]))
get_pred_zeros(E_count=E_count)
get_pred_zeros(Data_sp,sim=FALSE)


