# Abundance estimation for multiple species from count data 
# using spatial regression with prior distributions on detection probability at
# each location sampled.  In this version counts are disassociated from species; a
# prior distribution on confusion matrix parameters provides the link to species-specific 
# counts.

library( TMB )
#library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

data(Data2013)

# Compile
#setwd("C:/Users/paul.conn/git/BOSSst/")

TmbFile = "./BOSSst/src/BossSt"  #note: edit depending on your working directory
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


n_cells = length(Data2013$Area_hab)
t_steps = 30

#a few more data adjustments; don't model counts in cells where < 0.0001 of grid cell is sampled
area_thresh = 0.0002
Which_counts_model = which(Data2013$Area_trans>area_thresh)

# Settings
n_species = 4
n_samp = nrow(Data2013$DayHour)
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category
n_zeroes = sum(Data2013$Area_trans==0.9)  #these are 'pseudo-zeroes' - zero counts put in some places where ice = 0 to tell model there are no seals hauled out in open water

Thin_zero = rep(0.9,n_zeroes)
Thin = c(Data2013$Det_priors$p_sd,Thin_zero,Data2013$Det_priors$p_rn,Thin_zero,Data2013$Det_priors$p_bd,Thin_zero,Data2013$Det_priors$p_rd,Thin_zero)

Sigma_thin = vector("list",n_species)
Sigma_thin[[1]]=Data2013$Det_priors$Var_sd
Sigma_thin[[2]]=diag(n_zeroes)*0.0001
Sigma_thin[[3]]=Data2013$Det_priors$Var_rn
Sigma_thin[[4]]=diag(n_zeroes)*0.0001
Sigma_thin[[5]]=Data2013$Det_priors$Var_bd
Sigma_thin[[6]]=diag(n_zeroes)*0.0001
Sigma_thin[[7]]=Data2013$Det_priors$Var_sd  #set variance for ringed = the variance for spotted
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
MisID_pars = Data2013$Beta_psi
MisID_Sigma = Data2013$VC_psi

#assemble hot spot data into count data by cell surveyed
Dat = Data2013$Dat
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
X_day[,1] = Data2013$DayHour[,1]
X_day[,2] = (X_day[,1])^2 
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
X_day = kronecker(diag(n_species),X_day)

S_i = (Data2013$DayHour[,1]-1)*n_cells+Data2013$Mapping[,1]

X_s = Data2013$Hab_cov[,c(2,4,5,6,7,8,10,11,12)]

#take cubed root of depth
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2


Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Data2013$Area_trans,"A_s"=rep(Data2013$Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species,"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars,"Which_counts_model"=Which_counts_model-1,"Which_obs_sp"=c(0,0,0,1,1,1,2,2,2,3,3,3,1),"n_i_real"=n_samp-n_zeroes,"h_mean"=c(mean(Data2013$Det_priors$p_sd),mean(Data2013$Det_priors$p_rn),mean(Data2013$Det_priors$p_bd),mean(Data2013$Det_priors$p_rd)))

Beta_init = matrix(0,n_species,ncol(Data$X_s))
Params = list("log_N"=log(c(160000,40000,220000,180000)),"Beta"=Beta_init,"thin_beta_day"=c(1,-1,-.3,-.4,0,-.2,2,-1),"phi_log"=rep(1,4),"p_logit"=rep(-1,4),"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars)
#note: these initial values are based loosely on a simpler model fitted to data without variance components for thin_logit_i, MisID_pars
#in general this is a good strategy: start with a simple model, and use estimates as initial values for a more complex model, then use those estimates as initial values for a more complex one, etc.

# Random
Random = c("thin_logit_i","MisID_pars")

# Fix parameters
Map = list()
#Map[["thin_beta_day"]]=factor(rep(NA,length(Params$thin_beta_day)))  #uncomment to fix day effect = 0



# Make object
dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
#image(Obj$env$spHess(random=TRUE))  #look at covariance structure

# Perform estimation
start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Report = Obj$report()

Sys.time() - start.time

SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
Out = list(Report=Report,SD=SD)
save(Out,file="tweedie_ObsVar_Day_2013.RData")

