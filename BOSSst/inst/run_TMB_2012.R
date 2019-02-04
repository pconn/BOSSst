# Analyze seal abundance from aerial survey flights in the eastern Bering Sea in spring of 2012
library( TMB )
#library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

data(Data2012)

# Compile
TmbFile = "./BOSSst/src/BossTMB_noST_tweedie"  #adjust if needed based on working directory
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

n_cells = length(Data2012$Area_hab)
t_steps = 29

#a few more data adjustments; don't model counts in cells where < 0.0001 of grid cell is sampled
area_thresh = 0.0002
Which_counts_model = which(Data2012$Area_trans>area_thresh)

# Settings
n_species = 4
n_samp = nrow(Data2012$DayHour)
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category


Thin_zero = rep(0.9,100)
Thin = c(Data2012$Det_priors$p_sd,Thin_zero,Data2012$Det_priors$p_rn,Thin_zero,Data2012$Det_priors$p_bd,Thin_zero,Data2012$Det_priors$p_rd,Thin_zero)

Sigma_thin = vector("list",n_species)
Sigma_thin[[1]]=Data2012$Det_priors$Var_sd
Sigma_thin[[2]]=diag(100)*0.0001
Sigma_thin[[3]]=Data2012$Det_priors$Var_rn
Sigma_thin[[4]]=diag(100)*0.0001
Sigma_thin[[5]]=Data2012$Det_priors$Var_bd
Sigma_thin[[6]]=diag(100)*0.0001
Sigma_thin[[7]]=Data2012$Det_priors$Var_sd  #set variance for ringed = the variance for spotted
Sigma_thin[[8]]=diag(100)*0.0001
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
MisID_zero_cols = c(3,6,9,12)  #so matrix entries (1,3), (2,6), (3,9), (4,12) will be fixed to zero on multinomial logit scale
n_misID_par=length(MisID_pos_rows)
MisID_pos_cols = rep(c(1,2,4,5,7,8,10,11,13),4) 
MisID_pars = Data2012$Beta_psi
MisID_Sigma = Data2012$VC_psi

#assemble hot spot data into count data by cell surveyed
C_i = matrix(0,n_samp,n_obs_types)  #one entry per sample
Data2012$Dat[,"Group"]=as.numeric(Data2012$Dat[,"Group"])
for(i in 1:n_samp){
  Cur_dat = Data2012$Dat[which(Data2012$Dat$Transect==i),]
  if(nrow(Cur_dat)>0){
    for(j in 1:nrow(Cur_dat))C_i[i,Cur_dat[j,"Obs"]]=C_i[i,Cur_dat[j,"Obs"]]+Cur_dat[j,"Group"]
  }
}

# Options
Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?

#for extra day effect on availability...
X_day = matrix(0,n_samp,2)
X_day[,1] = Data2012$DayHour[,1]
X_day[,2] = (X_day[,1])^2 
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
X_day = kronecker(diag(n_species),X_day)

S_i = (Data2012$DayHour[,1]-1)*n_cells+Data2012$Mapping[,1]

X_s = Data2012$Hab_cov[,c(2,4,5,6,7,8,10,11,12)]
#Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Area_trans,"A_s"=rep(Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species)

#take cubed root of depth
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2

Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Data2012$Area_trans,"A_s"=rep(Data2012$Area_hab,t_steps),"S_i"=S_i-1,"X_s"=as.matrix(X_s),"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=Matrix(Sigma_logit_thin),"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"n_s" = n_cells, "n_sp" = n_species,"Which_counts_model"=Which_counts_model-1,"Which_obs_sp"=c(0,0,0,1,1,1,2,2,2,3,3,3,1),"n_i_real"=n_samp-100,"h_mean"=c(mean(Data2012$Det_priors$p_sd),mean(Data2012$Det_priors$p_rn),mean(Data2012$Det_priors$p_bd),mean(Data2012$Det_priors$p_rd)))


# Parameters / Initial values 
load('c:/users/paul.conn/git/BOSSst/tweedie_noObsVar_dayEff.RData')
Report=Out$Report
#Note: Initial parameter values obtained from fitting a simpler model without observation variance
Params = list("log_N"=log(c(300000,60000,300000,300000)),"Beta"=matrix(0,4,9),"thin_beta_day"=c(0,0,.5,.6,-.3,.3,1,-1),"phi_log"=log(c(3.6,1.6,1.6,1.7)),"p_logit"=c(-1,-2,-2,2),"thin_logit_i"=Thin_logit,"MisID_pars"=MisID_pars)
load('c:/users/paul.conn/git/BOSSst/tweedie_ObsVar_noDay_2012.RData')

# Random - declare observation errors as random effects to be integrated out with Laplace approx
Random = c("thin_logit_i")

# Fix MisID parameters to initial values
Map = list()
Map[["MisID_pars"]]=factor(rep(NA,length(Params$MisID_pars)))

# Make object
dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
#image(Obj$env$spHess(random=TRUE))  #look at covariance structure

# Run
start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
Report = Obj$report()

Sys.time() - start.time

SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )

Out = list(Report=Report,SD=SD)
save(Out,file="tweedie_ObsVar_2012.RData")
