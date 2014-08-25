#' Primary function for hierarchical, areal analysis of distance sampling data.  This function
#' pre-processes data and calls other functions to perform the analysis, and is the only function
#' the user needs to call themselves. 
#'
#' @param Dat 	A data frame with the following columns:
#' 		(1)numeric transect ID (1 to # of transects)
#' 		(2)photo obtained (0/1) ? 
#' 		(3) Observation type (integer - the max integer being 'unknown' if applicable) [NOTE: modeled as factor, but need to be input as integers to account for unknown species observations]
#' 		(4-x)(Observer covariates); (things like survey conditions or observer skill that affect misID; things that don't change during a transect.  Note these also need to be provided in Obs.cov)
#' 		(x+1-??)(Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+1); could also include observer certainty here
#' 		Note that column names can be used to tag covariates, and that object types (e.g. numeric, factor) will be preserved in analysis (note all covariates must be categorical right now!)
#' @param K    An S * n.knots matrix holding kernel weights for process convolution (can be sparse)
#' @param Area.hab   A vector giving the area of each geographical strata (default is equal area)
#' @param Mapping  A two column matrix where the first column gives the cell id number for each transect and the second gives the time step where observation was made
#' @param Area.trans	A vector giving the effective area covered by each transect as fraction of total area in the strata it is located
#' @param DayHour A (n.transect X 2) matrix providing row and column entries into the Thin array. Each row corresponds to an entry in Mapping
#' @param Thin An (n.species X n.days X n.hours X n.iter) array providing n.iter posterior samples of the thinning parameters
#' @param Prop.photo A vector giving the proportion of of the area sampled in each transect that is photographed
#' @param n.species An integer giving the true number of species
#' @param n.obs.cov	Number of observer covariates (e.g., visibility, etc.)
#' @param Hab.cov	A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates; # rows = S*T
#' @param Obs.cov  A (# of transects X # of observer covariates) size matrix giving observer covariate values for each transect flown
#' @param Hab.formula	A list of formulas giving the specific model for latent log abundance multinomial cell probabilities at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param Cov.prior.pdf	If individual covariates are provided, this character matrix gives the form of the prior pdfs for each species & covariate
#'		  current possibilities are "poisson", "pois1","poisson_ln","pois1_ln",uniform.disc","multinom","uniform.cont", or "normal".
#'		  "pois1" is 1+x where x~poisson; "poisson_ln" and "pois1_ln" are lognormal poisson models that incorporate overdispersion. 
#' @param Cov.prior.parms	A (s X k X n) array where s is the number of species, n is the number of individual covariates (other than distance), and
#' 		k is the maximum number of parameters considered for a single covariate (NAs can be used to fill this matrix
#'      out for covariate priors that have <k parameters).  If Cov.prior.fixed=1 for a given entry, the prior parameters supplied
#'      in each column apply to the prior pdf itself, and are treated as fixed.  If Cov.prior.fixed=0, the model will attempt
#'  	to estimate the posterior distribution of model parameters, given hyperpriors.  In this case, it is actually the hyperpriors
#'      that are being specified.  For "poisson", and "pois1", it is assumed that lambda~gamma(alpha,beta), so alpha
#' 		and beta must be supplied.  For "poisson_ln", and "pois1_ln", the model is lambda_i=exp(-sigma*Z_i+theta), so it is priors
#' 		for theta and sigma that are specified (in that order).  Theta is assumed to have a normal(mu,s^2) distribution,
#' 		and sigma is assumed to have a uniform(0,a) distribution; thus, priors are specified for these models as (mu,s, and a).
#' 		For the multinomial pdf, prior parameters of the dirichlet distribution must be specified if Cov.prior.fixed=1.
#' @param Cov.prior.fixed  An indicator matrix specifying which (if any) individual covariate distributions should be fixed during estimation
#' @param Cov.prior.n  An (# species X # indiv. covariates) matrix giving the number of parameters in each covariate pdf
#' @param spat.ind	If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE
#' @param Psi An array holding posterior samples of the confusion matrix (dim = #species,#obs,#mcmc.iter)
#' @param Thin An array holding posterior samples of thinning probabilities (dim = # species, # days surveyed, # hours in day, #iterations)
#' @param DayHour A matrix indicating which element of Thin each surveyed cell belongs to
#' @param grps 	If FALSE, detections are assumed to all be of individual animals
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, adapts MCMC proposals 
#'  "n.adapt": if adapt==TRUE, number of adapt iterations to employ (note: need to make burnin>n.adapt)
#'  "MH.omega" A matrix providing standard deviations for Langevin-Hastings proposals (rows: species; columns: number of time steps)
#'  "MH.N" vector of standard deviation for total abundance MH updates (1 for each species)
#'  "fix.tau.epsilon" If TRUE, fixes tau.epsilon to 100
#'  "species.optim" If TRUE, optimizes species updates [FALSE uses MH algorithm]
#' @param Inits	An (optional) list object providing initial values for model parameters, with the following objects:
#'  "hab": Initial values for habitat linear predictor parameters for poisson model;
#'	"G": Gives true group abundance (vector; 1 entry for each species)
#'	"tau.eta": If spat.ind==FALSE, precision for spatial ICAR model(s) for the Poisson component
#'	"tau.eps": Precision for exchangeable errors on classification probabilities 
#'  One need not specify an initial value for all parameter types (if less are specified, the others are generated randomly)
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'  "beta.tau" precision for Gaussian prior on regression parameters (default 0.1)
#'  "a.eps" alpha parameter for tau_epsilon~Gamma(alpha,beta) (default = 1.0)
#'  "b.eps" beta parameter for tau_epsilon~Gamma(alpha,beta) (default = 0.01)
#'  "a.eta" alpha param for gamma prior precision of spatio-temporal model
#'  "b.eta" beta param for gamma prior precision of spatio-temporal process
#'  "beta0.tau.rw2" prior precision for knot intercepts in rw2 model
#'  "beta1.tau.rw2" prior precision for knot slopes in rw2 model
#' @param post.loss If TRUE, calculates observed values and posterior predictions for detection data to use with posterior predictive loss functions
#' @param True.species if provided, this vector of true species values is used during estimation (for debugging only)
#' @param Omega.true If provided, this array gives true Omega values for debugging
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords areal model, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
#' @examples print("example analysis included in the script run_BOSS_sims.R")
hierarchical_boss_st<-function(Dat,K,Area.hab=1,Mapping,Area.trans,DayHour,Thin,Prop.photo=Prop.photo,Hab.cov,Obs.cov,Hab.formula,Cov.prior.pdf,Cov.prior.parms,Cov.prior.fixed,Cov.prior.n,n.species=1,n.obs.cov=0,spat.ind=FALSE,srr.tol=0.5,Psi,Inits=NULL,grps=FALSE,Control,adapt=TRUE,Prior.pars,post.loss=TRUE,True.species=NULL,Alpha.true=NULL,Omega.true=NULL,Eta.true=NULL,DEBUG=FALSE){
  require(mvtnorm)
	require(Matrix)
	require(truncnorm)
	require(mc2d)
	require(MCMCpack)
	
	S=nrow(K)
 	n.transects=length(Area.trans)
	n.ind.cov=ncol(Dat)-(3+n.obs.cov) #number of individual covariates 
  Obs.NArm=Dat[,3][-which(is.na(Dat[,3]))]
  n.obs.types=length(unique(Obs.NArm))
	
	#By no means exhaustive checking to make sure input values are internally consistent
  #More later...	
  
	#if(length(unique(Dat[,3]))>1)Dat[,3]=as.factor(Dat[,3])  #convert species to factors if not already
	cur.colnames=colnames(Dat)
  cur.colnames[1:3]=c("Transect","Photo","Species")
	if(grps==TRUE)cur.colnames[4+n.obs.cov]="Group"
	if(length(colnames(Dat))!=length(cur.colnames))cat("\n ERROR: mismatch between dimension of Dat and expected columns: check to make sure n.obs.cov, etc. correct")
	colnames(Dat)=cur.colnames
  
  CellTime=Mapping
  Mapping=(Mapping[,2]-1)*S+Mapping[,1]  #convert two dimensional mapping to one dimensional mapping
	t.steps=ceiling(max(Mapping)/S)  #currently assuming we're not projecting abundance ahead in time past the last survey
	
  if(is.null(True.species))Control$update.sp=TRUE
  else Control$update.sp=FALSE
  
  #convert character entries to factors
  for(i in 1:ncol(Dat))if(is.character(Dat[,i])&length(unique(Dat[,i]))>1)Dat[,i]=as.factor(Dat[,i])
  
  #for factors, determine levels, save labels, and convert to numeric
	factor.ind=sapply(Dat[1,],is.factor)
  which.factors=which(factor.ind==1)
  n.factors=sum(factor.ind)
  
  Factor.labels=vector("list",n.factors)
  if(n.factors>0){
    for(i in 1:n.factors){
      Factor.labels[[i]]=levels(Dat[,which.factors[i]])
    }
  }
  
	Dat.num=Dat
  #if(sum(Dat[,"Obs"]==0)>0)Dat.num[Dat[,"Obs"]==0,"Species"]=0  #if a missing obs, set species=0
  Levels=NULL
  if(n.factors>0){
    for(icol in which.factors){
		  Dat.num[,icol]=as.numeric((Dat[,icol]))
	  }
    Levels=vector("list",n.factors)
	  for(i in 1:n.factors){
	    Levels[[i]]=sort(unique(Dat.num[,which.factors[i]]))
	  }
	  names(Levels)=colnames(Dat[,which.factors])
  }
  #update observer covariate values to reflect new factor values going from 1,2,...
  if(n.obs.cov>0){
	  for(icov in 1:n.obs.cov){
	    if((icov+3)%in%which.factors){
        Obs.cov[,icov]=as.factor(as.numeric(as.factor(Obs.cov[,icov])))
	    }
	  }	
  }
	  
	#add an additional column for "True species" and fill
	True.sp=as.numeric(as.character(Dat[,"Species"]))
  n.missed=sum(is.na(True.sp))
  Which.photo=which(is.na(True.sp)==FALSE)
  if(n.missed>0){
	  True.sp.miss=which(is.na(True.sp)==1)
    if(n.species==1)True.sp[True.sp.miss[1]]=1
	  else True.sp[True.sp.miss[1]]=sample(c(1:(n.species-1)),1)
    if(n.missed>1){
      for(i in 2:length(True.sp.miss)){
        if(n.species==1)True.sp[True.sp.miss[i]]=1
        else True.sp[True.sp.miss[i]]=sample(c(1:n.species),1)
      }
    }
	  #fill group sizes for unphotographed animals 
	  for(imiss in 1:length(True.sp.miss))Dat.num[True.sp.miss[imiss],"Group"]=switch_sample(1,Cov.prior.pdf[True.sp[True.sp.miss[imiss]],1],Cov.prior.parms[True.sp[True.sp.miss[imiss]],,1])
  }

  for(iind in 1:length(Which.photo)){
    True.sp[Which.photo[iind]]=sample(c(1:n.species),1,prob=Psi[,True.sp[Which.photo[iind]],1])
  }
  
	  
	if(Control$update.sp==FALSE)True.sp=True.species #for debugging
	Dat.num=cbind(Dat.num[,1:3],True.sp,Dat.num[,4:ncol(Dat.num)])
	cur.colnames=colnames(Dat)
	cur.colnames[3]="Obs"
	cur.colnames[4]="Species"
  cur.colnames[5:ncol(Dat.num)]=colnames(Dat)[4:ncol(Dat)]
	colnames(Dat.num)=cur.colnames
	
	G.transect=matrix(0,n.species,n.transects)  #number of groups by transect; each row gives results for separate species
	N.transect=G.transect #total abundance by transect
  for(isp in 1:n.species){
    for(itrans in 1:n.transects){
      Cur.ind=which(Dat.num[,"Transect"]==itrans & Dat.num[,"Species"]==isp)
      if(length(Cur.ind)>0){
        G.transect[isp,itrans]=length(Cur.ind)
        N.transect[isp,itrans]=sum(Dat.num[Cur.ind,"Group"])
      }
    }
  }
  Dat=Dat.num


	N.hab.par=rep(0,n.species)
	DM.hab=vector('list',n.species)
	if(1==1){
		if(is.null(Hab.cov)|Hab.formula[[1]]==~1){
			DM.hab[[1]]=as.matrix(rep(1,S),ncol=1)
			colnames(DM.hab[[1]])="Intercept"
		}
		else DM.hab[[1]]=model.matrix(Hab.formula[[1]],data=Hab.cov)
	}
	N.hab.par[1]=ncol(DM.hab[[1]])
	if(n.species>1){
		for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab.pois1
			if(is.null(Hab.cov)|Hab.formula[[i]]==~1){
				DM.hab[[i]]=as.matrix(rep(1,S),ncol=1)
				colnames(DM.hab[[i]])="Intercept"
			}
			else DM.hab[[i]]=model.matrix(Hab.formula[[i]],data=Hab.cov)
			N.hab.par[i]=ncol(DM.hab[[i]])
		}
	}  

	Par=generate_inits_BOSSst(t.steps=t.steps,DM.hab=DM.hab,N.hab.par=N.hab.par,G.transect=G.transect,thin.mean=apply(Thin,1,'mean'),Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,spat.ind=spat.ind,grp.mean=Cov.prior.parms[,1,1])	
  Par$Psi=Psi[,,sample(dim(Psi)[3],1)]
  
  if(is.null(Inits)==FALSE){  #replace random inits with user provided inits for all parameters specified
		I.init=names(Inits)
		for(ipar in 1:length(I.init)){
			eval(parse(text=paste("Par$",names(Inits)[ipar],"=Inits$",names(Inits[ipar]))))
		}
	}
	#start Omega out at a compatible level
  for(isp in 1:n.species){
    if(ncol(Par$hab)<ncol(DM.hab[[isp]]))cat('Error: Abundance intensity model has parameter/formula mismatch')
    #Par$Omega[isp,]=DM.hab[[isp]]%*%Par$hab[isp,1:N.hab.par[isp]]+Par$Eta[isp,]+rnorm(ncol(Par$Omega),0,1/sqrt(Par$tau.eps[isp]))+log(Area.hab)
    Par$Omega[isp,]=DM.hab[[isp]]%*%Par$hab[isp,1:N.hab.par[isp]]+log(Area.hab)
  }
  if(length(Par$tau.eps)!=n.species)cat('Error: length of initial value vector for tau.eps should be equal to # of species')
  
	#get initial individual covariate parameter values
	Par$Cov.par=Cov.prior.parms 
	for(i in 1:n.ind.cov){	
		for(j in 1:n.species){
			if(Cov.prior.fixed[j,i]==1)Par$Cov.par[j,,i]=Cov.prior.parms[j,,i]
			else{
				temp=switch_sample_prior(Cov.prior.pdf[j,i],Cov.prior.parms[j,,i])
				Par$Cov.par[j,1:length(temp),i]=temp
			}
		}
	}
	

	n.hab.cov=ifelse(is.null(Hab.cov)==1 | length(Hab.cov)==1,0,ncol(Hab.cov))
		
	i.Covered=c(1:(S*t.steps))%in%Mapping
	Covered.area=rep(0,S*t.steps)
	for(i in 1:(S*t.steps)){
		if(i.Covered[i]==1){
			Covered.area[i]=sum(Area.trans[which(Mapping==i)])
		}
	}
  
	Meta=list(t.steps=t.steps,n.transects=n.transects,n.species=n.species,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
      DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Mapping=Mapping,Covered.area=Covered.area,
			factor.ind=factor.ind,Levels=Levels,CellTime=CellTime,
			G.transect=G.transect,N.transect=N.transect,grps=grps,n.ind.cov=n.ind.cov,
			Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,
			N.hab.par=N.hab.par,post.loss=post.loss)
  
	cur.iter=Control$adapt
  adapt=1

#	if(adapt==TRUE){
		#cat('\n Beginning adapt phase \n')
		#Out=mcmc_boss_st(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab=DM.hab,Prior.pars=Prior.pars,Meta=Meta,Omega.true=Omega.true,Eta.true=Eta.true)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_boss_st(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$iter,adapt=Control$adapt,Control=Control,DM.hab=DM.hab,Prior.pars=Prior.pars,Meta=Meta,Omega.true=Omega.true,Eta.true=Eta.true,Alpha.true=Alpha.true)
#	}
#	else{
#		cat('\n Beginning MCMC phase \n')
#		Out=mcmc_boss(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,Q=Q,Prior.pars=Prior.pars,Meta=Meta,Omega.true=Omega.true,Eta.true=Eta.true)
	#}
	Out	
}
