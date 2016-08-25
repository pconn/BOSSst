#' Function for MCMC analysis of BOSS surveys using CPIF spatio-temporal model
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"hab": a matrix giving the current iteration's linear model parameters for habitat suitability (Omega); 
#'   	         each row gives parameters for a particular species [note: not all columns need to be used for all species; see Meta$N.hab.par]
#' 		"Omega": a matrix giving the latent habitat suitability values (dims: species, (time-1)*S + sample unit );
#'    "Eta.pois": If Meta$spat.ind==FALSE, matrix latent spatio-temporal random effects (kappa in space-time paper) for Poisson abundance models; one for each cell and for each species
#'	  "tau.eta.pois": If Meta$spat.ind==FALSE, a vector of precisions for spatial-temporal process convolutin for the Poisson component
#'	  "tau.epsilon": Vector of precisions for Omega (exchangeable errors) [one for each species]
#' 		"G": a matrix giving the number of groups of animals in each spatio-temporal sampling unit (species on row, sampling unit on column) 
#' 		"N": a matrix giving the the number of animals in each spatio-temporal sampling unit (same as G)
#' 		"Psi": a matrix holding observation assignment probabilities (true species on rows and observation type on columns)
#' 		"Cov.par": an (n.species X n X n.ind.cov)  array holding parameters of individual covariate distributions.
#' @param Dat   A matrix with the first row giving the spatio-temporal "transect" (one for each spatio-temporal cell sampled),
#'          the second a binary indicator for whether the observation 
#'      had an accompanying photograph, the third gives observation type (if photographed), the fourth holds latent species values,
#'      the fifth to fifth+n.obs.cov give observer/survey condition covariates.  The final columns aftter this hold individual level
#'      covariates (starting with group size)
#' @param Psi An array holding posterior predictions of classification probabilities (from another analysis); dimension (n.species,n.obs.types,# posterior samples)
#' @param cur.iter   Number of iterations to run
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'  "adapt" if true, adapts MH proposals
#'  "n.adapt" If adapt==TRUE, number of iterations to apply adapt algorithm (note: burnin shous be > n.adapt)
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'  "n.files": number of files to spread posterior samples over
#'  "fname": base filename (no extension) to house posterior summaries
#'  "MH.G" A vector giving standard deviation for group abundance MH updates
#'  "species.optim": if TRUE, species updates are optimized for discrete covariates 
#'  "update.sp" If FALSE (default is TRUE), "True.species" must be provided and species updates won't be conducted
#'  "misID" If TRUE, update misID params
#'  "fix.tau.epsilon": if TRUE, fix tau_epsilon to 100
#'  "est.alpha": if TRUE, estimate full scale spatio-temporal random effects; if FALSE, just estimate intercept and slope for each knot 
#'  "post.loss"  If TRUE, observed and predicted detections are compiled for posterior predictive loss 
#'  "GOF"        If TRUE, calculate several posterior predictive checks
#'  "gIVH"       If TRUE, calculate generalized independent variable hull
#'  "start.file" An .RData file.  If present, load starting values from this file (useful for continuing MCMC)
#'  "save.file"  An .RData file.  If present, save the workspace at regular intervals for easier continuation
#'  "save.every" How often to save the workspace; e.g. if 100000, then save every 100000 iterations
#' @param DM.hab	A list design matrix for the fixed effects on abundance intensity (log scale; one for each species)
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'  "beta.tau" precision for Gaussian prior on regression parameters (default 0.1)
#'  "a.eps" alpha parameter for tau_epsilon~Gamma(alpha,beta) (default = 1.0)
#'  "b.eps" beta parameter for tau_epsilon~Gamma(alpha,beta) (default = 0.01)
#'  "a.eta" alpha param for gamma prior precision of spatio-temporal model
#'  "b.eta" beta param for gamma prior precision of spatio-temporal process
#'  "beta0.tau.rw2" prior precision for knot intercepts in rw2 model
#'  "beta1.tau.rw2" prior precision for knot slopes in rw2 model
#' @param Meta	A list object giving a number of other features of the dataset, including:
#'  "t.steps" Number of time steps
#' 	"n.transects"	Number of transects
#'  "n.species"     Number of species
#' 	"S"				Number of strata cells
#'  "K"  Spatial weights matrix for process convolution (spatial only at this stage; dimension S x #knots)
#'  "spat.ind"		Indicator for spatial dependence
#'  "Area.hab"		Vector giving relative area covered by each strata
#'  "Area.trans"	Vector giving fraction of area of relevant strata covered by each transect
#'  "DayHour" A (n.transect X 2) matrix providing row and column entries into the Thin array. Each row corresponds to an entry in Mapping
#'  "Thin" An (n.species X n.days X n.hours X n.iter) array providing n.iter posterior samples of the thinning parameters
#'  "Prop.photo"  Vector giving proportion of area in each transect that is photographed (if NULL, assume 1.0)
#'  "CellTime" A matrix with two columns giving each cell (col 1) and each time (col2) surveyed
#'  "Mapping" 		Vector mapping each transect into a parent strata
#'  "Covered.area"	Vector giving the fraction of each strata covered by transects
#'  "factor.ind"	Indicator vector specifying whether data columns are factors (1) or continuous (0)
#'  "Levels"		a list object, whose elements are comprised of detection model names; each element gives total # of levels in the combined dataset
#'  "G.transect"	vector holding current number of groups of animals present in area covered by each transect		
#'  "N.transect"    vector holding current number of animals present in covered area by each transect
#'  "grps"			indicator for whether observations are for groups rather than individuals
#'  "n.ind.cov" 	Number of individual covariates (distance is not included in this total, but group size is)
#'  "Cov.prior.pdf" character vector giving the probability density function associated with each individual covariate (type ? hierarchical_DS for more info)
#'  "Cov.prior.parms"	An (n.species X n X n.ind.cov) array providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
#'  "Cov.prior.fixed" indicator vector for whether parameters of each covariate distribution should be fixed within estimation routine
#'  "Cov.prior.n" 	(#species X #covariates) Matrix giving number of parameters in each covariate pdf 
#'  "N.hab.par"	    A vector specifying the number of parameters needed for each species' Poisson abundance model
#' @param Omega.true If provided, a vector of true Omega values for debugging
#' @return returns a list with the following objects: 
#' 	"MCMC": An 'mcmc' object (see 'coda' R package) containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
#'  "Obs.N":  Records latent abundance in each transect; dimension is (n.species X # samples X # transects)
#'  "Pred.N": Posterior predictive distribution for abundance in each transect; obtained by sampling a Poisson distribution given current parameter values (with possible zero inflation)
#'  "Obs.det":  if Meta$post.loss=TRUE, a matrix holding observed detection types for posterior predictive loss calculations dim = c(n.transects,n.obs.types) 
#'  "Pred.det": if Meta$post.loss=TRUE, an array holding predicted detection types for posterior predictive loss calculations dim = c(n.mcmc.iter,n.transects,n.obs.types)
#'  "GOF" : results of posterior predictive and sampled predictive GOF tests
#'  In addition, posterior samples of spatio-temporal abundance are saved in list objects 'Post' (strata specific group sizes ("Post$G") and abundance ("Post$N")) via .RData file(s) to save RAM    
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn

mcmc_boss_st<-function(Par,Dat,Psi,cur.iter,adapt,Control,DM.hab,Prior.pars,Meta,Omega.true=NULL,Eta.true=NULL,Alpha.true=NULL){	
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
  n.ST=Meta$S*Meta$t.steps
 	Lam.index=c(1:n.ST)
  thin.pl=sample(c(1:dim(Meta$Thin)[3]),1)
	grp.pl=NULL
	if(Meta$grps==TRUE)grp.pl=which(colnames(Dat)=="Group")
  if(is.null(Control$species.optim))Control$species.optim=TRUE
  
	
	#initialize G.obs (number of groups observed per transect)
  n.obs=nrow(Dat)
	g.tot.obs=colSums(Meta$G.transect) #total number of observations of animals seen at least once
	n.obs.cov=ncol(Dat)-Meta$n.ind.cov-4
 	n.obs.types=dim(Psi)[2]
 	
  n.knots=ncol(Meta$K)
  Q1=linear_adj_RW2(Meta$t.steps) #precision matrix
  Q=Q1
  for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix
  Q.t=t(Q)  
  
  #sort count data by time and cell so that some of the list to vector code will work right
  Order=order(Meta$CellTime[,2][Dat[,"Transect"]],Meta$CellTime[,1][Dat[,"Transect"]])
  Dat=Dat[Order,]
  row.names(Dat)=c(1:nrow(Dat))
  
  Sampled=unique(Meta$Mapping)
  n.unique=length(Sampled) 

  #declare index, function for accessing thinning probabilities for specific transects
  #Thin.pl=Meta$t.steps*(Meta$DayHour[,2]-1)+Meta$DayHour[,1]
  Thin.pl=c(1:Meta$n.transects)
 	get_thin<-function(Thin,Thin.pl,sp,iter){
 	  #Tmp.thin=Thin[sp,,,iter][Thin.pl]
 	  Tmp.thin=Thin[sp,,iter][Thin.pl]
 	  Tmp.thin
 	} 
  Thin.vec=rep(0,n.ST)  #needed in G, N updates
    
  if(is.null(Omega.true)==FALSE){
    for(isp in 1:Meta$n.species){
      Par$Omega[isp,]=as.numeric(Omega.true[isp,,])
    }
  }
  if(is.null(Omega.true)==TRUE)Par$Omega=Par$Omega+rnorm(length(Par$Omega),0,0.2)  #get rid of redundant values for pos-def sample covariance
  
	#initialize pi
  Log.area.adjust=log(Meta$Area.hab)
  One=matrix(1,Meta$t.steps,1)
  Y=Matrix(0,n.ST,Meta$t.steps)
  for(it in 1:Meta$t.steps)Y[((it-1)*Meta$S+1):(it*Meta$S),it]=1
  #Omega.exp=exp(t(t(Par$Omega)+Log.area.adjust))
  Omega.exp=exp(Par$Omega)
  Pi=Pi.mat=Pi.obs=Pi.obs.stnd=vector("list",n.species)
  Pi.obs.stnd.mat=matrix(0,Meta$n.species,length(Sampled))
  for(isp in 1:Meta$n.species){
    D=Diagonal(x=Omega.exp[isp,])
    Pi[[isp]]=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp[isp,]
    Pi.mat[[isp]]=Diagonal(x=as.numeric(Pi[[isp]]))%*%Y
    Pi.obs[[isp]]=Pi.mat[[isp]][Sampled,]   #don't include haulout, transect area here
    #Pi.obs[[isp]]=Pi.mat[[isp]][Meta$Mapping,]*Meta$Area.trans
    Pi.obs.stnd[[isp]]=Pi.obs[[isp]]%*%Diagonal(x=1/colSums(Pi.obs[[isp]]))   
  }
  grp.lam=rep(0,Meta$n.species)
  
  #calculate minimum number alive
  #if so, need to recalculate every time step
  ObsCellTime=Meta$CellTime[Dat[,"Transect"],]
  Ct=matrix(0,Meta$n.species,Meta$t.steps)
  Which_sampled=Which_unsampled=Which_day=vector('list',Meta$t.steps)  #list giving which of the S elements are sampled on each day
  for(it in 1:Meta$t.steps){
    for(isp in 1:Meta$n.species){
      Which.Ct.gt0=which(ObsCellTime[,2]==it & Dat[,"Species"]==isp)
      if(length(Which.Ct.gt0>0))Ct[isp,it]=length(Which.Ct.gt0)
    }
    Which_t = which(Meta$CellTime[,2]==it)
    Which_day[[it]]=Meta$S*(it-1)+c(1:Meta$S)
    if(length(Which_t)>0)Which_sampled[[it]]=Meta$CellTime[Which_t,1]
    Which_unsampled[[it]]=Which_day[[it]]
    if(length(Which_t)>0)Which_unsampled[[it]]=Which_unsampled[[it]][-Which_sampled[[it]]]
  }
  log.Mtp1=log(apply(Ct,1,'max')) #lower bound on log(abundance)
  log.G=log(Par$G.tot)   

 	DM.hab.sampled=DM.hab  #this stuff all needed so DMs remain matrices when intercept only model is used
  for(isp in 1:Meta$n.species)DM.hab.sampled[[isp]]=matrix(DM.hab[[isp]][Sampled,],ncol=ncol(DM.hab[[isp]]))
  DM.hab.unsampled=DM.hab
  for(isp in 1:Meta$n.species)DM.hab.unsampled[[isp]]=matrix(DM.hab[[isp]][-Sampled,],ncol=ncol(DM.hab[[isp]]))  
	Sampled.area.by.strata=rep(0,n.unique)
	for(i in 1:Meta$n.transects)Sampled.area.by.strata[Sampled==Meta$Mapping[i]]=Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]+Meta$Area.trans[i]
  
 	#statistics/matrices needed for regression parameter updates -  for sampled cells only
 	XpXinv=vector('list',Meta$n.species)
 	XpXinvXp=XpXinv
 	for(isp in 1:Meta$n.species){
 	  XpXinv[[isp]]=solve(crossprod(DM.hab[[isp]][Sampled,]))
 	  XpXinvXp[[isp]]=XpXinv[[isp]]%*%t(DM.hab[[isp]][Sampled,])
 	}
 	  
	#initialize MCMC, Acceptance rate matrices
  if(Control$fix.tau.epsilon==TRUE){for(isp in 1:Meta$n.species)Par$tau.eps[isp]=100}
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
  mcmc.length.buffer = mcmc.length/Control$n.files
  counter.buffer=1
  MCMC=list(Psi=array(0,dim=c(dim(Psi)[1:2],mcmc.length)),N.tot=matrix(0,Meta$n.species,mcmc.length),G.tot=matrix(0,Meta$n.species,mcmc.length),Hab=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab))),tau.eta=matrix(0,Meta$n.species,mcmc.length),tau.eps=matrix(0,Meta$n.species,mcmc.length),Cov.par=array(0,dim=c(Meta$n.species,mcmc.length,length(Par$Cov.par[1,,]))))

	Accept=list(G=rep(0,Meta$n.species),Omega=matrix(0,Meta$n.species,Meta$t.steps),Psi=0,thin=0)
	Pred.N=array(0,dim=c(Meta$n.species,mcmc.length.buffer,Meta$n.transects))
  Post=list(N=array(0,dim=c(Meta$n.species,mcmc.length.buffer,n.ST)))
  Post$G=Post$N
	Obs.N=Pred.N
  Old.accept.G=Accept$G
  Old.accept.omega=Accept$Omega


	Which.photo=which(Dat[,"Photo"]==1)
	Which.no.photo=which(Dat[,"Photo"]==0)
	n.photo=length(Which.photo)
	n.misID.updates=round(n.photo*0.5)
	n.no.photo=length(Which.no.photo) 
  Index.misID=c(1:n.misID.updates)
  Index.photo=c(1:n.photo)
    
  #establish framework for individual covariate contributions to species updates (all possible combos of individual covariates)
  if(Meta$n.ind.cov>0){
    if(Meta$n.ind.cov==1)Cov.pointer=matrix(unique(Dat[Which.photo,4+n.obs.cov+1]),length(unique(Dat[Which.photo,4+n.obs.cov+1])),1)
    else{
      Unique.tags=apply(Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)],1,'paste',collapse='')
      Unique.tags=which(!duplicated(Unique.tags))
      Cov.pointer=Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)][Unique.tags,]
    }
    colnames(Cov.pointer)=colnames(Dat)[(4+n.obs.cov+1):ncol(Dat)]   
    Cov.logL=matrix(0,nrow(Cov.pointer),Meta$n.species)  #this will hold log likelihoods for each species
    #fill log likelihoods
    for(isp in 1:Meta$n.species){
      for(icov in 1:Meta$n.ind.cov){
        Cov.logL[,isp]=Cov.logL[,isp]+sapply(Cov.pointer[,icov],'switch_pdf',pdf=Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,,icov],RE=0)       
      }
    }
    Cov.pl=apply(matrix(Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)],ncol=Meta$n.ind.cov),1,'get_place',Pointer=Cov.pointer)      
  }
  if(Control$species.optim==TRUE)Species.cell.probs=matrix(0,n.misID.updates,n.species)
   
  #setup process conv/RW2 space-time model
  n.knots=ncol(Meta$K)
  K=matrix(0,S*Meta$t.steps,n.knots*Meta$t.steps)
  #rearrage K to grab right elements of alpha
  cur.row=0
  for(it in 1:Meta$t.steps){
    for(iknot in 1:n.knots){
      K[(cur.row+1):(cur.row+S),(iknot-1)*Meta$t.steps+it]=Meta$K[,iknot]
    }
    cur.row=cur.row+S
  }
  K=Matrix(K)
  Alpha=matrix(0,n.species,nrow(Q))
  for(isp in 1:Meta$n.species){
    Alpha[isp,]=rrw(Par$tau.eta[isp]*Q) #initial values for space-time random effects
    if(is.null(Alpha.true)==FALSE)Alpha[isp,]=as.numeric(t(Alpha.true[isp,,]))
  }
  if(Control$est.alpha==FALSE)Alpha=Alpha*0
  #Eta=K%*%Alpha
  Eta=matrix(0,Meta$n.species,n.ST) #start Eta at 0
  if(is.null(Eta.true)==FALSE)Eta=Eta.true
  #K.obs=K[Which.obs,]
  #K.obs.t=t(K.obs)
  K.t=t(K)
  cross.K<-crossprod(K,K)
  Diag=diag(nrow(cross.K))*0.001
  X.rw2=matrix(0,n.knots*Meta$t.steps,2*n.knots)
  for(iknot in 1:n.knots){
    X.rw2[((iknot-1)*Meta$t.steps+1):((iknot-1)*Meta$t.steps+Meta$t.steps),iknot]=1
    X.rw2[((iknot-1)*Meta$t.steps+1):((iknot-1)*Meta$t.steps+Meta$t.steps),n.knots+iknot]=c(1:Meta$t.steps)
  }
  X.rw2.t=t(X.rw2)
  A=solve(crossprod(X.rw2),X.rw2.t) #for conditioning by kriging
  A.t=t(A)
  KX.rw2=K%*%X.rw2
  KX.rw2.t=t(KX.rw2)
  KXpX.rw2=crossprod(KX.rw2)
  Sigma.inv.rw2=Matrix(diag(c(rep(Prior.pars$beta0.tau.rw2,n.knots),rep(Prior.pars$beta1.tau.rw2,n.knots))))
  Beta.rw2=matrix(0,Meta$n.species,ncol(X.rw2))

  #summary statistics to quantify number of counts made in each time step
  Control$MH.omega=matrix(0,Meta$n.species,Meta$t.steps)
  Nt.obs=rep(0,Meta$t.steps)
  Time.indices=matrix(0,Meta$t.steps,2) #beginning and ending entries of cells visited for each time step
  cur.pl=1
  for(it in 1:Meta$t.steps){
    Nt.obs[it]=sum(Meta$CellTime[,2]==it)
    Time.indices[it,1]=cur.pl
    Time.indices[it,2]=(cur.pl+Nt.obs[it]-1)
    cur.pl=cur.pl+Nt.obs[it]
    Control$MH.omega[,it]=2.4^2/Nt.obs[it]
  }
  
  #make lists of certain parameters by species, year
  Omega.list=Omega.current=Prop.Sigma=Pi.list=vector("list",Meta$n.species)
  for(isp in 1:Meta$n.species){
    Omega.pred=DM.hab.sampled[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,Sampled]+Log.area.adjust[Sampled]
    Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)*Meta$Area.trans
    for(it in 1:Meta$t.steps){
      if(Nt.obs[it]>0){
        Omega.list[[isp]][[it]]=Par$Omega[isp,Sampled[which(Meta$CellTime[,2]==it)]] #define only for observed
        Pi.list[[isp]][[it]]=Pi.obs.stnd[[isp]][Time.indices[it,1]:Time.indices[it,2],it] 
        Pi.list[[isp]][[it]]=Pi.list[[isp]][[it]]*Cur.thin[Time.indices[it,1]:Time.indices[it,2]]
        Pi.list[[isp]][[it]]=Pi.list[[isp]][[it]]/sum(Pi.list[[isp]][[it]])
        Prop.Sigma[[isp]][[it]]=diag(Nt.obs[it])  #holds covariance matrix for MH proposals
        Omega.current[[isp]][[it]]=matrix(0,Nt.obs[it],100)
        Omega.current[[isp]][[it]][,1]=as.vector(Omega.list[[isp]][[it]])
      }
    }
  }  
  Omega.exp=exp(Par$Omega)
  
  Obs.det=NA
  Pred.det=NA
  Exp.det=NA
  #GOF.counts=NULL
  if(Control$GOF){ #note: some of these are hard wired for BOSS data
    #GOF.counts=array(0,dim=c(mcmc.length,Meta$n.transects,n.obs.types+1))
    if(is.null(Control$samp_pval_n))Control$samp_pval_n=(cur.iter-Control$burnin)/Control$thin
    samp.iter=round(runif(1,1,(cur.iter-Control$burnin)/Control$thin))*Control$thin+Control$burnin #which posterior sample to base sampled predictive p-value on
  }
  if(Control$gIVH){
  }
  if(is.null(Meta$Prop.photo))Meta$Prop.photo=rep(1.0,Meta$n.transects)
  if(Control$post.loss | Control$GOF){ #calculate observed counts of different detection types, initialize prediction arrays
    Obs.det=Exp.det=matrix(0,Meta$n.transects,n.obs.types+1) #col=n.obs.types+1 is for hotspots w/o accompanying photographs
    Pred.det=array(0,dim=c(mcmc.length,Meta$n.transects,n.obs.types+1))
    for(iobs in 1:n.photo)Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]=Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]+Dat[Which.photo[iobs],"Group"]
    Obs.det[,n.obs.types+1]=tabulate(Dat[Which.no.photo,"Transect"],nbins=Meta$n.transects)    
    apply_misID<-function(n,Cur.psi)rmultinom(1,n,prob=Cur.psi)  #define function for sampling observation types by transects
  }
  
	#initialize random effect matrices for individual covariates if required
	if(sum(1-Meta$Cov.prior.fixed)>0)RE.cov=matrix(0,n.obs,Meta$n.ind.cov)
  
  #make copy of Dat, with transect as a factor
  Dat2=as.data.frame(Dat) #establish levels so xtabs works when there are zeros 
  Dat2[,"Transect"]=factor(Dat[,"Transect"],levels=c(1:Meta$n.transects))  
  Dat2[,"Species"]=factor(Dat2[,"Species"],levels=c(1:Meta$n.species)) 
  
  if(is.null(Control$save.every))Control$save.every=100000000
  	
	PROFILE=FALSE  #outputs time associated with updating different groups of parameters
	DEBUG=FALSE
	if(DEBUG){
    #Par$Eta=0*Par$Eta
    #set initial values as in "spatpred" GLMM for comparison
    #Which.pos=which(Meta$G.transect[1,]>0)
    #Offset=Meta$Area.hab[Mapping]*Sampled.area.by.strata
    #Par$Nu[1,Mapping[Which.pos]]=log(Meta$G.transect[1,Which.pos]/Offset[Which.pos])
    #Par$Nu[1,-Mapping[Which.pos]]=min(Par$Nu[1,Mapping[Which.pos]])
    #Par$hab.pois[1,1]=mean(Par$Nu)
    #Par$tau.eta.pois=100
    #Par$tau.nu=100
    #Par$Eta.pois[1,]=rep(0,1299)
    #set.seed(12345)
	}
	if(is.null(Control$start.file)==FALSE){
	  Control_tmp = Control
	  iter_tmp = cur.iter
	  samp_tmp = samp.iter
	  MCMC_tmp=MCMC
	  load(Control$start.file)
	  cur.iter=iter_tmp
	  samp.iter=samp_tmp
	  MCMC=MCMC_tmp
	  MH.omega=Control$MH.omega
	  Control=Control_tmp
	  Control$MH.omega=MH.omega
	}
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
	  #if((iiter%%1000)==1)cat(paste("iteration ",iiter," of ",cur.iter,"\n"))
	  #if(DEBUG){
		#  Par$hab[,1]=c(10,1,7,4)
		#  Par$hab[,2]=c(-10,2,-4,-5)
		#}
		for(isp in 1:Meta$n.species){

		  #update total abundance
		  log.G.prop=log.G[isp]+rnorm(1,0,Control$MH.N[isp])
		  Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)*Meta$Area.trans
		  if(log.G.prop>log.Mtp1[isp]){
		    Cur.p=colSums(Pi.obs[[isp]]*Cur.thin)
		    G.old=round(exp(log.G[isp]))
		    G.prop=round(exp(log.G.prop))
        
        #old.lik=sum(dbinom(Ct[isp,],G.old,Cur.p,log=TRUE))
        #new.lik=sum(dbinom(Ct[isp,],G.prop,Cur.p,log=TRUE))
			  old.lik=Meta$t.steps*lfactorial(G.old)-sum(lfactorial(G.old-Ct[isp,])-(G.old-Ct[isp,])*log(1-Cur.p))
		    new.lik=Meta$t.steps*lfactorial(G.prop)-sum(lfactorial(G.prop-Ct[isp,])-(G.prop-Ct[isp,])*log(1-Cur.p))
		    if(runif(1)<exp(new.lik-old.lik)){
		      Accept$G[isp]=Accept$G[isp]+1
		      log.G[isp]=log.G.prop
          Par$G.tot[isp]=exp(log.G[isp])
		    }
		  }      

		  if(PROFILE==TRUE){
		    cat(paste("G: ", (Sys.time()-st),'\n'))
		    st=Sys.time()
		  } 

      
 			#update omega (sampled cells)
			Omega.pred=DM.hab.sampled[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,Sampled]+Log.area.adjust[Sampled]
			sd=sqrt(1/Par$tau.eps[isp])
			#simulate omega for unobserved times and places
			Par$Omega[isp,-Sampled]=rnorm(n.ST-n.unique,DM.hab.unsampled[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,-Sampled]+Log.area.adjust[-Sampled],sd)
			if(is.null(Omega.true)){
			  for(it in 1:Meta$t.steps){
			    if(Nt.obs[it]>0){
            Prop=rmvnorm(1,Omega.list[[isp]][[it]],Control$MH.omega[isp,it]*Prop.Sigma[[isp]][[it]])
			      Prop.exp=exp(Prop)
			      Pi.prop=Prop.exp/sum(Prop.exp)
            Cur.sum=sum(Pi.prop*Cur.thin[Time.indices[it,1]:Time.indices[it,2]])
			      post.new=sum(dnorm(Prop,Omega.pred[Time.indices[it,1]:Time.indices[it,2]],sd,log=TRUE))+sum(Meta$G.transect[isp,Time.indices[it,1]:Time.indices[it,2]]*(log(Pi.prop)))-sum(Meta$G.transect[isp,Time.indices[it,1]:Time.indices[it,2]])*log(Cur.sum)
			      Cur.sum=sum(Pi.list[[isp]][[it]]*Cur.thin[Time.indices[it,1]:Time.indices[it,2]])
            post.old=sum(dnorm(Omega.list[[isp]][[it]],Omega.pred[Time.indices[it,1]:Time.indices[it,2]],sd,log=TRUE))+sum(Meta$G.transect[isp,Time.indices[it,1]:Time.indices[it,2]]*(log(Pi.list[[isp]][[it]])))-sum(Meta$G.transect[isp,Time.indices[it,1]:Time.indices[it,2]])*log(Cur.sum)
			      if(runif(1)<exp(post.new-post.old)){
			        Omega.list[[isp]][[it]]=Prop
			        Pi.list[[isp]][[it]]=Pi.prop
			        Accept$Omega[isp,it]=Accept$Omega[isp,it]+1
			      }  
			    }
			  }
       
			  Par$Omega[isp,Sampled]=unlist(Omega.list[[isp]])
			  #Omega[Which.obs]=stack_list_vector(Omega.list) #convert back to vector format
			  #compute Pi.obs (need for abundance updates)
			  Omega.exp[isp,]=exp(Par$Omega[isp,])
			  D=Diagonal(x=Omega.exp[isp,])
			  Pi[[isp]]=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp[isp,]
			  Pi.mat[[isp]]=Diagonal(x=as.numeric(Pi[[isp]]))%*%Y
			  Pi.obs[[isp]]=Pi.mat[[isp]][Sampled,] 
			  Pi.obs.stnd[[isp]]=Pi.obs[[isp]]%*%Diagonal(x=1/colSums(Pi.obs[[isp]]))   
			  Pi.obs.stnd.mat[isp,]=rowSums(Pi.obs.stnd[[isp]])  #needed for species updates
			}
			if(is.null(Omega.true)==FALSE){
        Par$Omega[isp,]=as.numeric(Omega.true[isp,,])
        Omega.exp[isp,]=exp(Par$Omega[isp,])
        D=Diagonal(x=Omega.exp[isp,])
        Pi[[isp]]=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp[isp,]
        Pi.mat[[isp]]=Diagonal(x=as.numeric(Pi[[isp]]))%*%Y
        Pi.obs[[isp]]=Pi.mat[[isp]][Sampled,] 
        Pi.obs.stnd[[isp]]=Pi.obs[[isp]]%*%Diagonal(x=1/colSums(Pi.obs[[isp]]))   
        Pi.obs.stnd.mat[isp,]=rowSums(Pi.obs.stnd[[isp]])  #needed for species updates
			}	
      cur.j=iiter%%100
      if(cur.j==0)cur.j=100
      for(it in 1:Meta$t.steps){
        if(Nt.obs[it]>0)Omega.current[[isp]][[it]][,cur.j]=as.vector(Omega.list[[isp]][[it]])
      }
      
		  if(PROFILE==TRUE){
		    cat(paste("Omega: ", (Sys.time()-st),'\n'))
		    st=Sys.time()
		  } 

			#update tau_epsilon	 (precision for exchangeable errors on classification odds)
			if(Control$fix.tau.epsilon==FALSE){
			  Diff=Par$Omega[isp,Sampled]-Omega.pred
			  Par$tau.eps[isp] <- rgamma(1,n.unique/2 + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
			}
      #if(DEBUG)Par$tau.eps[isp]=100
        
			if(PROFILE==TRUE){
			  cat(paste("Tau epsilon: ", (Sys.time()-st),'\n'))
			  st=Sys.time()
			}

      #update fixed effect parameters
			Par$hab[isp,1:Meta$N.hab.par[isp]]=t(rmvnorm(1,XpXinvXp[[isp]]%*%(Par$Omega[isp,Sampled]-Par$Eta[isp,Sampled]-Log.area.adjust[Sampled]),XpXinv[[isp]]/(Par$tau.eps[isp]+Prior.pars$beta.tau)))

      
			if(PROFILE==TRUE){
			  cat(paste("fixed effects: ", (Sys.time()-st),'\n'))
			  st=Sys.time()
			}
			if(iiter==125000){
			  crap=1
			}
			#update kernel weights/REs for space-time model
			if(Meta$spat.ind==FALSE){
			  #first update mean and slope for each knot
			  Dat.minus.Exp=Par$Omega[isp,]-DM.hab[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]-Log.area.adjust-K%*%Alpha[isp,]
			  V.inv.rw2=KXpX.rw2*Par$tau.eps[isp]+Sigma.inv.rw2
			  Beta.rw2[isp,]=t(rmvnorm(1,solve(V.inv.rw2,KX.rw2.t%*%Dat.minus.Exp*Par$tau.eps[isp]),as.matrix(solve(V.inv.rw2))))
        #Beta.rw2[isp,]=0
        #now update alpha and correct for mean=0 and slope=0 constraints
			  Dat.minus.Exp=Par$Omega[isp,]-DM.hab[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]-Log.area.adjust-KX.rw2%*%Beta.rw2[isp,]
			  V.eta.inv <- cross.K*Par$tau.eps[isp]+Par$tau.eta[isp]*Q
			  M.eta <- solve(V.eta.inv,Par$tau.eps[isp]*K.t%*%Dat.minus.Exp)
        if(is.null(Alpha.true)==TRUE & Control$est.alpha==TRUE){
			    Alpha[isp,] <- as.numeric(M.eta + solve(chol(V.eta.inv), rnorm(length(M.eta),0,1)))
			    Alpha[isp,]=as.numeric(Alpha[isp,]-V.eta.inv %*% A.t %*% solve(A %*% V.eta.inv %*% A.t,A%*%Alpha[isp,]))  
        }
        Par$Eta[isp,]=as.numeric(K%*%(X.rw2%*%Beta.rw2[isp,]+Alpha[isp,]))
			  if(is.null(Alpha.true)==FALSE)Par$Eta[isp,]=as.numeric(K%*%(Alpha[isp,]))	
        if(is.null(Eta.true)==FALSE)Par$Eta[isp,]=Eta.true[isp,]
			  #update tau.eta
			  Par$tau.eta[isp] <- rgamma(1, (n.knots*(Meta$t.steps-2))*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha[isp,],Q %*% Alpha[isp,])*0.5) + Prior.pars$b.eta)    
			}
      
			if(PROFILE==TRUE){
			  cat(paste("spatio-temporal effects: ", (Sys.time()-st),'\n'))
			  st=Sys.time()
			}
		}
    

		if(Control$update.sp==TRUE){  
		  Lambda=Pi.obs.stnd.mat*Par$G.tot
		  if(length(Which.no.photo)>0){
		    ########## update species for observations without photos
		    Dat[Which.no.photo,"Species"]=sapply(Dat[Which.no.photo,"Transect"],sample_nophoto_sp,Lam=Lambda,n.sp=Meta$n.species)
		    Dat2[,"Species"]=Dat[,"Species"]
		    Meta$G.transect=xtabs(~Species+Transect,data=Dat2)  #recalculate number of groups per species/transect combo
		    ########## update ind. covariates for observations without photos  #############
		    if(Meta$n.ind.cov>0){
		      for(icov in 1:Meta$n.ind.cov){
		        if(Meta$Cov.prior.pdf[icov]=='poisson_ln' | Meta$Cov.prior.pdf[icov]=='pois1_ln')cur.RE=RE.cov[,icov]
		        else cur.RE=0
		        for(isp in 1:Meta$n.species){
		          I.cur=(Dat[,"Photo"]==0)*(Dat[,"Species"]==isp)
		          if(length(cur.RE>0))cur.RE=cur.RE[I.cur==1]
		          Dat[I.cur==1,4+n.obs.cov+icov]=switch_sample(n=sum(I.cur),pdf=Meta$Cov.prior.pdf[icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
		        }
		      }
		    }	
		  }
		  
		  ##### update true species for observed animals ######
		  if(Control$species.optim==FALSE){
		    Which.sampled=sample(n.photo,n.misID.updates) #replace needs to be false here or there's issues with using Old.sp
		    Old.sp=Dat[Which.photo[Which.sampled],"Species"]
		    Prop.sp=sample(c(1:Meta$n.species),n.misID.updates,replace=TRUE)
		    Cur.obs=Dat[Which.photo[Which.sampled],"Obs"]
		    #covariate contributions
		    MH=sapply(Index.misID,"get_mat_entries",Mat=Cov.logL,Row=Cov.pl[Which.sampled],Col=Prop.sp)-sapply(Index.misID,"get_mat_entries",Mat=Cov.logL,Row=Cov.pl[Which.sampled],Col=Old.sp)
		    #confusion matrix contributions
		    if(Meta$misID)MH=MH+log(sapply(Index.misID,"get_mat_entries",Mat=Par$Psi,Row=Prop.sp,Col=Dat[Which.photo[Which.sampled],"Obs"]))-log(sapply(Index.misID,"get_mat_entries",Mat=Par$Psi,Row=Dat[Which.photo[Which.sampled],"Species"],Col=Dat[Which.photo[Which.sampled],"Obs"]))
		    #Poisson abundance model contributions (new state)
        
		    MH=MH+log(sapply(Index.misID,"get_mat_entries",Mat=Lambda.trans,Row=Prop.sp,Col=Dat[Which.photo[Which.sampled],"Transect"]))
		    Rsamp=runif(n.misID.updates)
		    for(isamp in 1:n.misID.updates){
		      mh=MH[isamp]-log(Lambda.trans[Dat[Which.photo[Which.sampled[isamp]],"Species"],Dat[Which.photo[Which.sampled[isamp]],"Transect"]])
		      if(Rsamp[isamp]<exp(mh)){
		        Meta$G.transect[Old.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]=Meta$G.transect[Old.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]-1
		        Meta$G.transect[Prop.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]=Meta$G.transect[Prop.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]+1
		        Dat[Which.photo[Which.sampled[isamp]],"Species"]=Prop.sp[isamp]
		      }
		    }
		  }
		  
		  if(Control$species.optim==TRUE){
		    if(Meta$n.ind.cov>1)paste("ERROR: Control$species.optim==TRUE currently only implemented for one covariate (group size)")
		    Which.sampled=sample(n.photo,n.misID.updates) #replace needs to be false here or there's issues with using Old.sp
		    #lambda contribution
		    Species.cell.probs=Lambda[,Dat[Which.photo[Which.sampled],"Transect"]]
		    #confusion matrix contributions
		    Species.cell.probs=Species.cell.probs*Par$Psi[,Dat[Which.photo[Which.sampled],"Obs"]]
		    #covariate contributions: GROUP ONLY RIGHT NOW
		    for(isp in 1:Meta$n.species)Species.cell.probs[isp,]=Species.cell.probs[isp,]*exp(sapply(Index.misID,"get_mat_entries",Mat=Cov.logL,Row=Cov.pl[Which.sampled],Col=rep(isp,n.misID.updates)))
		    #Species.cell.probs=Species.cell.probs*t(exp(Cov.logL))[,Dat[Which.photo[Which.sampled],"Group"]]
		    Dat[Which.photo[Which.sampled],"Species"]=apply(Species.cell.probs,2,"sample_species",n.species=n.species)
		  }		
		  if(PROFILE==TRUE){
		    cat(paste("True species: ", (Sys.time()-st),'\n'))
		    st=Sys.time()
		  }
		}
		#if(iiter==1568)cat("\n made it 1 \n")
    
		########## recalculate abundance totals
    Dat2[,"Group"]=Dat[,"Group"]
    Dat2[,"Species"]=factor(Dat[,"Species"],levels=c(1:Meta$n.species))
    Meta$G.transect=xtabs(~Species+Transect,data=Dat2)  #recalculate number of groups per species/transect combo
    Meta$N.transect=xtabs(Group~Species+Transect,data=Dat2)

    
    for(isp in 1:Meta$n.species){
      ########## update group abundance at strata level using most recent updates to covered area abundance
      Par$G[isp,]=0
      Thin.vec[Sampled]=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)
      for(it in 1:Meta$t.steps){
        G.rem=Par$G.tot[isp]-sum(Meta$G.transect[isp,Meta$CellTime[,2]==it])  
        if(G.rem<0){
          cat("Error: unobserved group size < G.total group size \n")
        }
        if(G.rem>0)Par$G[isp,]=Par$G[isp,]+rmultinom(1,G.rem,(1-Meta$Covered.area*Thin.vec)*Pi.mat[[isp]][,it])  #unsurveyed or surveyed but not available
      }
      Par$G[isp,Sampled]=Par$G[isp,Sampled]+Meta$G.transect[isp,]
      grp.lam[isp]=ifelse(Meta$Cov.prior.pdf[isp,1] %in% c("pois1_ln","poisson_ln"),exp(Par$Cov.par[isp,1,1]+(Par$Cov.par[isp,2,1])^2/2),Par$Cov.par[isp,1,1])
      Par$N[isp,]=Par$G[isp,]+rpois(n.ST,grp.lam[isp]*Par$G[isp,]) #add the Par$G since number in group > 0 
    }		

		if(PROFILE==TRUE){
		  cat(paste("Updating N,G ", (Sys.time()-st),'\n'))
		  st=Sys.time()
		}
    
    ##### update thinning parameters in independence chain
    old.post=0
    new.post=0
    #new.ind=round(runif(1,0.5,dim(Thin)[4]+0.5))
    new.ind=round(runif(1,0.5,dim(Thin)[3]+0.5))
    for(isp in 1:Meta$n.species){
		  #G.sampled=rep(0,n.unique) #total number of groups currently in each sampled strata
		  #for(i in 1:Meta$n.transects)G.sampled[Sampled==Meta$Mapping[i]]=G.sampled[Sampled==Meta$Mapping[i]]+Meta$G.transect[isp,i]
      Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)
      New.thin=get_thin(Meta$Thin,Thin.pl,isp,new.ind)
      old.post=old.post+sum(dbinom(Meta$G.transect[isp,],Par$G[isp,Sampled],Meta$Area.trans*Cur.thin,log=TRUE))
      new.post=new.post+sum(dbinom(Meta$G.transect[isp,],Par$G[isp,Sampled],Meta$Area.trans*New.thin,log=TRUE))
      #old.post=old.post+sum(dbinom(G.sampled,Par$G[isp,Sampled],Meta$Area.trans*Cur.thin,log=TRUE))
      #new.post=new.post+sum(dbinom(G.sampled,Par$G[isp,Sampled],Meta$Area.trans*New.thin,log=TRUE))
    }
    if(runif(1)<exp(new.post-old.post)){
      Accept$thin=Accept$thin+1
      thin.pl=new.ind
    }
		if(PROFILE==TRUE){
		  cat(paste("Thinning pars: ", (Sys.time()-st),'\n'))
		  st=Sys.time()
		}
    
    ##### update misID parameters in independence chain    
    if(Control$misID){
      Psi.prop=Psi[,,sample(dim(Psi)[3],1)]
      mh=sum(log(sapply(Index.photo,"get_mat_entries",Mat=Psi.prop,Row=Dat[Which.photo,"Species"],Col=Dat[Which.photo,"Obs"]))-log(sapply(Index.photo,"get_mat_entries",Mat=Par$Psi,Row=Dat[Which.photo,"Species"],Col=Dat[Which.photo,"Obs"])))
      if(runif(1)<exp(mh)){
        Par$Psi=Psi.prop
        Accept$Psi=Accept$Psi+1
      }
      if(PROFILE==TRUE){
        cat(paste("misID pars: ", (Sys.time()-st),'\n'))
        st=Sys.time()
      }
    }
    
		#####update parameters of individual covariate distributions (if fixed=0)
    Cov.logL=Cov.logL*0
		for(isp in 1:Meta$n.species){
      if(sum(Meta$G.transect[isp,])>0){
        for(icov in 1:Meta$n.ind.cov){
          Cur.cov=Dat[Dat[,"Species"]==isp & Dat[,"Photo"]==1,4+n.obs.cov+icov] 
          if(Meta$Cov.prior.fixed[isp,icov]==0){
            if(Meta$Cov.prior.pdf[isp,icov]=="normal")cat("\n Warning: hyper-priors not yet implemented for normal dist. \n")
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson"){
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="pois1"){
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)-length(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson_ln" | Meta$Cov.prior.pdf[isp,icov]=="pois1_ln"){
              cat('RE on covariates not currently implemented for BOSS analysis')
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="multinom"){
              Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov]=rdirichlet(1,Meta$Cov.prior.parms[isp,1:Meta$Cov.prior.n[isp,icov],icov]+tabulate(factor(Cur.cov)))
            }
          }
        }
      }
      #update log likelihood contributions for combinations of different covariate values
      for(icov in 1:Meta$n.ind.cov){
        Cov.logL[,isp]=Cov.logL[,isp]+sapply(Cov.pointer[,icov],'switch_pdf',pdf=Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,,icov],RE=0)       
      }    
		}

    #update Mtp1, etc    
    for(it in 1:Meta$t.steps){
      for(isp in 1:Meta$n.species){
        Which.Ct.gt0=which(ObsCellTime[,2]==it & Dat[,"Species"]==isp)
        if(length(Which.Ct.gt0>0))Ct[isp,it]=length(Which.Ct.gt0)
      }
    }
    log.Mtp1=log(apply(Ct,1,'max')) #lower bound on log(abundance)


    if(PROFILE==TRUE){
			cat(paste("Ind cov pars: ", (Sys.time()-st),'\n'))
			st=Sys.time()
		}
		#if(iiter==1568)cat("\n made it 4 \n")
		
		
		if(Control$adapt & iiter<=Control$n.adapt){
			if(iiter%%10==0){
				for(isp in 1:Meta$n.species){
          if((Accept$G[isp]-Old.accept.G[isp])<3)Control$MH.N[isp]=Control$MH.N[isp]*0.95
          if((Accept$G[isp]-Old.accept.G[isp])>4)Control$MH.N[isp]=Control$MH.N[isp]*1.053
					#for(i in 1:Meta$t.steps){
					#  if((Accept$Omega[isp,i]-Old.accept.omega[isp,i])<2)Control$MH.omega[isp,i]=Control$MH.omega[isp,i]*.95
					#  if((Accept$Omega[isp,i]-Old.accept.omega[isp,i])>3)Control$MH.omega[isp,i]=Control$MH.omega[isp,i]*1.053
				  #}          
				}
				#Old.accept.omega=Accept$Omega
        Old.accept.G=Accept$G
			}
		}
    if(Control$adapt & iiter%%100==0){
      #adapt Omega proposal covariance matrix using Alg 2 of Shaby and Wells (2010)
      R=(Accept$Omega-Old.accept.omega)/100
      for(isp in 1:Meta$n.species){
        for(it in 1:Meta$t.steps){
          if(Nt.obs[it]>0){
            #Omega.mean=rowMeans(Omega.current[[isp]][[it]])
            #Sample.Cov=1/99*crossprod(Omega.current[[isp]][[it]]-Omega.mean)
            Sample.Cov=cov(t(Omega.current[[isp]][[it]]))
            gamma.1=(iiter/100-1)^-0.8
            if(iiter==100)gamma.1=0
            if(iiter>100)Control$MH.omega[isp,it]=sqrt(Control$MH.omega[isp,it]^2*exp(gamma.1*(R[isp,it]-0.234)))
            Prop.Sigma[[isp]][[it]]=Prop.Sigma[[isp]][[it]]+gamma.1*(Sample.Cov-Prop.Sigma[[isp]][[it]])
            if(iiter==200)Prop.Sigma[[isp]][[it]]=Prop.Sigma[[isp]][[it]]+0.1*diag(Nt.obs[it])  #add initial diagonal nugget to make nonsingular
          }
        }
      }
      Old.accept.omega=Accept$Omega      
    }
		
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
      if(is.null(Par$Psi)==FALSE)MCMC$Psi[,,(iiter-Control$burnin)/Control$thin]=Par$Psi
      tmp.pl=((iiter-Control$burnin)/Control$thin)%%mcmc.length.buffer
      for(isp in 1:Meta$n.species){
				MCMC$G.tot[isp,(iiter-Control$burnin)/Control$thin]=Par$G.tot[isp]
				MCMC$N.tot[isp,(iiter-Control$burnin)/Control$thin]=sum(Par$N[isp,])/Meta$t.steps
				MCMC$Hab[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab[isp,]
				MCMC$tau.eta[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta[isp]
				MCMC$tau.eps[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eps[isp]
				MCMC$Cov.par[isp,(iiter-Control$burnin)/Control$thin,]=Par$Cov.par[isp,,]
        Obs.N=NULL
        Pred.N=NULL
				#Obs.N[isp,(iiter-Control$burnin)/Control$thin,]=Meta$N.transect[isp,]
				#Temp.G=Meta$Area.hab[Meta$Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab.pois[[isp]]%*%Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]+Par$Eta.pois[isp,])[Meta$Mapping],sqrt(1/Par$tau.nu[isp])))
        #Pred.N[isp,(iiter-Control$burnin)/Control$thin,]=Temp.G+rpois(Meta$n.transects,grp.lam[isp]*Temp.G)	
        if(tmp.pl!=0){
          Post$G[isp,tmp.pl,]=Par$G[isp,]
				  Post$N[isp,tmp.pl,]=Par$N[isp,]
        }
        else{
          Post$G[isp,mcmc.length.buffer,]=Par$G[isp,]
          Post$N[isp,mcmc.length.buffer,]=Par$N[isp,]
        }				
      }
      if(tmp.pl==0){
				save(Post,file=paste0(Control$fname,counter.buffer,".RData"))
				counter.buffer=counter.buffer+1
      }	
      #omnibus Bayesian p-value GOF tests
      #condition on Omega for these
      #if(Control$GOF){
      #  for(isp in 1:Meta$n.species){
      #    Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)*Meta$Area.trans
      ##    G.sim=rep(0,n.ST)
      #    for(it in 1:Meta$t.steps)G.sim=G.sim+rmultinom(1,Par$G.tot[isp],Pi.mat[[isp]][,it])
      #    G.det=rbinom(Meta$n.transects,G.sim[Sampled],Cur.thin)
      #    Cur.no.photo=rbinom(Meta$n.transects,G.det,1-Meta$Prop.photo)
      #    GOF.counts[(iiter-Control$burnin)/Control$thin,,n.obs.types+1]=GOF.counts[(iiter-Control$burnin)/Control$thin,,n.obs.types+1]+Cur.no.photo
      #    G.det=G.det-Cur.no.photo
      #    GOF.counts[(iiter-Control$burnin)/Control$thin,,1:n.obs.types]=GOF.counts[(iiter-Control$burnin)/Control$thin,,1:n.obs.types]+t(sapply(G.det,"apply_misID",Cur.psi=Par$Psi[isp,]))
      #  }
      #}
      
      #posterior predictions of data given regression & misclassification parameters
      if(Control$post.loss | Control$GOF){
        Exp.det=0*Exp.det
        p_ft=0
        Obs.0 = length(which(Obs.det==0))
        for(isp in 1:Meta$n.species){
          #Omega.pred=DM.hab[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,]+Log.area.adjust
          #sd=sqrt(1/Par$tau.eps[isp])
          #Cur.omega.exp=exp(rnorm(length(Omega.pred),Omega.pred,sd))
          Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)*Meta$Area.trans
          #D=Diagonal(x=Cur.omega.exp)
          #Cur.pi=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Cur.omega.exp
          #Cur.pi.mat=Diagonal(x=as.numeric(Cur.pi))%*%Y
          Cur.pi.mat=Pi.mat[[isp]]
          G.sim=G.exp=rep(0,n.ST)
          for(it in 1:Meta$t.steps){
            G.sim=G.sim+rmultinom(1,Par$G.tot[isp],Cur.pi.mat[,it])
            G.exp=G.exp+Par$G.tot[isp]*Cur.pi.mat[,it]
          }
          G.det=rbinom(Meta$n.transects,G.sim[Sampled],Cur.thin)
          Cur.n=G.exp[Sampled]*Cur.thin
          Cur.psi = cbind(matrix(Par$Psi[isp,],nrow(Exp.det),ncol(Exp.det)-1,byrow=T),0)
          Exp.det=Exp.det+Cur.n*Cur.psi
          Cur.no.photo=rbinom(Meta$n.transects,G.det,1-Meta$Prop.photo)
          Pred.det[(iiter-Control$burnin)/Control$thin,,n.obs.types+1]=Pred.det[(iiter-Control$burnin)/Control$thin,,n.obs.types+1]+Cur.no.photo
          G.det=G.det-Cur.no.photo
          Pred.det[(iiter-Control$burnin)/Control$thin,,1:n.obs.types]=Pred.det[(iiter-Control$burnin)/Control$thin,,1:n.obs.types]+t(sapply(G.det,"apply_misID",Cur.psi=Par$Psi[isp,]))
        }
      }
      if(iiter==samp.iter & Control$GOF){ #calculate sampled predictive p-value
        ft_y=sum((sqrt(Obs.det)-sqrt(Exp.det))^2)
        for(irep in 1:Control$samp_pval_n){
          Pred=0*Pred.det[1,,]
          for(isp in 1:Meta$n.species){
            #Omega.pred=DM.hab[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,]+Log.area.adjust
            #sd=sqrt(1/Par$tau.eps[isp])
            #Cur.omega.exp=exp(rnorm(length(Omega.pred),Omega.pred,sd))
            Cur.thin=get_thin(Meta$Thin,Thin.pl,isp,thin.pl)*Meta$Area.trans
            #D=Diagonal(x=Cur.omega.exp)
            #Cur.pi=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Cur.omega.exp
            #Cur.pi.mat=Diagonal(x=as.numeric(Cur.pi))%*%Y
            Cur.pi.mat=Pi.mat[[isp]]
            G.sim=rep(0,n.ST)
            for(it in 1:Meta$t.steps)G.sim=G.sim+rmultinom(1,Par$G.tot[isp],Cur.pi.mat[,it])
            G.det=rbinom(Meta$n.transects,G.sim[Sampled],Cur.thin)
            Cur.no.photo=rbinom(Meta$n.transects,G.det,1-Meta$Prop.photo)
            Pred[,n.obs.types+1]=Pred[,n.obs.types+1]+Cur.no.photo
            G.det=G.det-Cur.no.photo
            Pred[,1:n.obs.types]=Pred[,1:n.obs.types]+t(sapply(G.det,"apply_misID",Cur.psi=Par$Psi[isp,]))
          }  
          ft=sum((sqrt(Pred)-sqrt(Exp.det))^2)
          p_ft = p_ft + (ft_y<ft)
        }
        p_ft=p_ft/Control$samp_pval_n
      }
		}		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
    if(iiter%%Control$save.every==0)save(list=ls(),file=Control$save.file,envir=environment())
	}
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))

	#convert Out$MCMC into mcmc object for use with coda, S3 methods
	Hab.names=vector("list",Meta$n.species)
	for(isp in 1:Meta$n.species)Hab.names[[isp]]=colnames(DM.hab[[isp]])
	Cov.names=vector("list",Meta$n.ind.cov)
	Cov.par.n=0
	if(Meta$n.ind.cov>0){
		for(icov in 1:Meta$n.ind.cov){
			Par.name=switch(Meta$Cov.prior.pdf[icov],pois1_ln=c("mean.minus.1","sd"),poisson_ln=c("mean","sd"),multinom=paste("prop.cell.",c(1:(Meta$Cov.prior.n[isp,icov]-1)),sep=''),normal="mean",pois1="mean.minus.1",poisson="mean")
			Cov.names[[icov]]=paste(Meta$stacked.names[Meta$dist.pl+icov],".",Par.name,sep='')
			Cov.par.n=Cov.par.n+length(Cov.names[[icov]])
		}
	}

	#compute some post hoc GOF stuff
	if(Control$GOF){
	  n.comp = (cur.iter-Control$burnin)/Control$thin
	  I0=(Obs.det==0)
	  obs0=apply(I0,2,'sum')
	  sum.count=apply(Obs.det,1,'sum')
	  obs.tail = quantile(sum.count,0.95)
	  p_0 =rep(0,ncol(Obs.det))
	  p_tail=0
	  for(icomp in 1:n.comp){
	    I0=(Pred.det[icomp,,]==0)
	    pred0=apply(I0,2,'sum')
	    Sum.pred=apply(Pred.det[icomp,,],1,'sum')
	    pred.tail=quantile(Sum.pred,0.95)
	    p_0 = p_0 + (obs0<pred0) + 0.5*(obs0==pred0)
	    p_tail = p_tail + (obs.tail<pred.tail) + 0.5*(obs.tail==pred.tail)
	  }
	}
	if(is.null(Control$save.file)==FALSE)save(list = ls(),file=Control$save.file,envir=environment())
	if(Meta$n.species>1)Post$Psi=MCMC$Psi
	MCMC=convert.BOSSst.to.mcmc(MCMC=MCMC,N.hab.par=Meta$N.hab.par,Cov.par.n=Cov.par.n,Hab.names=Hab.names,Cov.names=Cov.names,fix.tau.eps=Meta$fix.tau.eps,spat.ind=Meta$spat.ind)
	Out=list(MCMC=MCMC,Accept=Accept,Control=Control,Obs.N=Obs.N,Pred.N=Pred.N,Obs.det=Obs.det,Pred.det=Pred.det,DM.hab=DM.hab,Hab.names=Hab.names)
	if(Control$GOF)Out$GOF=list(p_sampled_ft=p_ft,p_0=p_0/n.comp,p_tail=p_tail/n.comp)
	Out
}

