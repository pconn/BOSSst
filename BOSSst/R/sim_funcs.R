#' function to simulate latent abundance with a closed pop ideal free superpopulation model with misID.  This version
#' uses a slightly different spatio-temporal specification than the estimation model (initial alpha is RSR, then AR1 time series after that)
#' @param S Number of cells (must be a square number)
#' @param n.species
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#'    Adj - Adjacency matrix for use in ICAR modeling
#' @param Sim.pars A list holding parameters that describe evolution of the spatio-temporal process.  Included are 
#'         Hab - a list of spatial regression parameters for each species
#'         grp.sizes - vector giving expected group sizes by species
#'         G - vector giving total number of animal groups for each species
#'         tau.eta - vector giving precision for spatial random effects on initial abundance for each species
#'         rho.ar1 - vector giving correlation of AR1 process for knot weights for each species
#'         sig.ar1 - vector giving standard deviation for AR1 process
#'         tau.epsilon - a vector of precision of random error (one for each species)
#'         Psi - a classification matrix with true species on columns and observation types on rows
#' @param hab.formula A list vector of formula object holding the regression model formulation for each species
#' @param Q.knot A structure matrix for reduced rank rsr model for spatio-temporal random effects at the first time step
#' @param K.cpif A matrix holding S by k weights associated with process convolution
#' @param Area.adjust   If provided, a vector allowing for differences in suitable habitat for each cell.  Can be used for different grid cell sizes or different
#'        proportions of suitable habitat (e.g., 1.0 = 100% of habitat is suitable for the focal species)
#' @param n.transects Number of transects to simulate at each time step (set to NULL if effort information is provided instead)
#' @param line.width Proportional diameter covered by a transect passing through a cell
#' @return Output A list including the following objects:
#'       Effort: data.frame holding the following columns:  "Cell", "Time", and "AreaSurveyed"
#'       Obs: data.frame holding following columns: "Cell", "Time", "Obs", "Group"
#'       G.true: an array giving the true number of groups by species, sample unit, time
#'       N.true: " " animals by species, sample unit, time
#'       Mapping: A matrix
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_CPIF_misID<-function(S,n.species,Data,Sim.pars,hab.formula,Q.knot,K.cpif,Area.adjust=NULL,n.transects,line.width){
  require(aster)
  DEBUG=TRUE
  
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  Log.area.adjust=log(Area.adjust)
  t.steps=length(Data$Grid)
  n.knots=ncol(K.cpif)
  #generate true abundance by group
  G.true=array(0,dim=c(n.species,S,t.steps))
  Omega.true=G.true
  Eta.true=array(0,dim=c(n.species,S*t.steps))
  for(isp in 1:n.species){
    if(DEBUG)Sim.pars$tau.epsilon[isp]=100
    X=model.matrix(hab.formula[[isp]],Data$Grid[[1]]@data)  
    Alpha=matrix(0,n.knots,t.steps)
    Cell.probs=matrix(0,S,t.steps)
    Q=Matrix(Sim.pars$tau.eta[isp]*Q.knot)
    Alpha[,1]=rrw(Q)
    Omega.true[isp,,1]=Log.area.adjust+X%*%Sim.pars$Hab[[isp]]+K.cpif%*%Alpha[,1]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp]))
    Cell.probs[,1]=exp(Omega.true[isp,,1])
    G.true[isp,,1]=rmultinom(1,Sim.pars$G[isp],Cell.probs[,1])
    Eta.true[isp,1:S]=K.cpif%*%Alpha[,1]
    #plot_N_map(1,as.matrix(Lambda[,1],ncol=1),Grid=Data$Grid,highlight=c(1,2),cell.width=1,leg.title="Covariate")
    for(it in 2:t.steps){
      X=model.matrix(hab.formula[[isp]],Data$Grid[[it]]@data)  
      Alpha[,it]=Sim.pars$rho.ar1[isp]*Alpha[,it-1]+rnorm(n.knots,0,Sim.pars$sig.ar1[isp]) 
      #if(DEBUG)Alpha[,it]=0
      Eta.true[isp,((it-1)*S+1):(it*S)]=K.cpif%*%Alpha[,it]
      Omega.true[isp,,it]=Log.area.adjust+X%*%Sim.pars$Hab[[isp]]+K.cpif%*%Alpha[,it]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp]))
      Cell.probs[,it]=exp(Omega.true[isp,,it])
      G.true[isp,,it]=rmultinom(1,Sim.pars$G[isp],Cell.probs[,it])
    }
  }
  #now simulate survey track effort info
  Coords.y=unique(coordinates(Data$Grid[[1]])[,1])
  Effort=data.frame(Cell=rep(NA,n.transects*t.steps*sqrt(S)),Time=rep(1:t.steps,each=n.transects*sqrt(S)),AreaSurveyed=rep(line.width,t.steps*n.transects*sqrt(S)))
  cur.pl=1
  for(it in 1:t.steps){
    Cur.y=sample(Coords.y,n.transects)
    Sampled=which(coordinates(Data$Grid[[1]])[,1]%in%Cur.y)
    Effort[cur.pl:(cur.pl+length(Sampled)-1),1]=Sampled
    cur.pl=cur.pl+length(Sampled)
  }
  #simulate observed counts
  
  Obs=matrix(0,10000,4)
  colnames(Obs)=c("Transect","Photo","Obs","Group")
  N.true=array(0,dim=c(n.species,S,t.steps))
  True.sp=0
  counter=0
  for(icell in 1:nrow(Effort)){
    for(isp in 1:n.species){
      if(G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]>0){
        cur.g=rbinom(1,G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]],Effort[icell,"AreaSurveyed"])
        if(cur.g>0){
          for(ig in 1:cur.g){
            True.sp=c(True.sp,isp)
            counter=counter+1
            iphoto=round(runif(1))  # 50% photo'd is hardwired
            grp.size=rktp(1,0,Sim.pars$grp.sizes[isp]) #zero truncated poisson
            obs.sp=which(rmultinom(1,1,Sim.pars$Psi[isp,])==1)
            N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]=N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]+grp.size
            if(iphoto==1)Obs[counter,]=c(icell,1,obs.sp,grp.size)
            else Obs[counter,]=c(icell,0,NA,NA)
          }
          cur.tot=0
        }
        remain=G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]-cur.g
        if(remain>0){
          cur.g.size=rktp(1,0,Sim.pars$grp.sizes[isp],xpred=remain)
          N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]=N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]+cur.g.size
        }
      }
    }
  }
  Obs=Obs[1:counter,]
  
  #now simulate abundance for unsampled cells
  Mapping=matrix(0,n.species,S*t.steps)
  Mapping[isp,]=(Effort$Time-1)*S+Effort$Cell
  Not.sampled=c(1:(S*t.steps))
  Not.sampled=Not.sampled[-Mapping]
  for(isp in 1:n.species){
    for(icell in 1:length(Not.sampled)){ 
      cur.time=ceiling(Not.sampled[icell]/S)
      cur.cell=Not.sampled[icell]-S*(cur.time-1)
      cur.g=G.true[isp,cur.cell,cur.time]
      if(cur.g>0)N.true[isp,cur.cell,cur.time]=rktp(1,0,Sim.pars$grp.sizes[isp],xpred=cur.g)
    }     
  }
  Out=list(Mapping=Mapping,Effort=Effort,N.true=N.true,G.true=G.true,Eta.true=Eta.true,Obs=Obs,True.sp=True.sp[-1],Omega.true=Omega.true)
}

#' function to simulate latent abundance with a closed pop ideal free superpopulation model with misID.
#' This version uses the same spatio-temporal formulation as the estimation model
#' @param S Number of cells (must be a square number)
#' @param n.species
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#'    Adj - Adjacency matrix for use in ICAR modeling
#' @param Sim.pars A list holding parameters that describe evolution of the spatio-temporal process.  Included are 
#'         Hab - a list of spatial regression parameters for each species
#'         grp.sizes - vector giving expected group sizes by species
#'         G - vector giving total number of animal groups for each species
#'         tau.eta - vector giving precision for spatial random effects on initial abundance for each species
#'         rho.ar1 - vector giving correlation of AR1 process for knot weights for each species
#'         sig.ar1 - vector giving standard deviation for AR1 process
#'         tau.epsilon - a vector of precision of random error (one for each species)
#'         Psi - a classification matrix with true species on columns and observation types on rows
#' @param hab.formula A list vector of formula object holding the regression model formulation for each species
#' @param Q.knot A structure matrix for reduced rank rsr model for spatio-temporal random effects at the first time step
#' @param K.cpif A matrix holding S by k weights associated with process convolution
#' @param Area.adjust   If provided, a vector allowing for differences in suitable habitat for each cell.  Can be used for different grid cell sizes or different
#'        proportions of suitable habitat (e.g., 1.0 = 100% of habitat is suitable for the focal species)
#' @param n.transects Number of transects to simulate at each time step (set to NULL if effort information is provided instead)
#' @param line.width Proportional diameter covered by a transect passing through a cell
#' @return Output A list including the following objects:
#'       Effort: data.frame holding the following columns:  "Cell", "Time", and "AreaSurveyed"
#'       Obs: data.frame holding following columns: "Cell", "Time", "Obs", "Group"
#'       G.true: an array giving the true number of groups by species, sample unit, time
#'       N.true: " " animals by species, sample unit, time
#'       Mapping: A matrix
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_CPIF_misID2<-function(S,n.species,Data,Sim.pars,hab.formula,Q.knot,K.cpif,Area.adjust=NULL,n.transects,line.width){
  require(aster)
  DEBUG=TRUE
  
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  Log.area.adjust=log(Area.adjust)
  t.steps=length(Data$Grid)
  n.knots=ncol(K.cpif)
  #generate true abundance by group
  G.true=array(0,dim=c(n.species,S,t.steps))
  Omega.true=G.true
  Alpha.true=array(0,dim=c(n.species,n.knots,t.steps))
  Eta.true=array(0,dim=c(n.species,S*t.steps))
  Q1=linear_adj_RW2(t.steps) #precision matrix
  Q=Q1
  for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix
  
  for(isp in 1:n.species){
    if(DEBUG)Sim.pars$tau.epsilon[isp]=100
    X=model.matrix(hab.formula[[isp]],Data$Grid[[1]]@data)  
    Alpha=t(matrix(rrw(Q*Sim.pars$tau.eta[isp]),t.steps,n.knots))
    Cell.probs=matrix(0,S,t.steps)
    #Q=Matrix(Sim.pars$tau.eta[isp]*Q)
    #Alpha[,1]=rrw(Q)
    Alpha.true[isp,,]=Alpha
    Omega.true[isp,,1]=Log.area.adjust+X%*%Sim.pars$Hab[[isp]]+K.cpif%*%Alpha[,1]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp]))
    Cell.probs[,1]=exp(Omega.true[isp,,1])
    G.true[isp,,1]=rmultinom(1,Sim.pars$G[isp],Cell.probs[,1])
    Eta.true[isp,1:S]=K.cpif%*%Alpha[,1]
    #plot_N_map(1,as.matrix(Lambda[,1],ncol=1),Grid=Data$Grid,highlight=c(1,2),cell.width=1,leg.title="Covariate")
    for(it in 2:t.steps){
      X=model.matrix(hab.formula[[isp]],Data$Grid[[it]]@data)  
      #Alpha[,it]=Sim.pars$rho.ar1[isp]*Alpha[,it-1]+rnorm(n.knots,0,Sim.pars$sig.ar1[isp]) 
      #if(DEBUG)Alpha[,it]=0
      Eta.true[isp,((it-1)*S+1):(it*S)]=K.cpif%*%Alpha[,it]
      Omega.true[isp,,it]=Log.area.adjust+X%*%Sim.pars$Hab[[isp]]+K.cpif%*%Alpha[,it]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp]))
      Cell.probs[,it]=exp(Omega.true[isp,,it])
      G.true[isp,,it]=rmultinom(1,Sim.pars$G[isp],Cell.probs[,it])
    }
  }
  #now simulate survey track effort info
  Coords.y=unique(coordinates(Data$Grid[[1]])[,1])
  Effort=data.frame(Cell=rep(NA,n.transects*t.steps*sqrt(S)),Time=rep(1:t.steps,each=n.transects*sqrt(S)),AreaSurveyed=rep(line.width,t.steps*n.transects*sqrt(S)))
  cur.pl=1
  for(it in 1:t.steps){
    Cur.y=sample(Coords.y,n.transects)
    Sampled=which(coordinates(Data$Grid[[1]])[,1]%in%Cur.y)
    Effort[cur.pl:(cur.pl+length(Sampled)-1),1]=Sampled
    cur.pl=cur.pl+length(Sampled)
  }
  #simulate observed counts
  
  Obs=matrix(0,500000,4)
  colnames(Obs)=c("Transect","Photo","Obs","Group")
  N.true=array(0,dim=c(n.species,S,t.steps))
  True.sp=0
  counter=0
  for(icell in 1:nrow(Effort)){
    for(isp in 1:n.species){
      if(G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]>0){
        cur.g=rbinom(1,G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]],Effort[icell,"AreaSurveyed"])
        if(cur.g>0){
          for(ig in 1:cur.g){
            True.sp=c(True.sp,isp)
            counter=counter+1
            iphoto=round(runif(1))  # 50% photo'd is hardwired
            grp.size=rktp(1,0,Sim.pars$grp.sizes[isp]) #zero truncated poisson
            obs.sp=which(rmultinom(1,1,Sim.pars$Psi[isp,])==1)
            N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]=N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]+grp.size
            if(iphoto==1)Obs[counter,]=c(icell,1,obs.sp,grp.size)
            else Obs[counter,]=c(icell,0,NA,NA)
          }
          cur.tot=0
        }
        remain=G.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]-cur.g
        if(remain>0){
          cur.g.size=rktp(1,0,Sim.pars$grp.sizes[isp],xpred=remain)
          N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]=N.true[isp,Effort[icell,"Cell"],Effort[icell,"Time"]]+cur.g.size
        }
      }
    }
  }
  Obs=Obs[1:counter,]
  
  #now simulate abundance for unsampled cells
  Mapping=matrix(0,n.species,S*t.steps)
  Mapping[isp,]=(Effort$Time-1)*S+Effort$Cell
  Not.sampled=c(1:(S*t.steps))
  Not.sampled=Not.sampled[-Mapping]
  for(isp in 1:n.species){
    for(icell in 1:length(Not.sampled)){ 
      cur.time=ceiling(Not.sampled[icell]/S)
      cur.cell=Not.sampled[icell]-S*(cur.time-1)
      cur.g=G.true[isp,cur.cell,cur.time]
      if(cur.g>0)N.true[isp,cur.cell,cur.time]=rktp(1,0,Sim.pars$grp.sizes[isp],xpred=cur.g)
    }     
  }
  Out=list(Mapping=Mapping,Effort=Effort,N.true=N.true,G.true=G.true,Eta.true=Eta.true,Alpha.true=Alpha.true,Obs=Obs,True.sp=True.sp[-1],Omega.true=Omega.true)
}

  