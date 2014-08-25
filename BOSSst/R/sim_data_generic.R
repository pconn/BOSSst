#' function to simulate spatio-temporal count data over generic simulated landscapes.  This version
#' calls sim_CPIF_misID
#' @param n.species Number of species
#' @param S Number of cells (must be a square number)
#' @param t.steps Number of time steps
#' @param n.transects Number of transects to simulate at each time step
#' @param line.width Proportional of a cell's diameter that is covered by a transect when it surveys a cell
#' @param delta Expected proportion increase/decrease for habitat covariate at each time step
#' @param buffer If TRUE, adds a 10 cell buffer around study area to simulate dynamics (i.e. so pop not closed)
#' @return A list object composed of 
#'          "Data" A list consisting of
#'               "Grid" is a list vector holding covariate data by day,
#'               "Meta" some metadata
#'               "Adj" an adjacency matrix for initial ICAR modeling
#'               "Knot.locations" self explanatory...
#'          "Mapping" a vector providing the space-time locations for each cell sampled
#'          "Effort" which holds simulated effort information 
#'          "Obs" which holds simulated observation data
#'          "N.true" True abundance across space and time
#'          "G.true" Truen # of groups across space and time
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic<-function(n.species,S,t.steps,n.transects,line.width,delta=0,buffer=FALSE){
  source('./BOSSst/R/util_funcs.R')
  source('./BOSSst/R/sim_funcs.R')
  
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
   
  x.len=sqrt(S)
  x.len2=x.len
  if(buffer)x.len2=x.len+10
  S2=x.len2^2
  tau.eta.hab=15
  
  #parameters for matern habiatat covariate (the matern stuff isn't used anymore now that we have a space-time auto process for the covariate)
  kappa=12
  r=.25
  mu=1000
  Dat.matern=rMatClust(kappa,r,mu)  
  X=round((x.len2)*Dat.matern$x+0.5)
  Y=round((x.len2)*Dat.matern$y+0.5)
  X[X<1]=1
  Y[Y<1]=1
  X[X>x.len2]=x.len2
  Y[Y>x.len2]=x.len2
  Grid=matrix(0,x.len2,x.len2)
  for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
  Grid=Grid/max(as.vector(Grid))
  Grid.topo=GridTopology(c(0,0),c(1,1),c(x.len2,x.len2))
  Grid.SpG=SpatialGrid(Grid.topo)
  Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
  Matern.df=as.vector(Grid)
  Matern.df=as.data.frame(cbind(Matern.df,Matern.df^2))
  colnames(Matern.df)=c("matern","matern2")
  Grid.SpPDF=SpatialPolygonsDataFrame(Grid.SpP,data=Matern.df,match.ID=FALSE)
  laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                         "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(Grid.SpPDF)=CRS(laea_180_proj)
  Data=list(Meta="Simulated space-time dataset")  #define generic data structure for space-time estimation
  Data$Grid=vector("list",t.steps)
  Data$Effort=NULL
  for(it in 1:t.steps)Data$Grid[[it]]=Grid.SpPDF

  #Determine which cells should be included, removed when going from S+10 by S+10 grid to S by S grid (if buffer=TRUE)
  if(buffer){
    Which.remove=c(1:(x.len2*5),(S2-x.len2*5+1):S2) #ends
    for(icell in 1:S2){
      if(icell%%x.len2<=5)Which.remove=c(Which.remove,icell)
      if(icell%%x.len2>(x.len2-5))Which.remove=c(Which.remove,icell)
    }
    Which.remove=sort(unique(Which.remove))
    Which.include=rep(1,S2)
    Which.include[Which.remove]=0
    Which.include=which(Which.include==1)
  }
  
  #set up some parameters for covariate s-t process
  tau.cov=0.15
  rho.ar1=0.9
  #sig.ar1=(1-rho.ar1)/tau.cov #to make initial variance approx equal to subsequent variance
  sig.ar1=0.4 
  beta0=0.5
  beta1=log(1+delta)  
  
  #Define knot centers for evolution of habitat covariate using 6 by 6 grid over the larger 30 by 30 grid
  Which.knot=NA
  cur.pl=1
  for(i in 1:x.len2){
    for(j in 1:x.len2){
      if(i%%5==3 & j%%5==3)Which.knot=c(Which.knot,cur.pl)
      cur.pl=cur.pl+1
    }
  }
  Knot.centers=gCentroid(Data$Grid[[1]][Which.knot[-1],],byid=TRUE)
  n.knots=length(Knot.centers)
  #calculate kernel densities at grid cell centroids 
  Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
  Distances=gDistance(Knot.centers,Cell.centroids,byid=TRUE)
  K=matrix(dnorm(Distances,0,5),S2,n.knots)  #knot sd=5 
  K=K/rowSums(K)
  #initial log kernel weights are set via an ICAR(10) prcoess
  Knot.Adj=rect_adj(sqrt(n.knots),sqrt(n.knots))
  Q.knot=-Knot.Adj
  diag(Q.knot)=apply(Knot.Adj,2,'sum')
  Q.knot=Matrix(Q.knot)
  Alpha=matrix(0,n.knots,t.steps)
  Alpha[,1]=rrw(tau.cov*Q.knot)
  #future log kernel weights are set via a random walk process
  for(it in 2:t.steps)Alpha[,it]=rho.ar1*Alpha[,it-1]+rnorm(n.knots,0,sig.ar1)
  #Covariate values are obtained through process convolution
  #Data$Grid[[1]]@data[,1]=exp(beta0+K%*%Alpha[,1])
  #cov.mult=1/max(Data$Grid[[1]]@data[,1])
  cov.mult=1
  for(it in 1:t.steps){
    Data$Grid[[it]]@data[,1]=beta0+beta1*(it-1)+K%*%Alpha[,it]
    Data$Grid[[it]]@data[,1]=Data$Grid[[it]]@data[,1]/max(Data$Grid[[it]]@data[,1])
    Data$Grid[[it]]@data[,1][Data$Grid[[it]]@data[,1]<0]=0
    Data$Grid[[it]]@data[,2]=(Data$Grid[[it]]@data[,1])^2
  }
  #pdf(file="sim_knots.pdf")
  #plot(Data$Grid[[1]],xlab="Easting",ylab="Northing",col='gray')
  #points(Knot.centers,col="blue",pch=20)
  #dev.off()

  #plot_N_map(1,as.matrix(Data$Grid[[2]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")

  
  #Go ahead and switch to S by S for RS2closed, CPIF
  Cur.S=S2
  if(buffer){
    for(it in 1:t.steps)Data$Grid[[it]]=Data$Grid[[it]][Which.include,]
    Cur.S=S
  }
  Data$Adj=rect_adj(sqrt(Cur.S),sqrt(Cur.S)) 
  Data$Knot.locations=Knot.centers
  
  hab.formula=vector("list",n.species)
  for(isp in 1:n.species)hab.formula[[isp]]=~0+matern+matern2
  
  #now simulate data using CPIF model
  K.cpif=K
  if(buffer)K.cpif=K[Which.include,]
  Hab=vector("list",n.species)
  Hab[[1]]=c(10,-10)
  Hab[[2]]=c(1,2)
  Hab[[3]]=c(7,-4)
  Hab[[4]]=c(4,-5)
  Psi=0.8*diag(4)
  Psi[1,2]=0.1
  Psi[1,3]=Psi[1,4]=0.05
  Psi[2,1]=0.1
  Psi[2,3]=Psi[2,4]=0.05
  Psi[3,1]=Psi[3,2]=0.05
  Psi[3,4]=0.1
  Psi[4,1]=Psi[4,2]=0.05
  Psi[4,3]=0.1
  Sim.pars=list(Hab=Hab,grp.sizes=c(1,2,1,3),G=c(20000,5000,15000,10000),tau.eta=c(5,10,15,20),rho.ar1=c(0.95,0.95,0.8,0.5),
                sig.ar1=c(0.1,0.15,0.2,0.1),tau.epsilon=c(2,5,10,20),Psi=Psi)              
  Sim.out=sim_CPIF_misID(S=Cur.S,n.species=n.species,line.width=line.width,n.transects=n.transects,Data=Data,Sim.pars=Sim.pars,hab.formula=hab.formula,Q.knot=Q.knot,K.cpif=K.cpif,Area.adjust=NULL)
  Sim.out$Data=Data      
  return(Sim.out)
}






#' function to simulate spatio-temporal count data over generic simulated landscapes.  This version
#' calls sim_CPIF_misID2
#' @param n.species Number of species
#' @param S Number of cells (must be a square number)
#' @param t.steps Number of time steps
#' @param n.transects Number of transects to simulate at each time step
#' @param line.width Proportional of a cell's diameter that is covered by a transect when it surveys a cell
#' @param delta Expected proportion increase/decrease for habitat covariate at each time step
#' @param buffer If TRUE, adds a 10 cell buffer around study area to simulate dynamics (i.e. so pop not closed)
#' @return A list object composed of 
#'          "Data" A list consisting of
#'               "Grid" is a list vector holding covariate data by day,
#'               "Meta" some metadata
#'               "Adj" an adjacency matrix for initial ICAR modeling
#'               "Knot.locations" self explanatory...
#'          "Mapping" a vector providing the space-time locations for each cell sampled
#'          "Effort" which holds simulated effort information 
#'          "Obs" which holds simulated observation data
#'          "N.true" True abundance across space and time
#'          "G.true" Truen # of groups across space and time
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic2<-function(n.species,S,t.steps,n.transects,line.width,delta=0,buffer=FALSE){
  source('./BOSSst/R/util_funcs.R')
  source('./BOSSst/R/sim_funcs.R')
  
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
  
  x.len=sqrt(S)
  x.len2=x.len
  if(buffer)x.len2=x.len+10
  S2=x.len2^2
  tau.eta.hab=15
  
  #parameters for matern habiatat covariate (the matern stuff isn't used anymore now that we have a space-time auto process for the covariate)
  kappa=12
  r=.25
  mu=1000
  Dat.matern=rMatClust(kappa,r,mu)  
  X=round((x.len2)*Dat.matern$x+0.5)
  Y=round((x.len2)*Dat.matern$y+0.5)
  X[X<1]=1
  Y[Y<1]=1
  X[X>x.len2]=x.len2
  Y[Y>x.len2]=x.len2
  Grid=matrix(0,x.len2,x.len2)
  for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
  Grid=Grid/max(as.vector(Grid))
  Grid.topo=GridTopology(c(0,0),c(1,1),c(x.len2,x.len2))
  Grid.SpG=SpatialGrid(Grid.topo)
  Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
  Matern.df=as.vector(Grid)
  Matern.df=as.data.frame(cbind(Matern.df,Matern.df^2))
  colnames(Matern.df)=c("matern","matern2")
  Grid.SpPDF=SpatialPolygonsDataFrame(Grid.SpP,data=Matern.df,match.ID=FALSE)
  laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                         "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(Grid.SpPDF)=CRS(laea_180_proj)
  Data=list(Meta="Simulated space-time dataset")  #define generic data structure for space-time estimation
  Data$Grid=vector("list",t.steps)
  Data$Effort=NULL
  for(it in 1:t.steps)Data$Grid[[it]]=Grid.SpPDF
  
  #Determine which cells should be included, removed when going from S+10 by S+10 grid to S by S grid (if buffer=TRUE)
  if(buffer){
    Which.remove=c(1:(x.len2*5),(S2-x.len2*5+1):S2) #ends
    for(icell in 1:S2){
      if(icell%%x.len2<=5)Which.remove=c(Which.remove,icell)
      if(icell%%x.len2>(x.len2-5))Which.remove=c(Which.remove,icell)
    }
    Which.remove=sort(unique(Which.remove))
    Which.include=rep(1,S2)
    Which.include[Which.remove]=0
    Which.include=which(Which.include==1)
  }
  
  #set up some parameters for covariate s-t process
  tau.cov=0.15
  rho.ar1=0.9
  #sig.ar1=(1-rho.ar1)/tau.cov #to make initial variance approx equal to subsequent variance
  sig.ar1=0.4 
  beta0=0.5
  beta1=log(1+delta)  
  
  #Define knot centers for evolution of habitat covariate using 6 by 6 grid over the larger 30 by 30 grid
  Which.knot=NA
  cur.pl=1
  for(i in 1:x.len2){
    for(j in 1:x.len2){
      if(i%%5==3 & j%%5==3)Which.knot=c(Which.knot,cur.pl)
      cur.pl=cur.pl+1
    }
  }
  Knot.centers=gCentroid(Data$Grid[[1]][Which.knot[-1],],byid=TRUE)
  n.knots=length(Knot.centers)
  #calculate kernel densities at grid cell centroids 
  Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
  Distances=gDistance(Knot.centers,Cell.centroids,byid=TRUE)
  K=matrix(dnorm(Distances,0,5),S2,n.knots)  #knot sd=5 
  K=K/rowSums(K)
  #initial log kernel weights are set via an ICAR(10) prcoess
  Knot.Adj=rect_adj(sqrt(n.knots),sqrt(n.knots))
  Q.knot=-Knot.Adj
  diag(Q.knot)=apply(Knot.Adj,2,'sum')
  Q.knot=Matrix(Q.knot)
  Alpha=matrix(0,n.knots,t.steps)
  Alpha[,1]=rrw(tau.cov*Q.knot)
  #future log kernel weights are set via a random walk process
  for(it in 2:t.steps)Alpha[,it]=rho.ar1*Alpha[,it-1]+rnorm(n.knots,0,sig.ar1)
  #Covariate values are obtained through process convolution
  #Data$Grid[[1]]@data[,1]=exp(beta0+K%*%Alpha[,1])
  #cov.mult=1/max(Data$Grid[[1]]@data[,1])
  cov.mult=1
  for(it in 1:t.steps){
    Data$Grid[[it]]@data[,1]=beta0+beta1*(it-1)+K%*%Alpha[,it]
    Data$Grid[[it]]@data[,1]=Data$Grid[[it]]@data[,1]/max(Data$Grid[[it]]@data[,1])
    Data$Grid[[it]]@data[,1][Data$Grid[[it]]@data[,1]<0]=0
    Data$Grid[[it]]@data[,2]=(Data$Grid[[it]]@data[,1])^2
  }
  #pdf(file="sim_knots.pdf")
  #plot(Data$Grid[[1]],xlab="Easting",ylab="Northing",col='gray')
  #points(Knot.centers,col="blue",pch=20)
  #dev.off()
  
  #plot_N_map(1,as.matrix(Data$Grid[[2]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  
  
  #Go ahead and switch to S by S for RS2closed, CPIF
  Cur.S=S2
  if(buffer){
    for(it in 1:t.steps)Data$Grid[[it]]=Data$Grid[[it]][Which.include,]
    Cur.S=S
  }
  Data$Adj=rect_adj(sqrt(Cur.S),sqrt(Cur.S)) 
  Data$Knot.locations=Knot.centers
  
  hab.formula=vector("list",n.species)
  for(isp in 1:n.species)hab.formula[[isp]]=~0+matern+matern2
  
  #now simulate data using CPIF model
  K.cpif=K
  if(buffer)K.cpif=K[Which.include,]
  
  Hab=vector("list",n.species)
  Hab[[1]]=c(10,-10)
  Hab[[2]]=c(1,2)
  Hab[[3]]=c(7,-4)
  Hab[[4]]=c(4,-5)
  Psi=0.8*diag(4)
  Psi[1,2]=0.1
  Psi[1,3]=Psi[1,4]=0.05
  Psi[2,1]=0.1
  Psi[2,3]=Psi[2,4]=0.05
  Psi[3,1]=Psi[3,2]=0.05
  Psi[3,4]=0.1
  Psi[4,1]=Psi[4,2]=0.05
  Psi[4,3]=0.1
  Sim.pars=list(Hab=Hab,grp.sizes=c(1,2,1,3),G=c(20000,5000,15000,10000),tau.eta=c(5,10,15,20),rho.ar1=c(0.95,0.95,0.8,0.5),
                sig.ar1=c(0.1,0.15,0.2,0.1),tau.epsilon=c(2,5,10,20),Psi=Psi)              
  Sim.out=sim_CPIF_misID2(S=Cur.S,n.species=n.species,line.width=line.width,n.transects=n.transects,Data=Data,Sim.pars=Sim.pars,hab.formula=hab.formula,Q.knot=Q.knot,K.cpif=K.cpif,Area.adjust=NULL)
  Sim.out$Data=Data      
  return(Sim.out)
}

