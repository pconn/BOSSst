library('RPostgreSQL')
library('sf')


load("c:/users/paul.conn/git/haulout/test_ribbon.Rdata")
load("c:/users/paul.conn/git/haulout/test_bearded.Rdata")
load("c:/users/paul.conn/git/haulout/test_spotted.Rdata")
load('SampledCells_BOSS2012.RData')  #sampled grid cell centroids and date time in UTC

con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))

grid.sf <- sf::st_read_db(con, 
                          query = "SELECT * FROM base.geo_analysis_grid", 
                          geom_column = "geom")  

grid.wx.2012 <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2012-04-10 00:00:00' AND '2012-05-10 00:00:00'")

# 2013:  4/7 - 5/6
grid.wx.2013 <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2013-04-07 00:00:00' AND '2013-05-08 00:00:00'")

n = nrow(Sampled_cells_ST)
Covs = data.frame(matrix(0,n,6))
SolarT = SolarT = solaR::local2Solar(Sampled_cells_ST$dt,lon=Sampled_cells_ST$lon)
colnames(Covs)=c("day","hour","precip","temp","pressure","wind")

#determine closest grid cell # to each location surveyed
Distances = st_distance(Sampled_cells_ST,st_centroid(grid.sf))
Which_closest = rep(0,n)
for(i in 1:n)Which_closest[i]=which(Distances[i,]==min(Distances[i,]))
Cell_num = grid.sf[Which_closest,]$cell

for(i in 1:n){
  Cur_dat = grid.wx.2012[which(grid.wx.2012$cell==Cell_num[i]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-Sampled_cells_ST[i,]$dt)
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i,3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i,"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}
#transform covariates to scale used in GLMPMs
Covs[,"temp"]=(Covs[,"temp"]-270)/27
Covs[,"pressure"]=(Covs[,"pressure"]-100000)/10000
Covs[,"wind"]=Covs[,"wind"]/10
#now convert day, hour into format used in GLMPMs
Covs$day= (lubridate::yday(SolarT)-120)/10
Covs$day2=Covs$day^2
Covs$day3=Covs$day^3
Covs$hour = factor(lubridate::hour(SolarT),levels=c(0:23))

Hour = lubridate::hour(SolarT)
Sin1 = sin(pi*Hour/12)
Cos1 = cos(pi*Hour/12)
Sin2 = sin(pi*Hour/6)
Cos2 = cos(pi*Hour/6)
Sin3 = sin(pi*Hour/4)
Cos3 = cos(pi*Hour/4)
AS.vec = c("ADULT.F","ADULT.M","SUB","YOY")
L_list = vector("list",4)  
Covs$temp2 = Covs$temp
Covs$Dry = rep(0,n)  #needed for model.matrix
Covs$sin1 = Sin1
Covs$sin2 = Sin2
Covs$sin3 = Sin3
Covs$cos1 = Cos1
Covs$cos2 = Cos2
Covs$cos3 = Cos3
Covs$Northing = -st_coordinates(Sampled_cells_ST)[,2]/3152522  #denominator is mean of spatial points used in analysis; see process_haulout_data4.R in haulout code/paper


# ------------------------------------------------------------------------------
#                   RIBBON SEALS
# ------------------------------------------------------------------------------

FE = test.ribbon$fixed.effects
npar=nrow(test.ribbon$covb)
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.ribbon$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}
  
Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.ribbon$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.36,0.37,0.13,0.14),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_ribbon = M %*% Pred
Var_ho_ribbon = M %*% tcrossprod(Pred_var,M)


# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------

npar=nrow(test.bearded$covb)
FE = test.bearded$fixed.effects
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.bearded$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.bearded$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.27,0.23,0.38,0.12),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_bearded = M %*% Pred
Var_ho_bearded = M %*% tcrossprod(Pred_var,M)
# ------------------------------------------------------------------------------
#                   Spotted SEALS
# ------------------------------------------------------------------------------

npar=nrow(test.spotted$covb)
FE = test.spotted$fixed.effects
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.spotted$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.spotted$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.31,0.24,0.33,0.12),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_spotted = M %*% Pred
Var_ho_spotted = M %*% tcrossprod(Pred_var,M)




#output means, var-cov matrices for HO predictions
HO_out = list(Mu_bd = Pred_ho_bearded,Mu_sd=Pred_ho_spotted,Mu_rn=Pred_ho_ribbon,
              Var_bd = Var_ho_bearded, Var_sd = Var_ho_spotted,Var_rn = Var_ho_ribbon)
save(HO_out,file="Haulout_dists.RData")



#########################
#   2013
#########################

load('SampledCells_BOSS2013.RData')  #sampled grid cell centroids and date time in UTC

n = nrow(Sampled_cells_ST)
Covs = data.frame(matrix(0,n,6))
SolarT = SolarT = solaR::local2Solar(Sampled_cells_ST$dt,lon=Sampled_cells_ST$lon)
colnames(Covs)=c("day","hour","precip","temp","pressure","wind")

#determine closest grid cell # to each location surveyed
Distances = st_distance(Sampled_cells_ST,st_centroid(grid.sf))
Which_closest = rep(0,n)
for(i in 1:n)Which_closest[i]=which(Distances[i,]==min(Distances[i,]))
Cell_num = grid.sf[Which_closest,]$cell

for(i in 1:n){
  Cur_dat = grid.wx.2013[which(grid.wx.2013$cell==Cell_num[i]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-Sampled_cells_ST[i,]$dt)
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i,3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i,"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}
#transform covariates to scale used in GLMPMs
Covs[,"temp"]=(Covs[,"temp"]-270)/27
Covs[,"pressure"]=(Covs[,"pressure"]-100000)/10000
Covs[,"wind"]=Covs[,"wind"]/10
#now convert day, hour into format used in GLMPMs
Covs$day= (lubridate::yday(SolarT)-120)/10
Covs$day2=Covs$day^2
Covs$day3=Covs$day^3
Covs$hour = factor(lubridate::hour(SolarT),levels=c(0:23))

Hour = lubridate::hour(SolarT)
Sin1 = sin(pi*Hour/12)
Cos1 = cos(pi*Hour/12)
Sin2 = sin(pi*Hour/6)
Cos2 = cos(pi*Hour/6)
Sin3 = sin(pi*Hour/4)
Cos3 = cos(pi*Hour/4)
AS.vec = c("ADULT.F","ADULT.M","SUB","YOY")
L_list = vector("list",4)  
Covs$temp2 = Covs$temp
Covs$Dry = rep(0,n)  #needed for model.matrix
Covs$sin1 = Sin1
Covs$sin2 = Sin2
Covs$sin3 = Sin3
Covs$cos1 = Cos1
Covs$cos2 = Cos2
Covs$cos3 = Cos3
Covs$Northing = -st_coordinates(Sampled_cells_ST)[,2]/3152522  #denominator is mean of spatial points used in analysis; see process_haulout_data4.R in haulout code/paper


# ------------------------------------------------------------------------------
#                   RIBBON SEALS
# ------------------------------------------------------------------------------

FE = test.ribbon$fixed.effects
npar=nrow(test.ribbon$covb)
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.ribbon$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.ribbon$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.36,0.37,0.13,0.14),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_ribbon = M %*% Pred
Var_ho_ribbon = M %*% tcrossprod(Pred_var,M)


# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------

npar=nrow(test.bearded$covb)
FE = test.bearded$fixed.effects
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.bearded$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.bearded$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.27,0.23,0.38,0.12),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_bearded = M %*% Pred
Var_ho_bearded = M %*% tcrossprod(Pred_var,M)
# ------------------------------------------------------------------------------
#                   Spotted SEALS
# ------------------------------------------------------------------------------

npar=nrow(test.spotted$covb)
FE = test.spotted$fixed.effects
for(iage in 1:4){
  #design matrix
  Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(test.spotted$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%test.spotted$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.31,0.24,0.33,0.12),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_spotted = M %*% Pred
Var_ho_spotted = M %*% tcrossprod(Pred_var,M)




#output means, var-cov matrices for HO predictions
HO_out = list(Mu_bd = Pred_ho_bearded,Mu_sd=Pred_ho_spotted,Mu_rn=Pred_ho_ribbon,
              Var_bd = Var_ho_bearded, Var_sd = Var_ho_spotted,Var_rn = Var_ho_ribbon)
save(HO_out,file="Haulout_dists2013.RData")




