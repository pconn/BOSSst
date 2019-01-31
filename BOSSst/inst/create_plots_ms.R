### create plots and results table for manuscript
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')


Table.df = data.frame(Model=c("Day effect","No day effect","Day effect","No day effect"),
                      LogL = rep(0,4),DeltaAIC=rep(0,4),Nbd=rep("",4),Nrn=rep("",4),Nsd=rep("",4))


log_CI <- function(N,SE){
  C = exp(1.96*sqrt(log(1+(SE/N)^2)))
  CI = c(N/C,N*C)
  CI
}

############          
#  2012
############
load("tweedie_ObsVar_2012.RData")  #

t_steps = 29
load('BOSS_data_TMB_2013.Rda')
load('Data2013.RData')
DayHour13=DayHour
n_i13=Data$n_i_real
load('BOSS_data_TMB_2012.Rda')
load('Data2012.RData')

## day effects on detection probability
Species = c("Spotted","Ribbon","Bearded","Ringed")
Models = c("Baseline p","Day shift")
Day.DF = data.frame("Year"=c(rep("2012",Data$n_i_real*4*2),rep("2013",n_i13*4*2)),
                    "Species"=c(rep(Species,each=Data$n_i_real*2),rep(Species,each=n_i13*2)),
                    "Day"=c(rep(DayHour[1:Data$n_i_real,1]+5,8),rep(DayHour13[1:n_i13,1],8)),
                    "Model"= c(rep(rep(c(Models[1],Models[2]),each=Data$n_i_real),4),rep(rep(c(Models[1],Models[2]),each=n_i13),4)),
                    "Det.prob" = 0
)

## day effects on detection probability
n_samp = length(Data$P_i)
for(isp in 1:4){
  Day.DF[which(Day.DF$Species==Species[isp] & Day.DF$Model==Models[2] & Day.DF$Year==2012),"Det.prob"]=Out$Report$Thin_trans[isp,1:Data$n_i_real]
  Day.DF[which(Day.DF$Species==Species[isp] & Day.DF$Model==Models[1] & Day.DF$Year==2012),"Det.prob"]=plogis(Data$thin_mu_logit[(isp-1)*n_samp+c(1:Data$n_i_real)])
}

##  Spatial maps
Report=Out$Report
Coords2=Coords
Coords2[,1] = round(Coords2[,1])
Coords2[,2] = round(Coords2[,2])
n_cells = ncol(Report$Z_s)/t_steps
Nmaps_2012 = vector("list",3)
Day = c(1,15,29)
Species = c("Spotted","Ribbon","Bearded")
Dates=c("April 10","April 25","May 8")
Upper = c(5000,3006,2300)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
axis_incr = 0.5*(max(Coords2[,1])-min(Coords2[,1])+min(Coords2[,2])-max(Coords2[,2]))
y_lim = c(min(Coords2[,2])-axis_incr,max(Coords2[,2])+axis_incr)

#bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))

for(iday in 1:3){
  Nmaps_2012[[iday]]=vector("list",3)  #three species
  Cur.date=Dates[iday]
  for(isp in 1:3){
    Nmaps_2012[[iday]][[isp]]=plot_N_map_xy(N=Report$Z_s[isp,n_cells*(Day[iday]-1)+1:n_cells],XY=Coords2,leg.title="Abundance")+ylim(y_lim)+
      labs(title=paste0(Species[isp],", ",Cur.date))+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=8))+
      scale_fill_gradientn(limits=c(0,Upper[isp]),colours=orangePalette(100))
  }
}
#get legends then replot without
legend1 = cowplot::get_legend(Nmaps_2012[[1]][[1]])
legend2 = cowplot::get_legend(Nmaps_2012[[1]][[2]])
legend3 = cowplot::get_legend(Nmaps_2012[[1]][[3]])
for(iday in 1:3){
  Nmaps_2012[[iday]]=vector("list",3)  #three species
  Cur.date=Dates[iday]
  for(isp in 1:3){
    Nmaps_2012[[iday]][[isp]]=plot_N_map_xy(N=Report$Z_s[isp,n_cells*(Day[iday]-1)+1:n_cells],XY=Coords2,leg.title="Abundance")+ylim(y_lim)+
      labs(title=paste0(Species[isp],", ",Cur.date))+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=8))+
      scale_fill_gradientn(limits=c(0,Upper[isp]),colours=orangePalette(100))+ theme(legend.position="none")
  }
}
pdf(file="Nmaps_2012.pdf")
cowplot::plot_grid(Nmaps_2012[[1]][[1]],
                   Nmaps_2012[[2]][[1]],
                   Nmaps_2012[[3]][[1]],
                   legend1,
                   Nmaps_2012[[1]][[2]],
                   Nmaps_2012[[2]][[2]],
                   Nmaps_2012[[3]][[2]],
                   legend2,
                   Nmaps_2012[[1]][[3]],
                   Nmaps_2012[[2]][[3]],
                   Nmaps_2012[[3]][[3]],
                   legend3,nrow=3,ncol=4,rel_widths=c(1,1,1,.4,1,1,1.4,1,1,1,.4,1,1,1,.4))
dev.off()


##GOF
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
Labels1 = factor(Labels1,levels=Labels1)
Resids.df$Labels=rep(Labels1,each=nrow(Resids))

myplot=ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:20)/20)+facet_wrap(~Labels)+ylab('Frequency')+xlab('Randomized quantile residual')
pdf('GOF2012.pdf')
myplot
dev.off()


# results table

Table.df[,"Nbd"] = Table.df[,"Nrn"] =Table.df[,"Nsd"] = as.character(Table.df[,"Nbd"])
Table.df[1,"LogL"]=Report$jnll
CI = round(log_CI(Report$N[3],Out$SD$sd[3]))
Table.df[1,"Nbd"]=paste0(round(Report$N[3]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[2],Out$SD$sd[2]))
Table.df[1,"Nrn"]=paste0(round(Report$N[2]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[1],Out$SD$sd[1]))
Table.df[1,"Nsd"]=paste0(round(Report$N[1]),' (',CI[1],',',CI[2],')')

load("tweedie_ObsVar_noDay_2012.RData")  #
Report=Out$Report
Table.df[2,"LogL"]=Report$jnll
Table.df[2,"DeltaAIC"]= -2*Table.df[1,"LogL"]+2*Report$jnll-16  #16 is 2.0 * 8 additional params
CI = round(log_CI(Report$N[3],Out$SD$sd[3]))
Table.df[2,"Nbd"]=paste0(round(Report$N[3]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[2],Out$SD$sd[2]))
Table.df[2,"Nrn"]=paste0(round(Report$N[2]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[1],Out$SD$sd[1]))
Table.df[2,"Nsd"]=paste0(round(Report$N[1]),' (',CI[1],',',CI[2],')')

############          
#  2013
############
t_steps=28
load('BOSS_data_TMB_2013.Rda')
load("tweedie_ObsVar_Day_2013.RData")  #
load("Data2013.RData")
Report=Out$Report
n_cells = ncol(Report$Z_s)/t_steps
Table.df[3,"LogL"]=Report$jnll
CI = round(log_CI(Report$N[3],Out$SD$sd[3]))
Table.df[3,"Nbd"]=paste0(round(Report$N[3]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[2],Out$SD$sd[2]))
Table.df[3,"Nrn"]=paste0(round(Report$N[2]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[1],Out$SD$sd[1]))
Table.df[3,"Nsd"]=paste0(round(Report$N[1]),' (',CI[1],',',CI[2],')')

##  Spatial maps
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
Coords2=Coords
Coords2[,1] = round(Coords2[,1])
Coords2[,2] = round(Coords2[,2])
Nmaps_2013 = vector("list",3)
Day = c(1,15,28)
Species = c("Spotted","Ribbon","Bearded")
Dates=c("April 7","April 22","May 4")
#Upper = c(4000,4000,2000)
orangePalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
#bluePalette <- colorRampPalette(rev(brewer.pal(8, "Blues")))

for(iday in 1:3){
  Nmaps_2013[[iday]]=vector("list",3)  #three species
  Cur.date=Dates[iday]
  for(isp in 1:3){
    Nmaps_2013[[iday]][[isp]]=plot_N_map_xy(N=Report$Z_s[isp,n_cells*(Day[iday]-1)+1:n_cells],XY=Coords2,leg.title="Abundance")+ylim(y_lim)+
      labs(title=paste0(Species[isp],", ",Cur.date))+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=8))+
      scale_fill_gradientn(limits=c(0,Upper[isp]),colours=orangePalette(100))
  }
}
#get legends then replot without
legend1 = cowplot::get_legend(Nmaps_2013[[1]][[1]])
legend2 = cowplot::get_legend(Nmaps_2013[[1]][[2]])
legend3 = cowplot::get_legend(Nmaps_2013[[1]][[3]])
for(iday in 1:3){
  Nmaps_2013[[iday]]=vector("list",3)  #three species
  Cur.date=Dates[iday]
  for(isp in 1:3){
    Nmaps_2013[[iday]][[isp]]=plot_N_map_xy(N=Report$Z_s[isp,n_cells*(Day[iday]-1)+1:n_cells],XY=Coords2,leg.title="Abundance")+ylim(y_lim)+
      labs(title=paste0(Species[isp],", ",Cur.date))+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=8))+
      scale_fill_gradientn(limits=c(0,Upper[isp]),colours=orangePalette(100))+ theme(legend.position="none")
  }
}
pdf(file="Nmaps_2013.pdf")
cowplot::plot_grid(Nmaps_2013[[1]][[1]],
                   Nmaps_2013[[2]][[1]],
                   Nmaps_2013[[3]][[1]],
                   legend1,
                   Nmaps_2013[[1]][[2]],
                   Nmaps_2013[[2]][[2]],
                   Nmaps_2013[[3]][[2]],
                   legend2,
                   Nmaps_2013[[1]][[3]],
                   Nmaps_2013[[2]][[3]],
                   Nmaps_2013[[3]][[3]],
                   legend3,nrow=3,ncol=4,rel_widths=c(1,1,1,.4,1,1,1.4,1,1,1,.4,1,1,1,.4))
dev.off()


## Residual plots


## day effects on detection probability
n_samp = length(Data$P_i)
Species = c("Spotted","Ribbon","Bearded","Ringed")
for(isp in 1:4){
  Day.DF[which(Day.DF$Species==Species[isp] & Day.DF$Model==Models[2] & Day.DF$Year==2013),"Det.prob"]=Out$Report$Thin_trans[isp,1:Data$n_i_real]
  Day.DF[which(Day.DF$Species==Species[isp] & Day.DF$Model==Models[1] & Day.DF$Year==2013),"Det.prob"]=plogis(Data$thin_mu_logit[(isp-1)*n_samp+c(1:Data$n_i_real)])
}
Pplot = ggplot(Day.DF)+geom_point(aes(x=Day,y=Det.prob,color=Model,shape=Model),alpha=.1)+facet_grid(Species~Year)+
  ylab('Detection probability')+xlab('Survey Day (from 5 April)')+scale_color_manual(values=c("black","orange"))+
  guides(colour = guide_legend(override.aes = list(alpha=1)))+theme(text = element_text(size=14))
pdf("Day_effect.pdf")
 Pplot
dev.off()

##GOF
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
Labels1 = factor(Labels1,levels=Labels1)
Resids.df$Labels=rep(Labels1,each=nrow(Resids))

myplot=ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:20)/20)+facet_wrap(~Labels)+ylab('Frequency')+xlab('Randomized quantile residual')
pdf('GOF2013.pdf')
 myplot
dev.off()


# results table
load("tweedie_ObsVar_noDay_2013.RData")  #
Report=Out$Report
Table.df[4,"LogL"]=Report$jnll
Table.df[4,"DeltaAIC"]= -2*Table.df[3,"LogL"]+2*Report$jnll-16  #16 is 2.0 * 8 additional params
CI = round(log_CI(Report$N[3],Out$SD$sd[3]))
Table.df[4,"Nbd"]=paste0(round(Report$N[3]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[2],Out$SD$sd[2]))
Table.df[4,"Nrn"]=paste0(round(Report$N[2]),' (',CI[1],',',CI[2],')')
CI = round(log_CI(Report$N[1],Out$SD$sd[1]))
Table.df[4,"Nsd"]=paste0(round(Report$N[1]),' (',CI[1],',',CI[2],')')

library(xtable)
xtable(Table.df)


#############
# covariate effects
n.points = 100
Plot_df = data.frame(Covariate = rep("0",n.points*36))
Plot_df$Cov.val = Plot_df$Response=rep(0,n.points*36)
Plot_df$Year = c(rep(2012,n.points*18),rep(2013,n.points*18))
Plot_df$Species = rep(rep(c("Spotted","Ribbon","Bearded"),each=100),12)
Plot_df$Covariate = rep(rep(c("Dist.land","Dist.break","Depth","Sea.ice","Dist.contour","Dist.edge"),each=n.points*3) ,2)
Plot_df$Year = as.factor(as.character(Plot_df$Year))

#2013 first
load("tweedie_ObsVar_Day_2013.RData")  #
load('BOSS_data_TMB_2013.Rda')
load('Data2013.RData')

Report=Out$Report
Cov_means = colMeans(Data$X_s)
Cov_means[7]=Cov_means[4]^2
Cov_means[3]=sqrt(Cov_means[8])
Cov_means[9]=sqrt(Cov_means[6])
cur.pl = n.points*3*6+1
for(icov in 1:6){
  X_mean = t(matrix(Cov_means,9,n.points))
  min.max = c(min(Data$X_s[,icov]),max(Data$X_s[,icov]))
  Cov.vals = min.max[1]+(min.max[2]-min.max[1])/(n.points-1)*c(0:(n.points-1))
  Plot_df[cur.pl:(cur.pl+3*n.points-1),"Cov.val"]=rep(Cov.vals,3)
  #Cur_X = X_mean
  X_mean[,icov]=Cov.vals
  if(icov==3){
    X_mean[,8]=Cov.vals^2
  }
  if(icov==4){
    X_mean[,7]=Cov.vals^2
  }
  if(icov==6){
    X_mean[,9]=sqrt(Cov.vals)
  }
  for(isp in 1:3){
    Plot_df[cur.pl:(cur.pl+n.points-1),"Response"]=exp(X_mean%*%Report$Beta[isp,])/sum(exp(X_mean%*%Report$Beta[isp,]))
    cur.pl=cur.pl+n.points
  }
}

#now 2012
load("tweedie_ObsVar_Day_2012.RData")
load('BOSS_data_TMB_2012.Rda')
load('Data2012.RData')

Report=Out$Report
Cov_means = colMeans(Data$X_s)
Cov_means[7]=Cov_means[4]^2
Cov_means[3]=sqrt(Cov_means[8])
Cov_means[9]=sqrt(Cov_means[6])
cur.pl = 1
for(icov in 1:6){
  X_mean = t(matrix(Cov_means,9,n.points))
  min.max = c(min(Data$X_s[,icov]),max(Data$X_s[,icov]))
  Cov.vals = min.max[1]+(min.max[2]-min.max[1])/(n.points-1)*c(0:(n.points-1))
  Plot_df[cur.pl:(cur.pl+3*n.points-1),"Cov.val"]=rep(Cov.vals,3)
  #Cur_X = X_mean
  X_mean[,icov]=Cov.vals
  if(icov==3){
    X_mean[,8]=Cov.vals^2
  }
  if(icov==4){
    X_mean[,7]=Cov.vals^2
  }
  if(icov==6){
    X_mean[,9]=sqrt(Cov.vals)
  }
  for(isp in 1:3){
    Plot_df[cur.pl:(cur.pl+n.points-1),"Response"]=exp(X_mean%*%Report$Beta[isp,])/sum(exp(X_mean%*%Report$Beta[isp,]))
    cur.pl=cur.pl+n.points
  }
}




Cov_eff_plot = ggplot(Plot_df)+geom_line(aes(y=Response,x=Cov.val,colour=Species,linetype=Species))+facet_grid(Covariate~Year,scales='free_y')+xlab('Standardized covariate value')+ylab('Relative probability')

pdf('Cov_eff_plots.pdf')
Cov_eff_plot
dev.off()




