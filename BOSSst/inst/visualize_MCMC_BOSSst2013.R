# take a look at initial results from 2012 BOSSst analysis (all cameras)

library(STabundance)
setwd('c:/users/paul.conn/git/BOSSst/output')
anal.iter = 720

load("BOSS2013_mcmc_workspace.RData")
MCMC_chain1 = MCMC$N.tot[,1:anal.iter]

load("BOSS2013_mcmc_workspace2.RData")
MCMC_chain2 = MCMC$N.tot[,1:anal.iter]

load("BOSS2012_MCMC_output2.RData")
MCMC_chain2 = MCMC

load("BOSS2012_eco_MCMC_output.RData")
MCMC_chain1_eco = MCMC

load("BOSS2012_eco_MCMC_output2.RData")
MCMC_chain2_eco = MCMC

load("BOSS2012_noST_MCMC_output.RData")
MCMC_chain1_noST = MCMC

load("BOSS2012_noST_MCMC_output2.RData")
MCMC_chain2_noST = MCMC


#construct data frame for plotting abundance by species and chain
n_iter=nrow(MCMC$MCMC)
N_df = matrix(0,4*6*n_iter,5) #rows need # species (4) * # input files, columns are N, species, chain, model,iteration
colnames(N_df)=c("N","species","chain","model","iteration")
N_df=data.frame(N_df)
N_df[,"model"]=rep(c("All.covs","Ecoregion","no.REs"),each=n_iter*4*2)
N_df[,"chain"]=rep(c("1","2","1","2","1","2"),each=4*n_iter)
N_df[,"species"]=rep(rep(c("Spotted","Ribbon","Bearded","Ringed"),each=n_iter),6)
N_df[,"iteration"]=rep(c(1:n_iter),4*6)

N_df[,"N"]=c(
  MCMC_chain1$MCMC[,"Abund.sp1"],
  MCMC_chain1$MCMC[,"Abund.sp2"],
  MCMC_chain1$MCMC[,"Abund.sp3"],
  MCMC_chain1$MCMC[,"Abund.sp4"],
  MCMC_chain2$MCMC[,"Abund.sp1"],
  MCMC_chain2$MCMC[,"Abund.sp2"],
  MCMC_chain2$MCMC[,"Abund.sp3"],
  MCMC_chain2$MCMC[,"Abund.sp4"],
  MCMC_chain1_eco$MCMC[,"Abund.sp1"],
  MCMC_chain1_eco$MCMC[,"Abund.sp2"],
  MCMC_chain1_eco$MCMC[,"Abund.sp3"],
  MCMC_chain1_eco$MCMC[,"Abund.sp4"],
  MCMC_chain2_eco$MCMC[,"Abund.sp1"],
  MCMC_chain2_eco$MCMC[,"Abund.sp2"],
  MCMC_chain2_eco$MCMC[,"Abund.sp3"],
  MCMC_chain2_eco$MCMC[,"Abund.sp4"],
  MCMC_chain1_noST$MCMC[,"Abund.sp1"],
  MCMC_chain1_noST$MCMC[,"Abund.sp2"],
  MCMC_chain1_noST$MCMC[,"Abund.sp3"],
  MCMC_chain1_noST$MCMC[,"Abund.sp4"],
  MCMC_chain2_noST$MCMC[,"Abund.sp1"],
  MCMC_chain2_noST$MCMC[,"Abund.sp2"],
  MCMC_chain2_noST$MCMC[,"Abund.sp3"],
  MCMC_chain2_noST$MCMC[,"Abund.sp4"]
)                

library(ggplot2)
my_plot <- ggplot(N_df,aes(iteration,N,colour=chain)) + geom_line(alpha=0.7)
my_plot <- my_plot + scale_colour_manual(values=c("red","blue")) + facet_grid(model ~ species)
my_plot 

pdf(file="2012_MCMC_chains.pdf")
  my_plot
dev.off()

#histograms of results 
hist(MCMC_chain2$MCMC[,"Abund.sp2"])

my_plot <



load('c:/users/paul.conn/git/BOSSst/AlaskaBeringData2012_2013_14Dec2015.Rdat')
S=1331
t_steps=29
source('c:/users/paul.conn/git/BOSSst/BOSSst/R/util_funcs.R')
load('BOSS2012_eco_mcmc_workspace.RData')
Cur.G=matrix(apply(Post$N[1,,],2,'median'),S,t_steps)
#for(i in 1:t.steps){
plot_N_map(1,Cur.G,Grid=Data$Grid$y2012)




