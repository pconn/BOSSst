load("c:/users/paul.conn/git/haulout/test_ribbon.Rdata")
load("c:/users/paul.conn/git/haulout/test_bearded.Rdata")
load("c:/users/paul.conn/git/haulout/test_spotted.Rdata")

set.seed(123456)
#load("c:/users/paul.conn/git/BOSS/BOSS/data/Effort_points.Rdat")

### we need to get a day, time, and location for each surveyed cell;
# Flt_table will get date
#
load('./Area_photographed.rda')
load('BOSS_Flt_dates.Rda')
load('./AlaskaBeringData2012_2013_14Dec2015.Rdat')


startDate <- "2006-04-10 01:00:00"
startDay <- (as.POSIXlt(as.POSIXct(startDate))$yday+1)/365
endDate <- "2006-05-08 23:00:00"
endDay <- (as.POSIXlt(as.POSIXct(endDate))$yday+1)/365
DateRange <- startDay + (0:28)/28*(endDay - startDay)
HourRange <-  as.character(0:23)
fracYear2POSIX(2006,DateRange)


# ------------------------------------------------------------------------------
#                   RIBBON SEALS
# ------------------------------------------------------------------------------



Fit.ribbon <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.ribbon.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.ribbon[i,j] <- glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept")*1 +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.ribbon.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.ribbon.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.ribbon=array(0,dim=c(dim(Fit.ribbon),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.ribbon[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.ribbon=1/(1+exp(-H.ribbon)) #back transform

# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------



Fit.bearded <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.bearded.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.bearded[i,j] <- glmmLDTS.bearded.fit$fixed.effect[
      glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept")*1 +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.bearded.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.bearded.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.bearded=array(0,dim=c(dim(Fit.bearded),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.bearded[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.bearded=1/(1+exp(-H.bearded)) #back transform

# ------------------------------------------------------------------------------
#                   Spotted SEALS
# ------------------------------------------------------------------------------



Fit.spotted <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.spotted.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.spotted[i,j] <- glmmLDTS.spotted.fit$fixed.effect[
      glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept")*1 +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.spotted.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.spotted.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.spotted=array(0,dim=c(dim(Fit.spotted),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.spotted[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.spotted=1/(1+exp(-H.spotted)) #back transform

Haulout.samples=list(spotted=H.spotted,bearded=H.bearded,ribbon=H.ribbon)
#output lookup table
save(Haulout.samples,file="c:/users/paul.conn/git/BOSSst/Haulout_samples.Rdat")
