#calculate ad hoc variation inflation factor for misID based on 2013 results
#basically taking ratio of SE with vs. without modeling misID as random effect

load("tweedie_ObsVar_Day_2013.RData")
d1=Out
load("tweedie_ObsVar_Day_2013_noMisIDRE.RData")
d2=Out

mean(((d1$SD$sd[1:4])/(d1$Report$N)) / ((d2$SD$sd[1:4])/(d2$Report$N)))

