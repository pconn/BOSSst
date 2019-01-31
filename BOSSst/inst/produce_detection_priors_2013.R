#formulate prior for detection, including availability, disturbance, and detection probability of
#thermal detection algorithms
load("SampledCells_BOSS2013.RData")
load("Haulout_dists2013.RData")  #from 'get_haulout_distributions_weather.R'

# ribbon  = no disturbance effect
n_surveyed = length(HO_out$Mu_rn)
Mu_ribbon = as.numeric(HO_out$Mu_rn) * Sampled_cells_ST$Ep
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_rn[i]*HO_out$Mu_rn[j]
    Y[i,j]= Sampled_cells_ST$Ep[i]*Sampled_cells_ST$Ep[j]
  }
}
Var_ribbon = Y*HO_out$Var_rn +
             X*diag(Sampled_cells_ST$Varp)-
             HO_out$Var_rn*diag(Sampled_cells_ST$Varp)

# bearded = no disturbance effect
Mu_bearded = as.numeric(HO_out$Mu_bd) * Sampled_cells_ST$Ep
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_bd[i]*HO_out$Mu_bd[j]
    Y[i,j]= Sampled_cells_ST$Ep[i]*Sampled_cells_ST$Ep[j]
  }
}
Var_bearded = Y*HO_out$Var_bd +
  X*diag(Sampled_cells_ST$Varp)-
  HO_out$Var_bd*diag(Sampled_cells_ST$Varp)

# spotted = with disturbance effect
mu_d = 23/24
var_d = mu_d*(1-mu_d)/24
Mu_dp  = Sampled_cells_ST$Ep * mu_d
Var_dp = mu_d^2*Sampled_cells_ST$Varp+Sampled_cells_ST$Ep*var_d-var_d*Sampled_cells_ST$Varp
Mu_spotted = as.numeric(HO_out$Mu_sd) * Mu_dp
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_sd[i]*HO_out$Mu_sd[j]
    Y[i,j]= Mu_dp[i]*Mu_dp[j]
  }
}
Var_spotted = Y*HO_out$Var_sd +
  X*diag(Var_dp)-
  HO_out$Var_sd*diag(Var_dp)

# ringed = with disturbance effect
Mu_ringed = 0.65 * Sampled_cells_ST$Ep * 23/29

Det_priors = list(p_bd = Mu_bearded, Var_bd=Var_bearded,p_sd = Mu_spotted,Var_sd=Var_spotted,
                  p_rn = Mu_ribbon, Var_rn = Var_ribbon, p_rd = Mu_ringed)

save(Det_priors,file="detection_priors2013.RData")


