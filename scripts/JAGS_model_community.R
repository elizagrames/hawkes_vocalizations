model{
  # Hawkes self-exciting point process model
  for (i in 1:nsites) {
    mu[i] ~ dgamma(mu.mu, mu.tau)T(0.0001, 2)
    for (j in 1:nobs) {
      t[i, j] ~ dpois(lambda[i, j])
      
      lambda[i, j] = mu[i] + gamma[i, j]
      
      
      for (n in 1:maxmemory) {
        
          gtmp[i, j, n] = ifelse(history[i, j, n] == 0, 0, alpha * (exp(1) ^ (-beta * history[i, j, n])))
        }
      
      gamma[i, j] = sum(gtmp[i, j,])
      
    }
  }
  
  # Priors

  alpha ~ dgamma(0.0001, 0.0001)T(0.0001, 1)
  
  beta ~ dgamma(0.0001, 0.0001)T(0.0001, 5)
  
  mu.mu ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
  tau.mu ~ dgamma(0.0001, 0.0001)T(0.0001, 10)
  # simulation hawkes
  
    for(i in 1:nsites){
      for(j in 1:nobs){
        sim_t[i,j] ~ dpois(lambda[i,j])
      }
    }
    
  
  
    for (i in 1:nsites) {
      for (j in 1:nobs) {
        t2[i, j] ~ dpois(lambda2[i])
      }
      lambda2[i] ~ dgamma(mu.lambda2, tau.lambda2)T(0.0001, 2)
    
    }
  
  mu.lambda2 ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
  tau.lambda2 ~ dgamma(0.0001, 0.0001)T(0.0001, 10)
  # Regular Poisson model
  
  # simulation of poisson
  
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_t2[i,j] ~ dpois(lambda2[i])
      
  }
}
} 
