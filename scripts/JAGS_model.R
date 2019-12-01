model{
  # Hawkes self-exciting point process model
  for (i in 1:nsites) {
    for (j in 1:nobs) {
      t[i, j] ~ dpois(lambda[i, j])
      
      lambda[i, j] = mu + gamma[i, j]
      
      
      for (n in 1:maxmemory) {
        
          gtmp[i, j, n] = ifelse(history[i, j, n] == 0, 0, alpha * (exp(1) ^ (-beta * history[i, j, n])))
        }
      
      gamma[i, j] = sum(gtmp[i, j,])
      
    }
  }
  
  # Priors
  mu ~ dgamma(0.0001, 0.0001)T(0.0001, 2)

  alpha ~ dgamma(0.0001, 0.0001)T(0.0001, 1)
  
  beta ~ dgamma(0.0001, 0.0001)T(0.0001, 5)
  
  # simulation hawkes
  
    for(i in 1:nsites){
      for(j in 1:nobs){
        sim_t[i,j] ~ dpois(lambda[i,j])
      }
    }
    
  
  
    for (i in 1:nsites) {
      for (j in 1:nobs) {
        t2[i, j] ~ dpois(lambda2)
      }
  
    
    }
  lambda2 ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
  # Regular Poisson model
  
  # simulation of poisson
  
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_t2[i,j] ~ dpois(lambda2)
      
  }
}
} 
