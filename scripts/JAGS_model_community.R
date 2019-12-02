model{
  # Hawkes self-exciting point process model
  for (i in 1:nsites) {
    mu[i] ~ dgamma(mu.shape, mu.rate)T(0.0001, 2)
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
  
  mu.shape <- pow(x.mu,2)/pow(sd.mu,2)
  mu.rate <- x.mu/pow(sd.mu,2)
  x.mu ~ dunif(0,100)
  sd.mu ~ dunif(0,100)

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
      lambda2[i] ~ dgamma(lambda2.shape, lambda2.rate)T(0.0001, 2)
    
    }
  
  lambda2.shape <- pow(x.lambda2,2)/pow(sd.lambda2,2)
  lambda2.rate <- x.lambda2/pow(sd.lambda2,2)
  x.lambda2 ~ dunif(0,100)
  sd.lambda2 ~ dunif(0,100)
  # Regular Poisson model
  
  # simulation of poisson
  
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_t2[i,j] ~ dpois(lambda2[i])
      
  }
}
} 
