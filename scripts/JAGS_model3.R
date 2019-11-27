model{
  # Hawkes self-exciting point process model
  for(k in 1:nchaos){
    
  for (i in 1:nsites) {
    for (j in 1:nobs) {
      t[i, j,k] ~ dpois(lambda[i, j,k])
      
      lambda[i, j,k] = mu[k] + gamma[i, j,k]
      

      
      for (n in 1:maxmemory) {
          if(n > maxmemories[i]){
            gtmp[i, j, n,k] = ifelse(history[i, j, n,k] == 0, 0, alpha[k] * (exp(1) ^ (-beta[k] *
                                                                                         history[i, j, n, k])))
                    } else{
                      gtmp[i,j,n,k] <- 0
                      
        }
      }
      
      gamma[i, j,k] = sum(gtmp[i, j,1:maxmemoy[i],k])
      
    }
  }
  
  # Priors
  mu[k] ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
  
  alpha[k] ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
  
  beta[k] ~ dgamma(0.0001, 0.0001)T(0.0001, 10)
  
  # simulation hawkes
  }
  for(k in 1:nchaos){
    for(i in 1:nsites){
      for(j in 1:nobs){
        sim_t[i,j,k] ~ dpois(lambda[i,j,k])
      }
    }
    
  }
  
  
  for(k in 1:nchaos){
    for (i in 1:nsites) {
      lambda2[i,k] ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
      for (j in 1:nobs) {
        t2[i, j,k] ~ dpois(lambda2[i,k])
      }
    }
    
  }
  # Regular Poisson model
  
  # simulation of poisson
  
for(k in 1:nchaos){
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_t2[i,j,k] ~ dpois(lambda2[i,k])
      
    }
  }
}
} 
