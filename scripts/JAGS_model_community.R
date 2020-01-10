model{
    
#### Hawkes self-exciting point process model ####
  for (i in 1:nsites) {
    for (j in 1:nobs) {
      t[i, j] ~ dpois(lambdaH[i, j])
      
      lambdaH[i, j] = mu + gamma[i, j]
      
      
      for (n in 1:maxmemory) {
        
          gtmp[i, j, n] = ifelse(history[i, j, n] == 0, 0, alpha * (exp(1) ^ (-beta * history[i, j, n])))
        }
      
      gamma[i, j] = sum(gtmp[i, j,])
      
    }
  }
  
  # Hawkes priors
  alpha ~ dgamma(0.0001, 0.0001)T(0.0001, 1)
  beta ~ dgamma(0.0001, 0.0001)T(0.0001, 5)
  mu ~ dgamma(mu.shape, mu.rate)T(0.0001, 2) # for consistency with poisson
  mu.shape <- pow(x.mu,2)/pow(sd.mu,2)
  mu.rate <- x.mu/pow(sd.mu,2)
  x.mu ~ dunif(0,100)
  sd.mu ~ dunif(0,100)

  # Simulation from Hawkes process
  
    for(i in 1:nsites){
      for(j in 1:nobs){
        sim_tH[i,j] ~ dpois(lambdaH[i,j])
      }
    }
    
  
##### Modified Poisson process model #### 
  
    for (i in 1:nsites) {
      for (j in 1:nobs) {
        tM[i, j] ~ dpois(lambdaM[i])
      }
      lambdaM[i] ~ dgamma(lambdaM.shape, lambdaM.rate)T(0.0001, 2)
    
    }
  
  # Priors for modified Poisson
  
  lambdaM.shape <- pow(x.lambdaM,2)/pow(sd.lambdaM,2)
  lambdaM.rate <- x.lambdaM/pow(sd.lambdaM,2)
  x.lambdaM ~ dunif(0,100)
  sd.lambdaM ~ dunif(0,100)
  
  # Simulation of modified Poisson
  
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_tM[i,j] ~ dpois(lambdaM[i])
      
  }
  }
  
#### Regular Poisson model #### 
  
  for (i in 1:nsites) {
    for (j in 1:nobs) {
      tP[i, j] ~ dpois(lambdaP)
    }
    
  }
  
  # Priors for Poisson model
  # note: the hyper-priors for lambdaP are to be consistent
  
  lambdaP ~ dgamma(lambdaP.shape, lambdaP.rate)T(0.0001, 2)
  lambdaP.shape <- pow(x.lambdaP,2)/pow(sd.lambdaP,2)
  lambdaP.rate <- x.lambdaP/pow(sd.lambdaP,2)
  x.lambdaP ~ dunif(0,100)
  sd.lambdaP ~ dunif(0,100)
  
  # Simulation of regular Poisson
  
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_tP[i,j] ~ dpois(lambdaP)
      
    }
  }
  
  # end jags
} 
