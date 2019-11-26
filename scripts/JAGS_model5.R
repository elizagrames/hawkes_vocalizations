model{
  # Hawkes self-exciting point process model
  for(k in 1:chaos){
    
  for (i in 1:nsites) {
    for (j in 1:nobs) {
      t[i, j,k] ~ dpois(lambda[i, j,k])
      
      lambda[i, j,k] = mu[k] + gamma[i, j,k]
      
      for(n in 1:maxmemory){
        tmp[i,j,n,k] = ifelse(n < maxmemory, 0, maxmemory)
        gtmp[i, j, n,k] = ifelse(tmp[i,j,n,k] == 0, 0, alpha[k] * (exp(1) ^ (-beta[k] *
                                                                         history[i, j, n, k])))
      }

      
      gamma[i, j,k] = sum(gtmp[i, j,,k])
      
    }
  }
  
  # Priors
  log(mu[k]) <-  logmu[k]
  logmu[k] ~ dnorm(mumu, mutau)
  
  log(alpha[k]) <-  logalpha[k]
  logalpha[k]    ~ dnorm(alphamu, alphatau)
  
  log(beta[k]) <-  logbeta[k]
  logbeta[k]    ~ dnorm(betamu, betatau)
  
  # simulation hawkes
  }
  
  mumu ~ dnorm(0, 1)
  mutau ~ dgamma(0.0001, 0.0001)
  
  alphamu ~ dnorm(0, 1)
  alphatau ~ dgamma(0.0001, 0.0001)
  
  betamu ~ dnorm(0, 1)
  betatau ~ dgamma(0.0001, 0.0001)
  
  for(k in 1:chaos){
    for(i in 1:nsites){
      for(j in 1:nobs){
        sim_t[i,j,k] ~ dpois(lambda[i,j,k])
      }
    }
    
  }
  
  
  for(k in 1:chaos){
    for (i in 1:nsites) {
      lambda2[i,k] ~ dgamma(0.0001, 0.0001)T(0.0001, 2)
      for (j in 1:nobs) {
        t2[i, j,k] ~ dpois(lambda2[i,k])
      }
    }
    
  }
  # Regular Poisson model
  
  # simulation of poisson
  
for(k in 1:chaos){
  for(i in 1:nsites){
    for(j in 1:nobs){
      sim_t2[i,j,k] ~ dpois(lambda2[i,k])
      
    }
  }
}
} 
