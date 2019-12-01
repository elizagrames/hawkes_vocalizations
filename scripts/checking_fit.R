files <- paste("./output/", list.files("./output/"), sep="")
i <- 0

i <- i+1
load(files[i])
params <- site_model$BUGSoutput$summary
round(head(params, 20),4)
type <- substr(rownames(params), 1, 5)
only_these <- which(type!="sim_t")
hist(params[only_these,8], main=files[i])
MCMCvis::MCMCtrace(site_model, params=c("alpha", "beta", "mu"), pdf=FALSE)

params[which(params[only_these,8]>1.1),8]
dev.off()


## BAD MODELS

ci_counts <- function(hawk, pois, t){
  nsims <- dim(hawk)[1]
  nsites <- dim(hawk)[2]
  
  tmp.h <- tmp.p <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      tmp.h[s,i] <- length(which(hawk[s,i,]>0))
      tmp.p[s,i] <- length(which(pois[s,i,]>0))
    }
  }
  
  hestimates <- apply(tmp.h, 2, mean)
  hest.ci <- apply(tmp.h, 2, quantile, c(0.025, 0.975))
  pestimates <- apply(tmp.p, 2, mean)
  pest.ci <- apply(tmp.p, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(t)
  tbl <- round(cbind(counts, hestimates, pestimates),2)
  colnames(tbl) <- c("obs", "hawkes", "poisson")
  return(tbl)
}

bad_files <- c()
sites <- sites[-bad_files]
files <- paste("./output/", list.files("./output/"), sep="")
files <- files[-bad_files]
i <- 0

i <- i+1
load(files[i])

hawk <- site_model$BUGSoutput$sims.list$sim_t
pois <- site_model$BUGSoutput$sims.list$sim_t2
t <- t(all_events[, which(colnames(all_events) == sites[i])])

table_1 <- ci_counts(hawk, pois, t)


