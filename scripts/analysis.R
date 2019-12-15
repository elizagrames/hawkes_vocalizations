#### -1: Functions ####
ci_counts <- function(hawk, pois, t){
  nsims <- dim(hawk)[1]
  nsites <- dim(hawk)[2]
  
  tmp.h <- tmp.p <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      tmp.h[s,i] <- length(which(hawk[s,i,]>0)) - sum(t[i,])
      tmp.p[s,i] <- length(which(pois[s,i,]>0)) - sum(t[i,])
    }
  }
  
  hestimates <- apply(tmp.h, 2, mean)
  hest.ci <- apply(tmp.h, 2, quantile, c(0.025, 0.975))
  pestimates <- apply(tmp.p, 2, mean)
  pest.ci <- apply(tmp.p, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(t)
  tbl <- round(cbind(counts, hest.ci[1,], hestimates, hest.ci[2,], pest.ci[1,], pestimates, pest.ci[2,]),2)
  colnames(tbl) <- c("obs", "h.low", "hawkes", "h.up", "p.low", "poisson", "p.up")
  return(tbl)
}

ci_counts2 <- function(hawk, pois, t){
  nsims <- dim(hawk)[1]
  nsites <- dim(hawk)[2]
  
  tmp.h <- tmp.p <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      tmp.h[s,i] <- (length(which(hawk[s,i,]>0)) - sum(t[i,]))/sum(t[i,])
      tmp.p[s,i] <- (length(which(pois[s,i,]>0)) - sum(t[i,]))/sum(t[i,])
    }
  }
  
  hestimates <- apply(tmp.h, 2, mean)
  hest.ci <- apply(tmp.h, 2, quantile, c(0.025, 0.975))
  pestimates <- apply(tmp.p, 2, mean)
  pest.ci <- apply(tmp.p, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(t)
  tbl <- round(cbind(counts, hest.ci[1,], hestimates, hest.ci[2,], pest.ci[1,], pestimates, pest.ci[2,]),2)
  colnames(tbl) <- c("obs", "h.low", "hawkes", "h.up", "p.low", "poisson", "p.up")
  return(tbl)
}

mcc <- function(tp, tn, fp, fn){
  ((tp*tn) - (fp*fn))/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
}

calculate_props <- function(hawk, pois, t){
  nsims <- dim(hawk)[1]
  nsites <- dim(hawk)[2]
  nobs <- dim(hawk)[3]
  
  tmp.h <- tmp.p <- diffph <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      obs.e <- which(t[i,]>0)
      window <- obs.e
      for(m in 1:length(obs.e)){
        window <- append(window, obs.e[m]+1)
        window <- append(window, obs.e[m]-1)
      }
      hawk.pos <- which(hawk[s,i,]>0)
      hawk.neg <- which(hawk[s,i,]==0)
      
      pois.pos <- which(pois[s,i,]>0)
      pois.neg <- which(hawk[s,i,]==0)
      
      obs.pos <- obs.e
      obs.neg <- which(t[i,]==0)
      
      tp.h <- sum(hawk.pos %in% window) 
      tn.h <- sum(hawk.neg %in% obs.neg)
      fp.h <- sum((hawk.pos %in% window)!=TRUE)
      fn.h <- sum(hawk.neg %in% obs.pos)
      
      tp.p <- sum(pois.pos %in% window) 
      tn.p <- sum(pois.neg %in% obs.neg)
      fp.p <- sum((pois.pos %in% window)!=TRUE)
      fn.p <- sum(pois.neg %in% obs.pos)
      
      tmp.h[s,i] <- mcc(tp.h, tn.h, fp.h, fn.h)
      tmp.p[s,i] <- mcc(tp.p, tn.p, fp.p, fn.p)
      diffph[s,i] <- tmp.h[s,i] - tmp.p[s,i]
      
    }
  }
  tmp.h[is.na(tmp.h)] <- NA
  tmp.p[is.na(tmp.p)] <- NA
  
  hestimates <- apply(tmp.h, 2, mean, na.rm=TRUE)
  
  hest.ci <- apply(tmp.h, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  pestimates <- apply(tmp.p, 2, mean, na.rm=TRUE)
  pest.ci <- apply(tmp.p, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  destimates <- apply(diffph, 2, mean, na.rm=TRUE)
  dest.ci <- apply(diffph, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  
  tbl <- cbind(hest.ci[1,], hestimates, hest.ci[2,], 
               pest.ci[1,], pestimates, pest.ci[2,],
               dest.ci[1,], destimates, dest.ci[2,])
  colnames(tbl) <- c("h.low", "hawkes", "h.up", 
                     "p.low", "poisson", "p.up",
                     "d.low", "diff", "d.up")
  return(tbl)
}

ci_condit <- function(param, lambda, mu=TRUE){
  nsims <- dim(param)[1]
  nsites <- dim(param)[2]
  
  tmp <- array(dim=c(nsims, nsites))
  if(mu==FALSE){
    for(s in 1:nsims){
      for(i in 1:nsites){
        if(mean(param[s,i,])!=0){
          tmp[s,i] <- mean(param[s,i,])/mean(lambda[s,i,])
        } else{tmp[s,i] <- NA}
      }
    }
    
  }
if(mu==TRUE){
  for(s in 1:nsims){
    for(i in 1:nsites){
      if(mean(param[s,i])!=0){
        tmp[s,i] <- mean(param[s,i])/mean(lambda[s,i,])
      } else{tmp[s,i] <- NA}
    }
  }
  
}
  gestimates <- apply(tmp, 2, mean)
  gest.ci <- apply(tmp, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  tbl <- as.data.frame(cbind(gest.ci[1,], gestimates, gest.ci[2,]))
  colnames(tbl) <- c("low", "mean", "up")
  return(tbl)
}

#### 0: Load in the data ####


files <- paste("./output/community/", list.files("./output/community/"), sep="")
bad_files <- c(1, 3, 18)
files <- files[-bad_files]
sites <- sites[-1] # oops i forgot to run babcock
sites <- sites[-bad_files]
condits <- events <- list(); length(condits) <- length(events) <- length(sites)
events2 <- events

for(i in 1:length(sites)){
  
load(files[i])
hawk <- site_model$BUGSoutput$sims.list$sim_t
pois <- site_model$BUGSoutput$sims.list$sim_t2
gamma <- site_model$BUGSoutput$sims.list$gamma
lambda <- site_model$BUGSoutput$sims.list$lambda
mu <- site_model$BUGSoutput$sims.list$mu
t <- t(all_events[, which(colnames(all_events) == sites[i])])

# Drop any observations where there is no conditional intensity
# why? because then we are just comparing two poisson processes, which is silly
gamma_means <- apply(gamma, 2, mean)
if(any(gamma_means==0)){
  hawk <- hawk[,-which(gamma_means==0),]
  pois <- pois[,-which(gamma_means==0),]
  t <- t[-which(gamma_means==0),]
  gamma <- gamma[,-which(gamma_means==0),]
  lambda <- lambda[,-which(gamma_means==0),]
  mu <- mu[,-which(gamma_means==0)]
}

if(sum(gamma_means>0)>1){
  
#### 1: Comparison of total number of estimated events ####

# are they both able to estimate observed values well?
# need to load in each extract mean + ci of total events for each iteration
# compare hawkes to poisson to observed values

events[[i]] <- ci_counts(hawk, pois, t)
events2[[i]] <- ci_counts2(hawk, pois, t)

##### 2: Comparison of correspondance ####

# is the predicted event timing accurate for both processes?
# how well does the event sequence match the observed process?
# use matthews correlation coefficient to check correspondance
# extract the mean + ci for the correlation 
# need to decide on an appropriate sliding window for positive detections
# for example, a false positive is only false if not within three seconds of an event
# and a true positive is true if it is within three seconds of any event
# how do we deal with true and false negatives though? or does it not matter if we only do a directional comparison?
# we could set positive and negative windows separately
# DYNAMIC TIME WARPING!!!!
# use asymmetric with a maximum slope of 1 for flexibility
# just kidding, this is a mess for interpretation; back to MCC

# corresp <- calculate_props(hawk, pois, t)
# or maybe just skip this part?

##### 3: Contribution of conditional intensity ####

# look at the relative contribution of mu and gamma within the Hawkes models

condits[[i]] <- ci_condit(param=gamma, lambda, mu=FALSE)}
}

##### Plots ####
col1 <- RColorBrewer::brewer.pal(11,"Spectral")[1]
col2 <- RColorBrewer::brewer.pal(11,"Spectral")[10]

#evsave <- events; evsav2 <- events2

events <- events[-13]
events2 <- events2[-13]

for(i in 1:length(events)){
  if(i == 1){
    all_calls <- events[[i]]
  }
  all_calls <- rbind(all_calls, events[[i]])
}

table_1 <- as.data.frame(all_calls)
table_1 <- table_1[-which(table_1$obs<7),]
rownames(table_1) <- seq(1, nrow(table_1))
means <- apply(table_1, 2, mean)
vars <- apply(table_1, 2, var)


midpoints <- as.numeric(rbind(table_1$hawkes, table_1$poisson, rep(NA, nrow(table_1))))
li <- as.numeric(rbind(table_1$h.low, table_1$p.low, rep(NA, nrow(table_1))))
ui <- as.numeric(rbind(table_1$h.up, table_1$p.up, rep(NA, nrow(table_1))))
plot(midpoints, col=c(col1, col2, "white"), pch=c(16,15,1), cex=1, las=1, 
     xlim=c(0,nrow(table_1)*3),
     ylim=c(min(li[!is.na(li)]), max(ui[!is.na(ui)])))
abline(h=0, lty=2, lwd=2)
arrows(x0 = seq(1, nrow(table_1)*3,1), y0 = li, x1 = seq(1, nrow(table_1)*3,1), y1 = ui,
       col=c(col1, col2, "white"), angle=90, length=.03, code = 3, lwd=1)

condits <- condits[-13]
for(i in 1:length(condits)){
  if(i == 1){
    all_gams <- condits[[i]]
  }
  all_gams <- rbind(all_gams, condits[[i]])
}

table_1 <- as.data.frame(all_gams)
table_1 <- rbind(table_1, apply(table_1, 2, mean))
sitecols <- rownames(table_1)

midpoints <- as.numeric(rbind(table_1$mean, rep(NA, nrow(table_1))))
li <- as.numeric(rbind(table_1$low, rep(NA, nrow(table_1))))
ui <- as.numeric(rbind(table_1$up, rep(NA, nrow(table_1))))
plot(midpoints, col=c(col1, "white"), pch=c(16,1), cex=1, las=1, 
     xlim=c(0,nrow(table_1)*2),
     ylim=c(min(li[!is.na(li)]), max(ui[!is.na(ui)])))
abline(h=0, lty=2, lwd=2)
arrows(x0 = seq(1, nrow(table_1)*2,1), y0 = li, x1 = seq(1, nrow(table_1)*2,1), y1 = ui,
       col=c(col1, "white"), angle=90, length=.03, code = 3, lwd=1)

