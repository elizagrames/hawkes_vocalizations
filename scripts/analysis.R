#### -1: Functions ####
ci_counts <- function(mod, truth){
  nsims <- dim(mod)[1]
  nsites <- dim(mod)[2]
  
  tmp.mod <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      tmp.mod[s,i] <- length(which(mod[s,i,]>0)) - sum(truth[i,])
    }
  }
  
  mod.mean <- apply(tmp.mod, 2, mean)
  mod.ci <- apply(tmp.mod, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(truth)
  tbl <- cbind(counts, mod.ci[1,], mod.mean, mod.ci[2,])
  colnames(tbl) <- c("obs", "li", "x", "ui")
  output <- list(tbl, tmp.mod)
  return(output)
}

mcc <- function(tp, tn, fp, fn){
  ((tp*tn) - (fp*fn))/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
}

calculate_props <- function(mod, truth, size){
  nsims <- dim(mod)[1]
  nsites <- dim(mod)[2]
  nobs <- dim(mod)[3]
  
  tmp <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      obs.e <- which(truth[i,]>0)
      window <- obs.e
      for(m in 1:length(obs.e)){
        window <- append(window, obs.e[m]+size)
        window <- append(window, obs.e[m]-size)
      }
      mod.pos <- which(mod[s,i,]>0)
      mod.neg <- which(mod[s,i,]==0)
      
      obs.pos <- obs.e
      obs.neg <- which(truth[i,]==0)
      
      tp <- sum(hawk.pos %in% window) 
      tn <- sum(hawk.neg %in% obs.neg)
      fp <- sum((hawk.pos %in% window)!=TRUE)
      fn <- sum(hawk.neg %in% obs.pos)
      
      tmp[s,i] <- mcc(tp, tn, fp, fn)
    }
  }
  tmp[is.na(tmp)] <- NA

  mod.mean <- apply(tmp.mod, 2, mean)
  mod.ci <- apply(tmp.mod, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(truth)
  tbl <- cbind(mod.ci[1,], mod.mean, mod.ci[2,])
  colnames(tbl) <- c("li", "x", "ui")
  output <- list(tbl, tmp)
  return(output)
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
  colnames(tbl) <- c("li", "x", "ui")
  output <- list(tbl, tmp)
  return(output)
}

#### 0: Load in the data ####

models <- list.files("./output/community/")
filenames <- paste("./output/community/", models, sep="")
filenames <- filenames[which(stringr::str_detect(filenames, "RData"))]
#filenames <- filenames[-13]

for(i in 1:length(models)){
  if(i==1){sites <- c()}
  sites[i] <- strsplit(models[i], "\\.")[[1]][1]
}

params <- c("mu",
            "alpha", 
            "beta", 
            "lambdaM",
            "lambdaP")

all_events <- read.csv("all_events.csv")

for(z in 2:ncol(all_events)){
name <- substring(colnames(all_events)[z], nchar(colnames(all_events)[z])-1, nchar(colnames(all_events)[z]))
if(grepl("\\.", name)){colnames(all_events)[z] <- substring(colnames(all_events)[z], 1, nchar(colnames(all_events)[z])-2)}
}

colnames(all_events) <- gsub("\\.", " ", colnames(all_events))
library(tm)
for(i in 1:length(filenames)){
  load(filenames[i])
  estimates <- site_model$BUGSoutput$summary
  types <- removePunctuation(removeNumbers(rownames(estimates)))
  rhats <- estimates[which(types %in% params), 'Rhat']
  if(any(rhats>1.1)){
    tag <- "bad"
  }else{
    tag <- "good"
  }
  write.csv(estimates, file=paste("./", sites[i], tag, ".csv", sep=""))

    if(tag=="good"){
    truth <- t(all_events[, which(colnames(all_events) == sites[i])])
    eventsH <- ci_counts(site_model$BUGSoutput$sims.list$sim_tH, truth)
    eventsM <- ci_counts(site_model$BUGSoutput$sims.list$sim_tM, truth)
    eventsP <- ci_counts(site_model$BUGSoutput$sims.list$sim_tP, truth)
    eventsW <- ci_counts(site_model$BUGSoutput$sims.list$sim_tW, truth)
    conditsH <- ci_condit(param=site_model$BUGSoutput$sims.list$gamma, site_model$BUGSoutput$sims.list$lambdaH, mu=FALSE)
    conditsW <- ci_condit(param=site_model$BUGSoutput$sims.list$gammaW, site_model$BUGSoutput$sims.list$lambdaW, mu=FALSE)
mccH1 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tH, truth, 1)
mccH2 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tH, truth, 2)
mccH3 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tH, truth, 3)

mccP1 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tP, truth, 1)
mccP2 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tP, truth, 2)
mccP3 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tP, truth, 3)

mccM1 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tM, truth, 1)
mccM2 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tM, truth, 2)
mccM3 <- calculate_props(site_model$BUGSoutput$sims.list$sim_tM, truth, 3)

    
    output <- list(eventsH, eventsM, eventsP, eventsW, conditsH, conditsW, mccH1, mccH2, mccH3, mccP1, mccP2, mccP3, mccM1, mccM2, mccM3)
   names( output) <- c("eventsH", "eventsM", "eventsP", "eventsW", "conditsH", "conditsW", "mccH1", "mccH2", "mccH3", "mccP1", "mccP2", "mccP3", "mccM1","mccM2", "mccM3")
    save(output, file=paste("./", sites[i], "output_for_plots.RData", sep=""))
    rm(truth, eventsH, eventsM, eventsP, eventsW, conditsH, conditsW)
    }
}


#### 1: Comparison of total number of estimated events ####

# are they both able to estimate observed values well?
# need to load in each extract mean + ci of total events for each iteration
# compare hawkes to poisson to observed values

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


