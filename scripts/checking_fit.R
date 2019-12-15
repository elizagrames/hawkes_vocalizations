files <- paste("./output/community/", list.files("./output/community/"), sep="")
i <- 0

par(mfrow=c(1,1))
i <- i+1
load(files[i])
params <- site_model$BUGSoutput$summary
round(head(params, 10),4)
type <- substr(rownames(params), 1, 5)
not_these <- which(type=="sim_t"| type=="gamma")
hist(params[-not_these,8], main=files[i])
MCMCvis::MCMCtrace(site_model, params=c("alpha", "beta", "mu", "lambda2"), pdf=FALSE)

params[which(params[-not_these,8]>1.1),8]
dev.off()


## BAD MODELS
# dropping BBL because mu is verging on bimodal
# dropping BSM because rhat is all over the place and chains have not converged visually for mu, alpha, or beta
# dropping SGL because wow those chains are all over the place
bad_files <- c(1, 3, 18)

##### ESTIMATED EVENTS ######
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

sites <- sites[-1] # oops i forgot to run babcock

sites <- sites[-bad_files]
files <- paste("./output/community/", list.files("./output/community/"), sep="")
files <- files[-bad_files]
i <- 0

i <- i+1

load(files[i])
hawk <- site_model$BUGSoutput$sims.list$sim_t
pois <- site_model$BUGSoutput$sims.list$sim_t2
t <- t(all_events[, which(colnames(all_events) == sites[i])])


table_1 <- as.data.frame(ci_counts(hawk, pois, t))
table_1 <- table_1[order(table_1$obs),]

table_1$obs <- rep(0,  nrow(table_1))
means <- apply(table_1, 2, mean)

col1 <- RColorBrewer::brewer.pal(11,"Spectral")[1]
col2 <- RColorBrewer::brewer.pal(11,"Spectral")[10]

midpoints <- as.numeric(rbind(table_1$hawkes, table_1$poisson, rep(NA, nrow(table_1))))
li <- as.numeric(rbind(table_1$h.low, table_1$p.low, rep(NA, nrow(table_1))))
ui <- as.numeric(rbind(table_1$h.up, table_1$p.up, rep(NA, nrow(table_1))))
plot(midpoints, col=c(col1, col2, "white"), pch=c(16,15,1), cex=2, las=1, xlim=c(0,nrow(table_1)*3+2),
     ylim=c(min(li[!is.na(li)]), max(ui[!is.na(ui)])))
abline(h=0, lty=2, lwd=3)
arrows(x0 = seq(1, nrow(table_1)*3,1), y0 = li, x1 = seq(1, nrow(table_1)*3,1), y1 = ui,
       col=c(col1, col2, "white"), angle=90, length=.07, code = 3, lwd=2)
points(nrow(table_1)*3+1,means['hawkes'], pch=16,cex=2, col=col1)
points(nrow(table_1)*3+2,means['poisson'], pch=15,cex=2, col=col2)

arrows(x0 = nrow(table_1)*3+1, y0 = means['h.low'], x1 = nrow(table_1)*3+1, 
       y1 = means['h.up'],
       col=col1, angle=90, length=.07, code = 3, lwd=2)

arrows(x0 = nrow(table_1)*3+2, y0 = means['p.low'], x1 = nrow(table_1)*3+2, 
       y1 = means['p.up'],
       col=col2, angle=90, length=.07, code = 3, lwd=2)

#abline(v=seq(3, 23.5, 3), lty=2)


##### CORRESPONDANCE ####
mcc <- function(tp, tn, fp, fn){
  ((tp*tn) - (fp*fn))/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
}

hawk <- site_model$BUGSoutput$sims.list$sim_t
pois <- site_model$BUGSoutput$sims.list$sim_t2
t <- t(all_events[, which(colnames(all_events) == sites[i])])


calculate_props <- function(hawk, pois, t){
  nsims <- dim(hawk)[1]
  nsites <- 2
  nobs <- dim(hawk)[3]
  
  tmp.h <- tmp.p <- diffph <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 2:nsites){
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

z3 <- calculate_props(hawk, pois, t)
table_1 <- as.data.frame(z3)
means <- apply(table_1, 2, mean)

midpoints <- as.numeric(rbind(table_1$hawkes, table_1$poisson, rep(NA, nrow(table_1))))
li <- as.numeric(rbind(table_1$h.low, table_1$p.low, rep(NA, nrow(table_1))))
ui <- as.numeric(rbind(table_1$h.up, table_1$p.up, rep(NA, nrow(table_1))))
plot(midpoints, col=c(col1, col2, "white"), pch=c(16,15, 1), 
     cex=2, las=1, xlim=c(0,nrow(table_1)*3+2), ylim=c(-.01,.5))
abline(h=0, lty=2, lwd=3)
arrows(x0 = seq(1, nrow(table_1)*3,1), y0 = li, x1 = seq(1, nrow(table_1)*3,1), y1 = ui,
       col=c(col1, col2, "white"), angle=90, length=.07, code = 3, lwd=2)
points(nrow(table_1)*3+1,means['hawkes'], pch=16,cex=2, col=col1)
points(nrow(table_1)*3+2,means['poisson'], pch=15,cex=2, col=col2)

arrows(x0 = nrow(table_1)*3+1, y0 = means['h.low'], x1 = nrow(table_1)*3+1, 
       y1 = means['h.up'],
       col=col1, angle=90, length=.07, code = 3, lwd=2)

arrows(x0 = nrow(table_1)*3+2, y0 = means['p.low'], x1 = nrow(table_1)*3+2, 
       y1 = means['p.up'],
       col=col2, angle=90, length=.07, code = 3, lwd=2)


############# MU vs GAMMA #####

# doing good at estimating the total number in both cases
# doing okay-ish with the times in some cases
# what exta info do we get from hawkes?
# for some sites, we can only fit poisson because gamma has to be fixed to zero when there is no history
# only makes sense to compare mu and gamma when gamma is not forcibly fixed by me

alpha <- site_model$BUGSoutput$sims.list$alpha
beta <- site_model$BUGSoutput$sims.list$beta
mu <- site_model$BUGSoutput$sims.list$mu
gamma <- site_model$BUGSoutput$sims.list$gamma
lambda <- site_model$BUGSoutput$sims.list$lambda

condit <- c()
for(p in 1:dim(mu)[2]){
  if(mean(gamma[,p,])!=0){
    condit[p] <- mean(gamma[,p,])/mean(lambda[,p,])
  } else{condit[p] <- NA}
}

plot(density(condit[!is.na(condit)]))

ci_condit <- function(gamma, lambda){
  nsims <- dim(gamma)[1]
  nsites <- dim(gamma)[2]
  
  tmp <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      if(mean(gamma[s,i,])!=0){
        tmp[s,i] <- mean(gamma[s,i,])/mean(lambda[s,i,])
      } else{tmp[s,i] <- NA}
    }
  }
  
  gestimates <- apply(tmp, 2, mean)
  gest.ci <- apply(tmp, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  tbl <- as.data.frame(cbind(gest.ci[1,], gestimates, gest.ci[2,]))
  colnames(tbl) <- c("low", "mean", "up")
  return(tbl)
}

table_c <- ci_condit(gamma, lambda)
table_c <- table_c[-which(is.na(table_c[,1])),]











