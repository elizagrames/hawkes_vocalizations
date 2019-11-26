sitedat <- read.csv("../Documents/Hawkes/hipmetada.csv")
sitedat$Date <- as.character(sitedat$Date)
sitedat$Date <- as.Date(lubridate::dmy(sitedat$Date))

normalize <- function(x){
  newx <- (x-mean(x))/sd(x)
  return(newx)
}


#sitedat$Date <- normalize(sitedat$Date)
sitedat$Size <- normalize(sitedat$Size)
sitedat$Leps <- normalize(sitedat$Leps)

Babcock <- site_model 
BBS <- site_model2
BishopSwampMd <- site_model10
Cockaponset <- site_model3
Fenton <- site_model11
Hubbard <- site_model4
MML <- site_model12
#SummerLane <- site_model5
#NathanHale <- site_model6
NyeHolmanMd <- site_model14
RockevilleMd <- site_model7
RockvilleMd <- RockevilleMd
SalmonLg <- site_model8
SalmonSm <- site_model16
SalmonSmall <- SalmonSm
WhettenW <- site_model10

#all_models <- list(Babcock, BBS, BishopSwampMd, Cockaponset, Fenton, Hubbard, MML, MoonRoad, SummerLane, 
#                   NathanHale, RockvilleMd, RockvilleSmall_corrected, SalmonLg, SalmonSmall, 
#                   SleepingGiant_corrected, ValleyFalls_corrected, WhettenW)
#corrected <- c(0,0,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0)

all_models <- list(MML)

nsims <- 5000

myparams <- matrix(nrow=nsims, ncol=4)
mu.params <- myparams
alpha.params <- myparams
beta.params <- myparams
gamma.params <- myparams

mus <- c()
alphas <- c()
betas <- c()
gammas <- c()

for(i in 1:length(all_models)){
  mus[i] <- all_models[[i]]$BUGSoutput$mean$mu
  alphas[i] <- all_models[[i]]$BUGSoutput$mean$alpha
  betas[i] <- all_models[[i]]$BUGSoutput$mean$beta
  gammas[i] <- mean(all_models[[i]]$BUGSoutput$sims.list$gamma)
}

alphas[as.logical(corrected)] <- 0.0000000000000001
betas[as.logical(corrected)] <- 0.0000000000000001
gammas[as.logical(corrected)] <- 0.0000000000000001

mydat <- cbind(sitedat[7,], mus, gammas, alphas, betas)

#mydat <- mydat[which(corrected==0),]
mydat <- mydat[order(mydat$Date),]
mydat$rate <- mydat$mus + mydat$gammas

mml_raw <- read.csv("../Documents/Hawkes/Ovenbird point count recordings THIS ONE - Timing data.csv", stringsAsFactors = FALSE)
mml_raw <- t

lambda_est <- MASS::fitdistr(t[1,], "poisson")
nsims <- 500

new_process <- array(dim=c(6,961,nsims))

for(i in 1:npoints){
  for(k in 1:nsims){
    new_process[i,,k] <- plot(rpois(961, as.numeric(lambda_est[1])) ~ seq(1,961,1), type="l")
  }
  
}



mydat$Class <- as.character(mydat$Class)
mydat$Class[mydat$Class=="small"] <- "_small"
mydat$Class[mydat$Class=="medium"] <- "_xmedium"


par(las=1, pty="s", xpd=TRUE, mar=c(5,15,5,5))
boxplot(mydat$rate ~ mydat$Class, names=c("","",""),
        ylab="", tcl=.5, ylim=c(0, 0.06), xlab="Forest Fragment Size", cex.lab=2, cex=1.5)
text(1, -.005, "Small", cex=1.5)
text(2, -.005, "Medium", cex=1.5)
text(3, -.005, "Large", cex=1.5)

text(-.5, .035, expression(lambda), cex=5)
text(-.5, .0275, "(overall\n rate)", cex=2)

plot(mydat$rate ~ mydat$Leps, ylab="", 
     xlab="Relative caterpillar abundance", pch=19, cex=1.5, xlim=c(-2, 2), tcl=.5, ylim=c(0, 0.06))
tmp <- glm(rate ~ Leps, data=mydat, family=Gamma(link="inverse"))
tmp2 <- lm(rate ~ Leps, data=mydat)
pred <- predict(tmp, se.fit = TRUE, type = "response")
segments(-2.15, coef(tmp2)[2]*-2.15+coef(tmp2)[1], 2.15, coef(tmp2)[2]*2.15+coef(tmp2)[1], lwd=2, lty=2)
# abline(lm(pred$fit ~ mydat$Leps), lwd=2, lty=2)
text(mydat$Leps[1]-3.5, .03, expression(lambda), cex=5)


ci <- matrix( c(pred$fit + 1.96 * pred$se.fit, pred$fit - 1.96 * 
                  pred$se.fit), ncol=2 ) 
plot(mydat$rate ~ mydat$Leps, ylim=c(0,.06), pch=20, cex=1.5, xlab="Relative caterpillar abundance", 
     ylab = c(expression(lambda)))
#lines(as.numeric(pred$fit) ~ mydat$Leps, lwd=2, lty=2)
polygon(c(mydat$Leps, rev(mydat$Leps)), c(ci[,1], rev(ci[,2])), col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Leps, lwd=2, lty=2)
points(mydat$rate ~ mydat$Leps, ylim=c(0,.06), pch=20, cex=2, xlab="Date", ylab = "mu")
box()



modmu <- glm(mus ~ Date, data=mydat, family=Gamma(link="inverse"))
modalpha <- glm(alphas ~ Date, data=mydat, family=Gamma(link="inverse"))
modbeta <- glm(betas ~ Date, data=mydat, family=Gamma(link="inverse"))
modgamma <- glm(gammas ~ Date, data=mydat, family=Gamma(link="inverse"))
modfull <- glm(rate ~ Date, data=mydat, family=Gamma(link="inverse"))
tmp <- glm(rate ~ Class, data=mydat, family=Gamma(link="inverse"))
summary(tmp)

summary(modmu)
summary(modalpha)
summary(modbeta)
summary(modgamma)
summary(modfull)

par(las=1, xpd=TRUE, mar=c(5,10,10,1))

pred <- predict(modmu, se.fit = TRUE, type = "response")
ci <- matrix( c(pred$fit + 1.96 * pred$se.fit, pred$fit - 1.96 * 
                  pred$se.fit), ncol=2 ) 
ci[ci>.0625] <- 0.0625
ci[ci<0] <- -0.0025
plot(mydat$mus ~ mydat$Date, ylim=c(0,.06),pch=20, cex=1.5, xlab="Date", cex.lab=2,
     ylab = "", tcl=.5)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
polygon(c(mydat$Date, rev(mydat$Date)), c(ci[,1], rev(ci[,2])), col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
points(mydat$mus ~ mydat$Date, ylim=c(0,.06), pch=20, cex=2.5, xlab="Date", ylab = "mu")
text(mydat$Date[1]-18, .035, expression(mu), cex=5)
text(mydat$Date[1]-18, .0275, "(background\nrate)", cex=2)
box()

pred <- predict(modgamma, se.fit = TRUE, type = "response")
ci <- matrix( c(pred$fit + 1.96 * pred$se.fit, pred$fit - 1.96 * 
                  pred$se.fit), ncol=2 ) 
ci[ci>.0625] <- 0.0625
ci[ci<0] <- -0.0025
plot(mydat$gammas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=1.5, xlab="Date",  cex.lab=2, ylab = "", tcl=.5)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
polygon(c(mydat$Date, rev(mydat$Date)), c(ci[,1], rev(ci[,2])), col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
points(mydat$gammas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=2)
text(mydat$Date[1]-18, .035, expression(gamma), cex=5)
text(mydat$Date[1]-18, .0275, "(conditional\nintensity)", cex=2)
box()


pred <- predict(modalpha, se.fit = TRUE, type = "response")
ci <- matrix( c(pred$fit + 1.96 * pred$se.fit, pred$fit - 1.96 * 
                  pred$se.fit), ncol=2 ) 
ci[ci>.0625] <- 0.0625
ci[ci<0] <- -0.0025
plot(mydat$alphas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=1.5, xlab="Date", ylab = "", tcl=.5, cex.lab=2)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
polygon(c(mydat$Date, rev(mydat$Date)), c(ci[,1], rev(ci[,2])), col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
points(mydat$alphas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=2)
text(mydat$Date[1]-19, .035, expression(alpha), cex=5)
text(mydat$Date[1]-19, .0275, "(self-excitement\nparameter)", cex=2)
box()



pred <- predict(modbeta, se.fit = TRUE, type = "response")
ci <- matrix( c(pred$fit + 2.576 * pred$se.fit, pred$fit - 2.576 * 
                  pred$se.fit), ncol=2 ) 
ci[ci>.0625] <- 0.0625
ci[ci<0] <- -0.0025
plot(mydat$betas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=1.5, xlab="Date", ylab = "", tcl=.5, cex.lab=2)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
polygon(c(mydat$Date, rev(mydat$Date)), c(ci[,1], rev(ci[,2])), col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Date, lwd=2, lty=2)
points(mydat$betas ~ mydat$Date, ylim=c(0,.06), pch=20, cex=2)
text(mydat$Date[1]-18, .035, expression(beta), cex=5)
text(mydat$Date[1]-18, .029, "(decay rate)", cex=2)
box()



pred <- predict(modfull, se.fit = TRUE, type="response")
ci <- matrix( c(pred$fit + 1.96 * pred$se.fit, pred$fit - 1.96 * 
                  pred$se.fit), ncol=2 ) 
ci[ci>.0625] <- 0.0625
ci[ci<0] <- -0.0025
plot(mydat$rate ~ mydat$Date, type="n", xlab="Date", ylim=c(0,.06), cex.lab=2, 
     ylab = "", tcl=.5)
lines(pred$fit ~ mydat$Date, lwd=2)
polygon(c(mydat$Date, rev(mydat$Date)), c(ci[,1], rev(ci[,2])), 
        col="#afeeee", border=NA)
lines(pred$fit ~ mydat$Date, lwd=2)
points(mydat$rate ~ mydat$Date, ylim=c(0,.06), pch=19, cex=2)
text(mydat$Date[1]-18, .035, expression(lambda), cex=5)
text(mydat$Date[1]-18, .03, "(overall rate)", cex=2)
box()


for(n in 1:nsims){
  
  tryCatch({
    mus <- c()
    alphas <- c()
    betas <- c()
    gammas <- c()
    
    for(i in 1:length(all_models)){
      mus[i] <- sample(all_models[[i]]$BUGSoutput$sims.list$mu,1)
      alphas[i] <- sample(all_models[[i]]$BUGSoutput$sims.list$alpha,1)
      betas[i] <- sample(all_models[[i]]$BUGSoutput$sims.list$beta,1)
      number <- sample(1:2400,1)
      gammas[i] <- mean(all_models[[i]]$BUGSoutput$sims.list$gamma[number,,])
    }
    
    alphas[as.logical(corrected)] <- 0.000000001
    betas[as.logical(corrected)] <- 0.0000000001
    gammas[as.logical(corrected)] <- 0.0000000001
    
    mydat <- cbind(sitedat, mus, gammas, alphas, betas)
    modmu <- glm(mus ~ Size + Leps + Date, data=mydat, family=Gamma(link="identity"))
    modalpha <- glm(alphas ~ Size + Leps + Date, data=mydat)
    modbeta <- glm(betas ~ Size + Leps + Date, data=mydat)
    modgamma <- glm(gammas ~ Size + Leps + Date, data=mydat)
    
    mu.params[n,] <- coefficients(modmu)
    alpha.params[n,] <- coefficients(modalpha)
    beta.params[n,] <- coefficients(modbeta)
    gamma.params[n,] <- coefficients(modgamma)
  }, error=function(e){})
  
  
  }

## MU #####
simparams <- mu.params
removals <- c()
for(i in 1:nrow(simparams)){
  if(any(is.na(simparams[i,]))){
    removals <- append(removals, i)  
  }
  
  if(i == nrow(simparams)){
    simparams <- simparams[-removals,]
  }
}

mu.params <- simparams
mu_estimates <- apply(mu.params, 2, mean)

for(i in 1:4){
  ci <- quantile(mu.params[,i], c(.025, .975))
  entry <- cbind(ci[1], mu_estimates[i], ci[2])
  if(i ==1){paramdata <- entry}
  if(i >1){paramdata <- rbind(paramdata, entry)}
  if(i == 4){
    paramdata <- as.data.frame(paramdata)
    rownames(paramdata) <- c("int", "size", "leps", "date")
    colnames(paramdata) <- c("lci", "est", "uci")
  }
}
mudata <- paramdata

## alpha #####
simparams <- alpha.params
removals <- c()
for(i in 1:nrow(simparams)){
  if(any(is.na(simparams[i,]))){
    removals <- append(removals, i)  
  }
  
  if(i == nrow(simparams)){
    simparams <- simparams[-removals,]
  }
}

alpha.params <- simparams
alpha_estimates <- apply(alpha.params, 2, mean)

for(i in 1:4){
  ci <- quantile(alpha.params[,i], c(.025, .975))
  entry <- cbind(ci[1], alpha_estimates[i], ci[2])
  if(i ==1){paramdata <- entry}
  if(i >1){paramdata <- rbind(paramdata, entry)}
  if(i == 4){
    paramdata <- as.data.frame(paramdata)
    rownames(paramdata) <- c("int", "size", "leps", "date")
    colnames(paramdata) <- c("lci", "est", "uci")
  }
}
alphadata <- paramdata


## beta #####
simparams <- mu.params
removals <- c()
for(i in 1:nrow(simparams)){
  if(any(is.na(simparams[i,]))){
    removals <- append(removals, i)  
  }
  
  if(i == nrow(simparams)){
    simparams <- simparams[-removals,]
  }
}

mu.params <- simparams
mu_estimates <- apply(mu.params, 2, mean)

for(i in 1:4){
  ci <- quantile(mu.params[,i], c(.025, .975))
  entry <- cbind(ci[1], mu_estimates[i], ci[2])
  if(i ==1){paramdata <- entry}
  if(i >1){paramdata <- rbind(paramdata, entry)}
  if(i == 4){
    paramdata <- as.data.frame(paramdata)
    rownames(paramdata) <- c("int", "size", "leps", "date")
    colnames(paramdata) <- c("lci", "est", "uci")
  }
}
mudata <- paramdata


plot(ksmooth(predict(modmu), as.numeric(names(predict(modmu))), "normal", bandwidth = 10), col = 2)

