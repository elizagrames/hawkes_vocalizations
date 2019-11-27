files <- paste("./output/", list.files("./output/"), sep="")

i <- i+1
load(files[i])
params <- site_model$BUGSoutput$summary
round(head(params, 20),4)
hist(params[,8])
MCMCvis::MCMCtrace(site_model, params=c("alpha", "beta", "mu"), pdf=FALSE)

#bad models
#beaver brook large
#rockville md is really bad
#rockville small is questionable beta and alpha spike to .4 sometimes but rhat okay
#something went wrong with rockville small gamma[1,everything]
#nye holman md is a bit wonky, but rhat okay

