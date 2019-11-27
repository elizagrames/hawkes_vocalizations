models <- list.files("./output/")
filenames <- paste("./output/", models, sep="")
filenames <- filenames[which(stringr::str_detect(filenames, "RData"))]
filenames <- filenames[-13]

for(i in 1:length(models)){
  if(i==1){sitenames <- c()}
  sitenames[i] <- strsplit(models[i], "\\.")[[1]][1]
}

for(i in 1:length(filenames)){
  assign(sitenames[i], site_model)
}
goodmods <- c()


load(filenames[i])
params <- site_model$BUGSoutput$summary
round(head(params, 20),4)
hist(params[,8])
MCMCvis::MCMCtrace(site_model, params=c("alpha", "beta", "mu"), pdf=FALSE)

goodmods <- append(goodmods, i)

i <- i +1


params <- params[,-3,]

which(stringr::str_detect(names(x), "lambda"))

i <- 7044
names(params[1,1,i])
plot(params[,1,i], type="l", col="pink")
lines(params[,2,i], col="turquoise")
mean(params[1,1,i])

round(head(params),4)
params[which(stringr::str_detect(rownames(params), "mu")),]

no1 <- MCMCvis::MCMCchains(site_model, chain_num = 1)
no2 <- MCMCvis::MCMCchains(site_model, chain_num = 2)

plot(alpha)

