library(R2jags)
options(warn=-1)
#### FUNCTIONS ####
clean_times <- function(dat, starttime = NULL) {
  for (m in 1:ncol(dat)) {
    x <- dat[, m]
    if (any(x == "")) {
      x <- x[-which(x == "")]
    }
    x <- strsplit(x, ":")
    times <- c()
    for (i in 1:length(x)) {
      times[i] <- as.numeric(x[[i]][1]) * 60 + as.numeric(x[[i]][2])
      if (!is.null(starttime)) {
        times[i] <- times[i] - starttime[m]
      }
    }
    length(times) <- nrow(dat)
    dat[, m] <- times
  }
  return(dat)
}

get_events <- function(x, int = .5, maxtime = 8 * 60) {
  times <- seq(1, maxtime, int)
  events <- as.numeric(times %in% x)
  return(events)
}
calculate_timediff <- function(events) {
  timestamps <- which(events > 0)
  timediff <- c()
  for (i in 1:length(timestamps)) {
    timediff[i] <- timestamps[i] - timestamps[i - 1]
  }
  if (any(is.na(timediff))) {
    timediff <- timediff[-is.na(timediff)]
  }
  length(timediff) <- length(events)
  timediff <- timediff[!is.na(timediff)]
  return(timediff)
}

clean_differences <- function(alltimediff) {
  for (i in 1:length(alltimediff)) {
    if (i == 1) {
      running_max <- length(alltimediff[[i]])
    }
    if (length(alltimediff[[i]]) > running_max) {
      running_max <- length(alltimediff[[i]])
    }
    
    if (i == length(alltimediff)) {
      for (k in 1:length(alltimediff)) {
        length(alltimediff[[k]]) <- running_max
        if (k == 1) {
          timediffs <- alltimediff[[k]]
        } else{
          timediffs <- cbind(timediffs, alltimediff[[k]])
        }
      }
    }
  }
  timediffs[is.na(timediffs)] <- 0
  
  return(timediffs)
}



cut_extra <- function(alltimediff) {
  cutoff <- as.numeric(quantile(alltimediff[alltimediff > 0], c(0, .5))[2])
  timecut <- alltimediff
  timecut[timecut == 0] <- NA
  timecut[timecut <= cutoff - 1] <- 1
  timecut[timecut > cutoff - 1] <- 0
  timecut[is.na(timecut)] <- 0
  return(timecut)
}


#### data clean up ####
# read in and filter data to only recordings with songs
setwd("~/hawkes_vocalizations/")
oven <-
  read.csv(
    "./data/Ovenbird point count recordings THIS ONE - Timing data.csv",
    stringsAsFactors = FALSE,
    header = FALSE
  )
oven <- oven[,-1]

oven <- oven[, which(oven[1,] != "")]
oven <- oven[, which(oven[8,] != "")]

sitename <- c()
for (i in 1:ncol(oven)) {
  sitename[i] <- paste(oven[c(2, 1), i], collapse = "_")
}

colnames(oven) <- sitename

sites <- names(which(table(sitename) > 4))

starts <- as.numeric(clean_times(oven[4,]))
oven_times <- clean_times(oven[-c(1:5),], starttime = starts)
all_events <- apply(oven_times, 2, get_events)
#alltimediff <- apply(all_events, 2, calculate_timediff)

# find the maximum number of timestamps in a recording

#alltimediff <- as.matrix(clean_differences(alltimediff))

# ignore time differences that are outside the expected range
# these are not included in the history
#timecut <- cut_extra(alltimediff)
#too_old <-
#  as.numeric(quantile(alltimediff[alltimediff > 0], c(0, .5))[2])

### turn the data from a site into useful stuff
### then run the jags model with it

##### running models and saving output ####

#clean the site data 

# make a big array of stuff

nfrags <- length(sites)
min_events <- 10
cutoff <- 60
maxmem <- 60

for(s in 1:nfrags){
  
  t <- t(all_events[, which(colnames(all_events) == sites[s])])
  i <- 1
  while(i < dim(t)[1]){
    if(sum(t[i,])<min_events){
      t <- t[-i,]
    }
    i <- i+1
  }
  
  history <- array(dim=c(dim(t)[1], dim(t)[2], maxmem))
  for(i in 1:dim(t)[1]){
    
    for(j in 1:dim(t)[2]){
      
      what_happened <- which(t[i,c(1:j)]>0)
      if(any(what_happened)>0){
        events <- j - what_happened
        memories <- events[which(events<cutoff)]
        }else{memories <- 0}
        length(memories) <- maxmem
        memories[is.na(memories)] <- 0
        
        history[i,j,] <- memories
    }
  }

jags.data <- list(
  t=t,
  tM=t,
  tP=t,
  history = history,
  maxmemory = maxmem,
  nobs = dim(t)[2],
  nsites = dim(t)[1]
)

site_model <- jags(
  data = jags.data,
  parameters.to.save = c(
    "mu",
    "lambdaH",
    "gamma",
    "alpha", 
    "beta", 
    "sim_tH",
    "lambdaM",
    "sim_tM",
    "lambdaP",
    "sim_tP"
  ),
  model.file = "./scripts/JAGS_model_community.R",
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 2000,
  n.thin = 3
)
filename <- paste("./output/Jan10/", gsub("/", "",sites[s]), ".RData", sep="")
save(site_model, file = filename)
#rm(t, history, maxmem, memories.plus, currentdiffs, memories, site_model)
  }

  # changed to .0001 per simulation and Congdon (2014)
 

# yay it runs


