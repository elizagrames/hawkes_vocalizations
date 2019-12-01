rm(list = ls())
library(R2jags)

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

find_marbles <-  function(history,differences,min_events = 19,cutoff = 100) {
    for (i in 1:nrow(history)) {
      if (i == 1) {
        memory_list <- list()
        length(memory_list) <- ncol(history)
        memories <- list()
        length(memories) <- nrow(history)
      }
      for (j in 1:ncol(history)) {
        temp <- history[i, c(1:j)]
        memory <- temp
        for (n in 1:length(temp)) {
          if (temp[n] != 0) {
            memory[n] <- j - temp[n]
          } else{
            memory[n] <- 0
          }
        }
        if (any(memory == 0)) {
          memory <- memory[-which(memory == 0)]
        }
        if (length(memory) == 0) {
          memory <- 0
        }
        if (any(rowSums(history) > min_events) == FALSE) {
          memory <- 0
        }
        mindiff <- min(differences[differences[, i] > 0, i])
        if (mindiff > cutoff) {
          memory <- 0
        }
        memory_list[[j]] <- memory
      }
      
      memories[[i]] <- memory_list
    }
    return(memories)
  }

find_max <- function(memories) {
  maxmemory <- 0
  for (i in 1:length(memories)) {
    memory_list <- memories[[i]]
    for (j in 1:length(memory_list)) {
      maxmemory <- max(append(maxmemory, length(memory_list[[j]])))
    }
  }
  return(maxmemory)
}

inflate_memories <- function(memories, maxmemory) {
  for (i in 1:length(memories)) {
    memory_list <- memories[[i]]
    for (j in 1:length(memory_list)) {
      if (length(memory_list[[j]]) < maxmemory) {
        length(memory_list[[j]]) <- maxmemory
      }
      memory_list[[j]][which(is.na(memory_list[[j]]))] <- 0
    }
    memories[[i]] <- memory_list
  }
  return(memories)
}

make_history <- function(history, memories, maxmemory) {
  newmat <- matrix(nrow = nrow(history), ncol = ncol(history))
  all_mats <- c()
  for (m in 1:maxmemory) {
    for (i in 1:length(memories)) {
      for (j in 1:ncol(newmat))
        newmat[i, j] <- memories[[i]][[j]][m]
    }
    all_mats <- append(all_mats, newmat)
  }
  
  
  myarray <- array(all_mats, dim = c(nrow(history), ncol(history), maxmemory))
  
  history <- myarray
  if (any(is.na(history))) {
    history[is.na(history)] <- 0
  }
  
  return(history)
}


#### data clean up ####
# read in and filter data to only recordings with songs
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
alltimediff <- apply(all_events, 2, calculate_timediff)



# find the maximum number of timestamps in a recording

alltimediff <- as.matrix(clean_differences(alltimediff))

# ignore time differences that are outside the expected range
# these are not included in the history
timecut <- cut_extra(alltimediff)
too_old <-
  as.numeric(quantile(alltimediff[alltimediff > 0], c(0, .5))[2])

### turn the data from a site into useful stuff
### then run the jags model with it

##### running models and saving output ####

#clean the site data 

# make a big array of stuff

nfrags <- length(sites)


for(s in 1:nfrags){
  
  t <- t(all_events[, which(colnames(all_events) == sites[s])])
  
  
  currentdiffs <-
    alltimediff[, which(colnames(oven_times) == sites[s])]
  
  memories <- find_marbles(t,
                           differences = currentdiffs,
                           min_events = 19,
                           cutoff = too_old)
  maxmem <- find_max(memories)

  memories.plus <- inflate_memories(memories, maxmem)
history <- make_history(t, memories.plus, maxmem)

jags.data <- list(
  t = t,
  t2=t,
  history = history,
  maxmemory = maxmem,
  nobs = dim(t)[2],
  nsites = dim(t)[1]
)

site_model <- jags.parallel(
  data = jags.data,
  parameters.to.save = c(
    "sim_t",
    "lambda",
    "lambda2",
    "sim_t2",
    "alpha",
    "beta",
    "gamma",
    "mu"
  ),
  model.file = "./scripts/JAGS_model.R",
  n.chains = 3,
  n.iter = 3000,
  n.burnin = 1000,
  n.thin = 3
)
filename <- paste("./output/", gsub("/", "",sites[s]), ".RData", sep="")
save(site_model, file = filename)
rm(t, history, maxmem, memories.plus, currentdiffs, memories, site_model)
  }

  # changed to .0001 per simulation and Congdon (2014)
 

# yay it runs


