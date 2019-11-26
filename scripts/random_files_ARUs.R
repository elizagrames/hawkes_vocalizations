babcock <- list.files("../Hawkes/Babcock/")
greathollow <- list.files("../Hawkes/GreatHollow/")

bab_files <- c()
for(i in 1:length(babcock)){
  folder <- babcock[i]
  x <- paste("../Hawkes/Babcock/", folder, sep="")
  internals <- list.files(x)
  filelist <- paste(x, "/", internals, sep="")
  all_files <- append(all_files, filelist)
}

all_files <- c()
for(i in 1:length(greathollow)){
  folder <- greathollow[i]
  x <- paste("../Hawkes/GreatHollow/", folder, sep="")
  internals <- list.files(x)
  filelist <- paste(x, "/", internals, sep="")
  all_files <- append(all_files, filelist)
}

nsample <- 50

sample(1:1000, 1)
# 530
set.seed(530)
file_samples <- sample(1:length(all_files), nsample, replace = FALSE)
hist(file_samples, 10)

files <- all_files[file_samples]

for(i in 1:length(files)){
  endfile <- paste(strsplit(files[i], "/")[[1]][c(3,5)], collapse="/")
  
  newpath <- paste("../Hawkes/sample/", endfile, sep="")
  file.copy(files[i], newpath)
}
