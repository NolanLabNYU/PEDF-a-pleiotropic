##rf hyperparameter search/tune
##George Crowley

start.time.total <- Sys.time()

set.seed(123)


##load required packages and functions
library(haven)
library(randomForest)
library(matrixStats)
library(foreach)
library(doParallel)
library(doRNG)
source("ham.dist.R")
source("uniqueEl.R")
source("parseSlurmNodeList.R")
source("fix.peak.labels.R")


#load data
##contains 693 variables for 30 cases
dat.qual <- read.csv("./dat/RF1 Data.csv", stringsAsFactors = FALSE, header = TRUE)
class.label <- "Case_Group"

#record which variables made it into qualified rf
met_info <- data.frame("all_mets" = setdiff(names(dat.qual), class.label))
met_info$qualified <- met_info$all_mets %in% setdiff(colnames(dat.qual), class.label)

##make formula for rf
rf.formula <- as.formula(paste("as.factor(Case_Group) ~", paste(setdiff(names(dat.qual), class.label), collapse = "+"), sep = ""))

##define tune parameters
ntrial <- 100
ntrial.per.param <- 10
n.qual <- sum(met_info$qualified)
ref.prof.size <- round(0.05*n.qual)

##set up parallel backend and register it with foreach package
nodes <- parseSlurmNodeList(Sys.getenv("SLURM_NODELIST"))
cpus.per.node <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
where.to.make.cluster <- unlist(lapply(nodes, rep, times=cpus.per.node))
cluster <- makePSOCKcluster(where.to.make.cluster)
registerDoParallel(cluster)

##grow all random forests using random hyperparameter search method
start.time.q <- Sys.time()
forest.results <- foreach(i=1:ntrial, .combine = 'c', .multicombine = TRUE, .packages = "randomForest", .options.rng = 456) %dorng% {
  forest.result <- list("mtry" = floor(sqrt(n.qual)),
                        "ntree" = round(5000*i),
                        "nodesize" = 1,
                        "sampsize" = c(-999, -999), ##set to -999 because the parameter is unused by call to randomForest()
                        "forests" = vector(mode = "list", length = ntrial.per.param))
  
  for (j in 1:ntrial.per.param) {
    forest <- randomForest(rf.formula, data = dat.qual, 
                           ntree = forest.result[["ntree"]], mtry = forest.result[["mtry"]],
                           #sampsize = forest.result[["sampsize"]],
                           nodesize = forest.result[["nodesize"]],
                           importance = TRUE, keep.forest = TRUE, keep.inbag = TRUE)
    forest.ref.prof <- rownames(forest$importance[order(forest$importance[,"MeanDecreaseAccuracy"], decreasing = TRUE),])
    forest.ref.prof.mda <- forest$importance[order(forest$importance[,"MeanDecreaseAccuracy"], decreasing = TRUE),"MeanDecreaseAccuracy"]
    forest.oob.error <- mean(forest$confusion[,"class.error"])
    forest.result[["forests"]][[j]] <- list("forest.ref.prof" = forest.ref.prof, 
                                            "forest.ref.prof.mda" = forest.ref.prof.mda,
                                            "forest.oob.error" = forest.oob.error,
                                            "randomForest.object" = NULL
                                            )
    if(5000*i == 345000){
      #saving the forest of what was determined to be the best forest by running 
      #this script once before determining which is the best iteration (to save entire forests), helps with peak memory usage
      forest.result[["forests"]][[j]][["randomForest.object"]] <- forest
    }
  }
  return(forest.result)
}
end.time.q <- Sys.time()
stopCluster(cluster)

#save results
saveRDS(forest.results, "./res/RandomHyperparamSearchQualified.RDS")

##get forest.results into nested structure that is compatible with code below
name <- c("mtry", "ntree", "nodesize", "sampsize", "forests")
forest.results.2 <- vector(mode = "list", length = ntrial)
for (i in 0:(ntrial-1)){
  forest <- vector(mode="list", length = 0)
  which.elems <- 1 + (5*i):(4+5*i)
  for (j in 1:5){
    forest[[name[[j]]]] <- forest.results[[which.elems[[j]]]]
  }
  forest.results.2[[i+1]] <- forest
}

forest.results <- forest.results.2

##set new seed so I can run the next section w/o running previous section
set.seed(789)
##get refined profiles and oob error from forest.results (all the rf results) and put them in arrays 
ref.prof.all <- array(dim = c(n.qual, ntrial.per.param, ntrial))
ref.prof.mda.all <- array(dim = c(n.qual, ntrial.per.param, ntrial))
oob.error.all <- array(dim = c(1,ntrial.per.param, ntrial))

for (i in 1:ntrial) {
  for (j in 1:ntrial.per.param){
    ref.prof.all[,j,i] <- forest.results[[i]][["forests"]][[j]][["forest.ref.prof"]]
    ref.prof.mda.all[,j,i] <- forest.results[[i]][["forests"]][[j]][["forest.ref.prof.mda"]]
    oob.error.all[,j,i] <- forest.results[[i]][["forests"]][[j]][["forest.oob.error"]]
  }
}



##calculate hamming distances
##calculate number of unique elements
##compute mean oob error rate for each set of parameters
ham.array <- array(dim = c(ntrial.per.param, ntrial.per.param, ntrial))
unique.array <- array(dim = c(ntrial.per.param, ntrial.per.param, ntrial))
oob.array <- array(dim = c(2, ntrial))

for (i in 1:ntrial) {
  ham.array[,,i] <- ham.dist(ref.prof.all[1:ref.prof.size,,i])
  unique.array[,,i] <- uniqueEl(ref.prof.all[1:ref.prof.size,,i])
  
  oob.array[1,i] <- mean(oob.error.all[1,,i])
  oob.array[2,i] <- sd(oob.error.all[1,,i])
}


##calculate average value of hamming distance and unique elements at each trial
average.n.unique.elem.per.param <- vector(mode = "numeric", length = ntrial)
average.hamming.dist.per.param <- vector(mode = "numeric", length = ntrial)
sd.n.unique.elem.per.param <- vector(mode = "numeric", length = ntrial)
sd.hamming.dist.per.param <- vector(mode = "numeric", length = ntrial)

for (i in 1:ntrial) {
  n.unique.elem.per.param <- unique.array[,,i]
  hamming.dist.per.param <- ham.array[,,i]
  average.n.unique.elem.per.param[i] <- mean(n.unique.elem.per.param[upper.tri(n.unique.elem.per.param, diag = FALSE)])
  average.hamming.dist.per.param[i] <- mean(hamming.dist.per.param[upper.tri(hamming.dist.per.param, diag = FALSE)])
  sd.n.unique.elem.per.param[i] <- sd(n.unique.elem.per.param[upper.tri(n.unique.elem.per.param, diag = FALSE)])
  sd.hamming.dist.per.param[i] <- sd(hamming.dist.per.param[upper.tri(hamming.dist.per.param, diag = FALSE)])
  
}

##also get tuning parameters at each trial (used for plotting results of tuning procedure)
mtry.per.trial <- vector(mode = "numeric", length = ntrial)
ntree.per.trial <- vector(mode = "numeric", length = ntrial)
nodesize.per.trial <- vector(mode = "numeric", length = ntrial)
sampsize.per.trial <- array(dim = c(length(forest.results[[1]][["sampsize"]]), ntrial))

for (i in 1:ntrial) {
  mtry.per.trial[i] <- forest.results[[i]][["mtry"]]
  ntree.per.trial[i] <- forest.results[[i]][["ntree"]]
  nodesize.per.trial[i] <- forest.results[[i]][["nodesize"]]
  sampsize.per.trial[,i] <- forest.results[[i]][["sampsize"]]
}

##save tune results in tabular format for plotting (in matlab probably)
tune.results <- data.frame(mtry.per.trial, ntree.per.trial, nodesize.per.trial, sampsize.per.trial,
                           average.hamming.dist.per.param, sd.hamming.dist.per.param,
                           average.n.unique.elem.per.param, sd.n.unique.elem.per.param)
write.csv(tune.results, "./res/tune.results.csv")

##find parameter set that minimizes ham.dist and/or unique.array

##for ham.dist
n.tune.params <- length(forest.results[[1]]) - 1

which.best.params.ham <- which(average.hamming.dist.per.param == min(average.hamming.dist.per.param))
which.best.params.n.unique <- which(average.n.unique.elem.per.param == min(average.n.unique.elem.per.param))

best.params.ham <- vector(mode = "list", length = length(which.best.params.ham))
best.params.n.unique <- vector(mode = "list", length = length(which.best.params.n.unique))

#have to run 2 separate loops b/c max indices of each outer loop may not be equal.
for (i in 1:length(which.best.params.ham)){
  params <- vector(mode = "list", length = n.tune.params)
  for (j in 1:n.tune.params){
    #best.params.ham[1,j,i] <- forest.results[[which.best.params.ham[[i]]]][[j]]
    params[[j]] <- forest.results[[which.best.params.ham[[i]]]][[j]]
  }
  best.params.ham[[i]] <- params
}

for (i in 1:length(which.best.params.n.unique)){
  params <- vector(mode = "list", length = n.tune.params)
  for (j in 1:n.tune.params){
    #best.params.n.unique[1,j,i] <- forest.results[[which.best.params.n.unique[[i]]]][[j]]
    params[[j]] <- forest.results[[which.best.params.n.unique[[i]]]][[j]]
  }
  best.params.n.unique[[i]] <- params
}

##concatenate best parameters into a single list and save it as an RDS
best.params <- list("best.params.ham"=best.params.ham,
                    "best.params.n.unique"=best.params.n.unique)
saveRDS(best.params, "./res/BestParamsQualified.RDS")

#now tune refined profile for classification accuracy.

#get refined profile metabolites from the best qualified tune result based on best.params.n.unique

##find which parameter set in which.best.params.n.unique
which.best.params <- which.best.params.n.unique

##save best forest
saveRDS(forest.results[[which.best.params[1]]][["forests"]][[which(oob.error.all[,,which.best.params[1]] == min(oob.error.all[,,which.best.params[1]]))[1]]][["randomForest.object"]], "./res/BestForestQualified.RDS")

##pick which qualified profile's forest to use within a given parameter set based on lowest oob.error
##Ties are broken by taking the first (by index).
ref.prof <- ref.prof.all[,
                         which(oob.error.all[,,which.best.params[1]] == min(oob.error.all[,,which.best.params[1]]))[1],
                         which.best.params[1]]
ref.prof.mda <- ref.prof.mda.all[,
                                 which(oob.error.all[,,which.best.params[1]] == min(oob.error.all[,,which.best.params[1]]))[1],
                                 which.best.params[1]]
ref.prof.df <- data.frame(ref.prof, ref.prof.mda)

##update met_info with refined profile info and save met_info
met_info <- merge(met_info, ref.prof.df, by.x = "all_mets", by.y = "ref.prof", all.x = TRUE)
write.csv(met_info, "./res/met_info.csv")

##get list of variables in refined profile by taking the top n=ref.prof.size variables by mda from qualified rf run
ref.prof.classifier <- ref.prof[1:ref.prof.size]

##save list of metabolites in refined profile
write.csv(ref.prof.classifier, "./res/list.ref.prof.metabolites.csv", row.names = FALSE)

#as well as a formatted version to use in identification
formatted.ref.prof.classifier <- fix.peak.labels(ref.prof.classifier)
write.csv(formatted.ref.prof.classifier, "./res/formatted.list.ref.prof.metabolites.csv", row.names = FALSE)

#make data frame of refined profile and save it
dat.ref <- dat.qual[,c(class.label, ref.prof.classifier)]
dat.ref[,class.label] <- as.numeric(as.factor(dat.ref[,class.label])) - 1
write.csv(dat.ref, "./res/dat.ref.csv", row.names = FALSE)

##create ref.formula
ref.formula <- as.formula(paste("as.factor(Case_Group) ~", paste(ref.prof.classifier, collapse = "+"), sep = ""))

##set refined profile tuning run parameters
set.seed(42)
ntrial.ref <- 30
ntrial.per.param.ref <- 10
inc.ntree.ref <- 500
forest.results.ref <- vector(mode = "list", length = ntrial.ref)

##set up parallel backend and register it with foreach package
cluster <- makePSOCKcluster(where.to.make.cluster)
registerDoParallel(cluster)

##grow all random forests using random hyperparameter search method
forest.results.ref <- foreach(i=1:ntrial.ref, .combine = 'c', .multicombine = TRUE, .packages = "randomForest", .options.rng = 999) %dorng% {
  forest.result <- list("mtry" = floor(sqrt(ref.prof.size)),
                        "ntree" = i*inc.ntree.ref,
                        "nodesize" = 1,
                        "sampsize" = c(-999, -999),
                        "forests" = vector(mode = "list", length = ntrial.per.param.ref))

  for (j in 1:ntrial.per.param.ref) {
    forest <- randomForest(ref.formula, data = dat.ref,
                           ntree = forest.result[["ntree"]],
                           mtry = forest.result[["mtry"]],
                           #sampsize = forest.result[["sampsize"]],
                           nodesize = forest.result[["nodesize"]],
                           importance = TRUE, keep.forest = TRUE, keep.inbag = TRUE)
    forest.ref.prof <- rownames(forest$importance[order(forest$importance[,"MeanDecreaseAccuracy"], decreasing = TRUE),])
    forest.ref.prof.mda <- forest$importance[order(forest$importance[,"MeanDecreaseAccuracy"], decreasing = TRUE),"MeanDecreaseAccuracy"]
    forest.oob.error <- forest$confusion
    forest.result[["forests"]][[j]] <- list("forest.ref.prof" = forest.ref.prof,
                                            "forest.ref.prof.mda" = forest.ref.prof.mda,
                                            "forest.oob.error" = forest.oob.error,
                                            "randomForest.object" = forest)
  }
  return(forest.result)
}
end.time.r <- Sys.time()
stopCluster(cluster)

##save refined profile results
saveRDS(forest.results.ref, "./res/RandomHyperparamSearchRefined.RDS")

##get forest.results.ref into nested structure that is compatible with code below
name.ref <- c("mtry", "ntree", "nodesize", "sampsize", "forests")
forest.results.2.ref <- vector(mode = "list", length = ntrial.ref)
for (i in 0:(ntrial.ref - 1)){
  forest <- vector(mode="list", length = 0)
  which.elems <- 1 + (5*i):(4+5*i)
  for (j in 1:5){
    forest[[name.ref[[j]]]] <- forest.results.ref[[which.elems[[j]]]]
  }
  forest.results.2.ref[[i+1]] <- forest
}

forest.results.ref <- forest.results.2.ref

##make array of all oob errors and find the smallest mean oob error

##get oob error from forest.results.ref (all the rf results) and put it in an array
oob.error.all.ref <- array(dim = c(1,ntrial.per.param.ref, ntrial.ref))

for (i in 1:ntrial.ref) {
  for (j in 1:ntrial.per.param.ref){
    oob.error.all.ref[,j,i] <- mean(forest.results.ref[[i]][["forests"]][[j]][["forest.oob.error"]][,"class.error"])
  }
}


##compute mean oob error rate for each set of parameters
oob.array.ref <- array(dim = c(2, ntrial.ref))

for (i in 1:ntrial.ref) {

  oob.array.ref[1,i] <- mean(oob.error.all.ref[1,,i])
  oob.array.ref[2,i] <- sd(oob.error.all.ref[1,,i])
}

#save results
saveRDS(oob.array.ref, "./res/oob.array.ref.RDS")

#save results in easy-to-plot format
tune.results.refined <- data.frame(t(oob.array.ref))
tune.results.refined$ntree <- seq(inc.ntree.ref, inc.ntree.ref*ntrial.ref, by = inc.ntree.ref)
names(tune.results.refined) <- c("mean", "sd", "ntree")
write.csv(tune.results.refined, "./res/tune.results.refined.csv", row.names = FALSE)

#plot
postscript("./res/ref.prof.oob.error.plot.eps")
plot(seq(inc.ntree.ref, inc.ntree.ref*ntrial.ref, by = inc.ntree.ref), oob.array.ref[1,])
arrows(oob.array.ref[1,], oob.array.ref[1,] - oob.array.ref[2,], oob.array.ref[1,] + oob.array.ref[2,])
dev.off()


#lowest number of trees to achieve minimum observed oob error
min.mean.oob.error <- min(oob.array.ref[1,])
which.min.oob <- which(oob.array.ref[1,] == min.mean.oob.error)[1]

#get importance for first forest in parameter set that achieved minimal oob
importance.ref <- forest.results.ref[[which.min.oob]][["forests"]][[1]][["forest.ref.prof.mda"]]

#calculate cpi for this forest
library(permimp)
best.forest.ref.cpi <- permimp(forest.results.ref[[which.min.oob]][["forests"]][[1]][["randomForest.object"]],
                          nperm = 1, OOB = TRUE, scaled = FALSE,
                          conditional = TRUE, threshold = 0, thresholdDiagnostics = FALSE,
                          progressBar = FALSE, do_check = FALSE)

#merge importance results and save
ref.prof.importances <- merge(data.frame(names(importance.ref), importance.ref),
                              data.frame(names(best.forest.ref.cpi$values), best.forest.ref.cpi$values),
                              by.x = "names.importance.ref.", by.y = "names.best.forest.ref.cpi.values.")

write.csv(ref.prof.importances, "./res/ref.prof.importances.csv", row.names = FALSE)

end.time.total <- Sys.time()
times <- list(start.time.total, start.time.q, end.time.q, start.time.r, end.time.r, end.time.total)
write.csv(times, "./res/times.csv")