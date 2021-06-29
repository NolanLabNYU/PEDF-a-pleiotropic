##cpi.R measure cpi using custom parallel permimp on hpc cluster

#load best forest from rf_qual.R
best.forest.qualified <- readRDS("./res/BestForestQualified.RDS")

#load data
dat.qual <- read.csv("./dat/RF1 Data.csv", stringsAsFactors = FALSE)

#define class label variable name used in dat.qual
class.label <- "Case_Group"

#create formula (rf.formula is what the formula is called in the best.forest.qualified randomForest object)
rf.formula <- as.formula(paste("as.factor(Case_Group) ~", paste(setdiff(names(dat.qual), class.label), collapse = "+"), sep = ""))

#calculate cpi using custom parallel cpi package
library(permimp)
best.forest.qualified.cpi <- permimp(best.forest.qualified,
                   nperm = 1, OOB = TRUE, scaled = FALSE,
                   conditional = TRUE, threshold = 0, thresholdDiagnostics = FALSE,
                   progressBar = FALSE, do_check = FALSE)

#save cpi
saveRDS(best.forest.qualified.cpi, "./res/best.forest.qualified.cpi.RDS")