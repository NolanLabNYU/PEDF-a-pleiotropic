#compile variable information and results, save for graphing, calculate correlation for variable importance rankings

#load libraries
library(haven)
library(Hmisc)

##compile met_info
#load results from tune process
met_info <- read.csv("./res/met_info.csv", header = TRUE, stringsAsFactors = FALSE)

#load variable metadata
met_info_metadata <- read_sav("./dat/Variable Metadata.sav")

#which columns to take from metadata and take them
which.cols.metadata <- c("Metabolite_Name_R", "Name_For_Graph", "Subpathway_Graph", "Pathway_Graph")
met_info_metadata <- met_info_metadata[,which.cols.metadata]

#merge metadata with classical pi scores
met_info <- merge(met_info, met_info_metadata, by.x = "all_mets", by.y = "Metabolite_Name_R")

#load cpi scores
best.forest.qualified.cpi <- readRDS("./res/best.forest.qualified.cpi.RDS")
best.forest.qualified.cpi <- data.frame(names(best.forest.qualified.cpi$values), best.forest.qualified.cpi$values)

#merge cpi scores onto metadata
met_info <- merge(met_info, best.forest.qualified.cpi, by.x = "all_mets", by.y = "names.best.forest.qualified.cpi.values.")


#make formatted compound name for graph
##make proper capitalization
met_info$compound_graph <- paste0(toupper(substr(met_info[,"Name_For_Graph"], 1, 1)),
                                          tolower(substr(met_info[,"Name_For_Graph"], 2, nchar(met_info[,"Name_For_Graph"])))
)

#which are refined profile based on mpi
met_info$in.refined.profile.mpi <- met_info$all_mets %in% met_info[order(met_info$ref.prof.mda, decreasing = TRUE)[1:35], "all_mets"]

#which are refined profile based on cpi
met_info$in.refined.profile.cpi <- met_info$all_mets %in% met_info[order(met_info$best.forest.qualified.cpi.values, decreasing = TRUE)[1:35], "all_mets"]

#add results from refined profile
ref.prof.importances <- read.csv("./res/ref.prof.importances.csv", stringsAsFactors = FALSE, header = TRUE)
met_info <- merge(met_info, ref.prof.importances, by.x = "all_mets", by.y = "names.importance.ref.")

#assess similarity of variable importance rankings
variable.rank.correlations <- rcorr(as.matrix(met_info[,c("ref.prof.mda", "importance.ref", "best.forest.ref.cpi.values")]), type = "spearman")

#save to sav for graphing
write_sav(met_info, "./res/met_info.sav")
