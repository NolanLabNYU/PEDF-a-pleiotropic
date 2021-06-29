# PEDF-a-pleiotropic
catalogue of files in repo

### Custom Functions
==> clean.summary.R <==  
clean.summary <- function(arg1, warn=TRUE){
  ##takes output of summary() and returns an html table with median(iqr) for continuous (numeric) variables

  ==> fix.peak.labels.R <==  
fix.peak.labels <- function(arg1) {
  #takes a vector of strings of metabolite/peak names (that were altered by R's makenames())
  #and fixes names to match the format in data sheets
  
==> ham.dist.R <==  
ham.dist <- function(arg1) {
  ##function ham.dist computes all pairwise hamming distances among the columns of arg1, a matrix

==> parseSlurmNodeList.R <==  
parseSlurmNodeList<- function(arg1) {
  #parses SLURM_NODELIST environment variable into a list of nodes (i.e. to be passed to makePSOCKcluster() or h2o.init())

==> uniqueEl.R <==  
uniqueEl <- function(arg1) {
  ##computes the number of unique elements in a pair of sets, for all pairs of columns in a matrix

==> permimp <==  
customized version of permimp to support parallel processing on SLURM
modifications were made to /R/doPermimp.R, DESCRIPTION, NAMESPACE, and MD5
to build, use: R CMD build permimp --no-build-vignettes

### Analysis Scripts
==> compile_met_info_and_results.R <==  
compile variable information and results, save for graphing, calculate correlation for variable importance rankings

==> cpi.R <==  
cpi.R measure cpi using custom parallel permimp on hpc cluster

==> plot_tune_results.R <==  
graph tune results

==> regressions.R <==  
regressions for validation cohort

==> rf_tune.R <==  
rf hyperparameter search/tune

==> clustergram_code.m <==  
code to reproduce clustergram

### Data Files
==> roc_coords.csv <==  
contains coordinates for ROC curves

==> met_info.sav <==  
contains classical permutation importance and conditional permutation importance as described in paper
can be read into R using haven::read_sav("met_info.sav")

==> tune.results.csv <==  
contains results of tuning process for random forests with all variables

==> tune.results.refined.csv <==  
contains results of tuning process for random forests with refined variable profile

==> lasso.results.csv <==  
contains results of lasso tuning procedure
