## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", 
                      fig.width = 7, fig.height = 5)

## ----eval = TRUE, echo = TRUE-------------------------------------------------
library("party", quietly = TRUE)
library("permimp")
set.seed(542863)
airq <- subset(airquality, !(is.na(Ozone) | is.na(Solar.R)))
cfAirq50 <- cforest(Ozone ~ ., data = airq,
                    control = cforest_unbiased(mtry = 2, ntree = 50,
                                              minbucket = 5, 
                                              minsplit = 10))

## ----eval = TRUE, echo = TRUE-------------------------------------------------
system.time(CPI_permimp <- permimp(cfAirq50, conditional = TRUE, progressBar = FALSE))
system.time(CPI_varimp <- varimp(cfAirq50, conditional = TRUE))
CPI_permimp
CPI_varimp

## ----eval = TRUE, echo = TRUE-------------------------------------------------
CPI_permimp <- permimp(cfAirq50, conditional = TRUE, threshold = .2, progressBar = FALSE)
CPI_permimp
CPI_varimp

## ----eval = TRUE, echo = TRUE-------------------------------------------------
CPI_varimp <- varimp(cfAirq50, conditional = TRUE, threshold = .95)
CPI_permimp
CPI_varimp

## ----eval = TRUE, echo = TRUE-------------------------------------------------
set.seed(542863)
system.time(CPI_asParty <- permimp(cfAirq50, conditional = TRUE, asParty = TRUE, progressBar = FALSE))
set.seed(542863)
system.time(CPI_varimp <- varimp(cfAirq50, conditional = TRUE))
CPI_asParty
CPI_varimp

## ----eval = TRUE, echo = TRUE-------------------------------------------------
## varimp returns a named numerical vector.
str(CPI_varimp)

## permimp returns a VarImp-object.
str(CPI_asParty)

## the results of permimp(asParty = TRUE) and varimp() are exactly the same.
all(CPI_asParty$values == CPI_varimp)

## ----eval = TRUE, echo = TRUE-------------------------------------------------
## Original Unconditional Permutation Importance
set.seed(542863)
PI_permimp <- permimp(cfAirq50, progressBar = FALSE, pre1.0_0 = TRUE)
set.seed(542863)
PI_varimp <- varimp(cfAirq50, pre1.0_0 = TRUE)
PI_permimp
PI_varimp

## Splitwise Unconditional Permutation Importance
set.seed(542863)
PI_permimp2 <- permimp(cfAirq50, progressBar = FALSE)
set.seed(542863)
PI_varimp2 <- varimp(cfAirq50)
PI_permimp2
PI_varimp2

## ----eval = TRUE, echo = TRUE-------------------------------------------------
## fit a new forest with 500 trees
set.seed(542863)
cfAirq500 <- cforest(Ozone ~ ., data = airq,
                     control = cforest_unbiased(mtry = 2, ntree = 500,
                                              minbucket = 5, 
                                              minsplit = 10))

## compute permutation importance
PI_permimp500 <- permimp(cfAirq500, progressBar = FALSE)

## different plots, all easy to make
## barplot
plot(PI_permimp500, type = "bar")

## barplot with visualization of the distribution: an
## interval between the .25 and .75 quantiles of the per 
## Tree values is added to the plot
plot(PI_permimp500, type = "bar", interval = "quantile")

## horizontal boxplot
plot(PI_permimp500, type = "box", horizontal = TRUE)

## unsorted dotplot
plot(PI_permimp500, type = "dot", sort = FALSE, 
     interval = "quantile")

## ----eval = TRUE, echo = TRUE-------------------------------------------------
library("randomForest", quietly = TRUE)
set.seed(542863)
rfAirq50 <- randomForest(Ozone ~ ., data = airq, mtry = 2, replace = FALSE, 
                         nodesize = 7, keep.forest = TRUE, keep.inbag = TRUE)

## ---- eval = TRUE, echo = FALSE-----------------------------------------------
CPI_permimpRF <- permimp(rfAirq50, conditional = TRUE, progressBar = FALSE, do_check = FALSE)
plot(CPI_permimpRF, horizontal = TRUE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  CPI_permimpRF <- permimp(rfAirq50, conditional = TRUE, progressBar = FALSE)
#  plot(CPI_permimpRF, horizontal = TRUE)

