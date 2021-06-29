##regressions
set.seed(123)
library(haven, quietly = TRUE)
library(pROC, quietly = TRUE)
library(boot, quietly = TRUE)
library(Hmisc, quietly = TRUE)
library(glmnet, quietly = TRUE)
source("../rfuncs/clean.summary.R")

##dat.metx contains cases in the metabolomics cohort
##dat.wc contains cases in the validation cohort

##variables medians and iqr
##for WTC-LI
dat.wc.worst <- dat.wc[dat.wc$W_C == 1,]
clean.summary(summary(dat.wc.worst, digits=15, quantile.type=6), warn=FALSE)
clean.summary(summary(dat.metx[dat.metx$W_C == 1,], digits=15, quantile.type=6), warn=FALSE)

#compute n(%) for exposure and race
n.worst <- nrow(dat.wc.worst)

#Exposure
exposures <- list("High", "Intermediate", "Low")
for (i in 1:3) {
  n <- sum(dat.wc.worst$Exposure_group_123.factor == i)
  print(paste(exposures[[i]],":", n,"(",n/n.worst,")"))
}

#Race
races <- list("Caucasian", "African American")
for (i in 1:2) {
  n <- sum(dat.wc.worst$RACE.factor == i)
  print(paste(races[[i]],":", n,"(",n/n.worst,")"))
}

##for Controls
dat.wc.control <- dat.wc[dat.wc$W_C == 0,]
clean.summary(summary(dat.wc.control, digits=15, quantile.type=6), warn=FALSE)
clean.summary(summary(dat.metx[dat.metx$W_C == 0,], digits=15, quantile.type=6), warn=FALSE)

#compute n(%) for exposure and race
n.control <- nrow(dat.wc.control)

#Exposure
for (i in 1:3) {
  n <- sum(dat.wc.control$Exposure_group_123.factor == i)
  print(paste(exposures[[i]],":", n, "(", n/n.control,")"))
}

#Race
for (i in 1:2) {
  n <- sum(dat.wc.control$RACE.factor == i)
  print(paste(races[[i]],":", n,"(",n/n.control,")"))
}

##significance testing
vars <- unique(colnames(dat.wc))
dat.wc <- dat.wc[,vars]
dat.metx <- dat.metx[,vars]
results.sig.wc.metx <- matrix(ncol = 2, nrow = length(vars))


##to compare W to C for each cohort (significance results for intracohort (ic) comparison)
dat <- list(dat.wc, dat.metx)
results.sig.ic <- matrix(nrow = length(vars), ncol = 3)
for (k in 1:3){
  for (i in 1:length(vars)){
  which.W <- dat[[k]]$W_C == 1
  which.C <- dat[[k]]$W_C == 0
  if (!is.numeric(unlist(dat[[k]][,i]))) {
    results.sig.ic[i,k] <- chisq.test(table(dat[[k]][c("W_C", paste(vars[[i]]))]), correct = TRUE)$p.value
  }
  else{
    results.sig.ic[i,k] <- wilcox.test(x = unlist(dat[[k]][which.C,i]), y = unlist(dat[[k]][which.W,i]))$p.value
  }
  }
}

results.sig.ic <- data.frame(results.sig.ic)
colnames(results.sig.ic) <- c("wvc.wc", "wvc.metx")
rownames(results.sig.ic) <- vars
results.sig.ic
for (i in 0:3){
  print(paste(i, "#################"))
  print(results.sig.ic[rowSums(results.sig.ic <= 0.05, na.rm = TRUE) == i,])
}

##wc.metx
for (i in 1:length(vars)){
  for (j in 0:1){
    if (!is.numeric(unlist(dat.wc[,i]))) {
      rows.wc <- cbind(rep("child", times = sum(dat.wc$W_C==j)), dat.wc[dat.wc$W_C==j,i])
      rows.metx <- cbind(rep("parent", times = sum(dat.metx$W_C==j)), dat.metx[dat.metx$W_C==j,i])
      colnames(rows.wc)[1] <- colnames(rows.metx)[1]
      results.sig.wc.metx[i,j+1] <- chisq.test(table(rbind(rows.wc, rows.metx)), correct = TRUE)$p.value
    }
    else {
      results.sig.wc.metx[i,j+1] <- wilcox.test(x = unlist(dat.wc[dat.wc$W_C==j,i]), y = unlist(dat.metx[dat.metx$W_C==j,i]))$p.value
    }
  }
}
results.sig.wc.metx <- data.frame(results.sig.wc.metx)
rownames(results.sig.wc.metx) <- vars
colnames(results.sig.wc.metx) <- c("C vs. C", "W vs. W")
results.sig.wc.metx
results.sig.wc.metx[rowSums(results.sig.wc.metx <= 0.05) >= 1,]

##also print sizes of cohorts
cohort.sizes <- matrix(nrow = 3, ncol = 2)

for (i in 1:length(dat)){
  cohort.sizes[,i] <- c(sum(dat[[i]]$W_C==1), sum(dat[[i]]$W_C==0), nrow(dat[[i]]))
}
cohort.sizes <- as.data.frame(cohort.sizes)
cohort.sizes.colnames <- c("dat.wc", "dat.metx")
cohort.sizes.rownames <- c("w.size", "c.size", "total.size")
colnames(cohort.sizes) <- cohort.sizes.colnames
rownames(cohort.sizes) <- cohort.sizes.rownames

cohort.sizes

#LASSO regression
x.matrix <- dat.wc[,c("PEDF", "MDC", "BP_SYS", "MIP_4", "GRO",
                      "MCP_1", "sIL_2Ra", "Amylin",
                      "sCD40L", "MMP1", "sVEGFR1", "EGF",
                      "Leptin", "Apo_AII", "MMP13", "BMI_WTCHP",
                      "age_on_911", "Exposure_group_123.factor", "BMI_SPE", "pre_911_fev1_pct_pred"
)]
x.matrix <- model.matrix(~.-1, x.matrix)
lasso.cv <- cv.glmnet(x = x.matrix,
                   y =  as.matrix(dat.wc[,c("W_C")]),
                   family = "binomial", type.measure = "auc", nfolds = 5)

#model coefficients
coef(lasso.cv, s=lasso.cv$lambda.min)

#model AUC
which.highest.auc <- which(lasso.cv$cvm == max(lasso.cv$cvm))
lasso.cv$lambda[which.highest.auc] == lasso.cv$lambda.min
lasso.cv$cvm[which.highest.auc]
lasso.cv$cvsd[which.highest.auc]

#save lasso plot
setEPS()
postscript("./fig/lasso_tune.eps")
plot(lasso.cv)
dev.off()

#save lasso results
lasso.results <- data.frame(log(lasso.cv$lambda), lasso.cv$cvm, lasso.cv$cvsd)
write.csv(lasso.results,"./res/lasso.results.csv", row.names = FALSE)

#make correlation matrix to identify significantly correlated variables
cormat <- rcorr(as.matrix(dat.wc[,c("PEDF", "MDC", "BP_SYS",
                                    "MIP_4", "GRO", "MCP_1")]), 
                type = "pearson")
which.sig.p.less.05 <- cormat$P < 0.05

##fit crude models
###to be reported in table 3
crude.model.variables <- list("PEDF", "MDC", "MIP_4", "GRO", "MCP_1", "BP_SYS",
                              "MMP1", "Apo_AII")
crude.models <- list()
roc.crude.models <- list()
for (i in 1:length(crude.model.variables)) {
crude.formula <- as.formula(paste("W_C ~ ", crude.model.variables[i], "+ age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred", sep = ""))
assign(paste("crude.", crude.model.variables[i], sep = ""), glm(crude.formula,
                  family = binomial,
                  data = dat.wc,
                  na.action = na.exclude,
                  maxit = 100,
                  trace = FALSE))

crude.models[[i]] <- get(paste("crude.", crude.model.variables[i], sep = ""))
df <- data.frame(crude.models[[i]]$y, crude.models[[i]]$linear.predictors)
roc.crude.models[[i]] <- roc(df[,1], df[,2])
}

for (i in 1:length(crude.models)){
  print(summary(crude.models[[i]]))
  print("***Odds ratios and 95% CI***")
  print(exp(crude.models[[i]]$coefficients))
  print(exp(confint(crude.models[[i]], level = 0.95)))
  print("***AUC_ROC and 95%CI")
  print(auc(roc.crude.models[[i]]))
  print(ci(roc.crude.models[[i]], of = "auc"))
  print("/////////////////////////////////////////////////////////////")
  }

#find cutpoints based on Youden's index
#all cutpoints below are based on Youden's index

##these cutpoints will be used to dichotomize variables
##some of these dichotomized variables will be entered in regression models to
##be reported in table 3, and entered simultaneously into the
##final model, reported in table 4.

###some continuous variables are correlated with each other.
###where this occurs (indicated in the code below), these variables 
###are transformed into a composite variable. the composite variable is then
###entered into a regression model without other serum biomarkers, and
###the results of this model are reported in table 3. finally, the composite
###variable is included simultaneously with other variables in the final model
###described in table 4.

##pedf cutpoint
pedf.roc <- roc(dat.wc$W_C, dat.wc$PEDF)
pedf.cp <- coords(pedf.roc, x = "best", best.method = "youden")[[1]]
dat.wc$PEDF.grtr.youden.cp <- dat.wc$PEDF < pedf.cp
pedf.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$PEDF.grtr.youden.cp != (dat.wc$PEDF < 3.94))

## mip4 cutpoint
mip4.roc <- roc(dat.wc$W_C, dat.wc$MIP_4)
mip4.cp <- coords(mip4.roc, x = "best", best.method = "youden")[[1]]
dat.wc$MIP_4.grtr.youden.cp <- dat.wc$MIP_4 < mip4.cp 
mip4.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$MIP_4.grtr.youden.cp != (dat.wc$MIP_4 < 368.70))

##pedf/mip4 combination, and its performance alone in a glm (crude model)
###results reported in table 3
dat.wc$pedf.mip4 <- dat.wc$PEDF.grtr.youden.cp & dat.wc$MIP_4.grtr.youden.cp

pedf.mip4.crude.model <- glm(W_C ~ pedf.mip4
                               + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                       family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

pedf.mip4.crude.model.df <- data.frame(pedf.mip4.crude.model$y, pedf.mip4.crude.model$linear.predictors)
pedf.mip4.crude.model.roc <- roc(pedf.mip4.crude.model.df$pedf.mip4.crude.model.y, pedf.mip4.crude.model.df$pedf.mip4.crude.model.linear.predictors)
auc(pedf.mip4.crude.model.roc)
ci(pedf.mip4.crude.model.roc, of = "auc")
summary(pedf.mip4.crude.model)

##Odds ratios for above model
exp(pedf.mip4.crude.model$coefficients)
exp(confint(pedf.mip4.crude.model, level = 0.95))

##mdc cutpoint (mdc is correlated with mcp1 and gro)
mdc.roc <- roc(dat.wc$W_C, dat.wc$MDC)
mdc.cp <- coords(mdc.roc, x = "best", best.method = "youden")[[1]]
dat.wc$MDC.grtr.youden.cp <- dat.wc$MDC < mdc.cp
mdc.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$MDC.grtr.youden.cp != (dat.wc$MDC < 1756.99))

##mcp1 cutpoint
mcp1.roc <- roc(dat.wc$W_C, dat.wc$MCP_1)
mcp1.cp <- coords(mcp1.roc, x = "best", best.method = "youden")[[1]]
dat.wc$MCP_1.grtr.youden.cp <- dat.wc$MCP_1 > mcp1.cp
mcp1.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$MCP_1.grtr.youden.cp != (dat.wc$MCP_1 > 279.07))

##gro cutpoint
gro.roc <- roc(dat.wc$W_C, dat.wc$GRO)
gro.cp <- coords(gro.roc, x = "best", best.method = "youden")[[1]]
dat.wc$GRO.grtr.youden.cp <- dat.wc$GRO > gro.cp
gro.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$GRO.grtr.youden.cp != (dat.wc$GRO > 580.28))

##mcp1/gro/mdc combination and its performance in a glm (crude model)
dat.wc$gro.mcp1.mdc <- dat.wc$GRO.grtr.youden.cp & dat.wc$MCP_1.grtr.youden.cp & dat.wc$MDC.grtr.youden.cp

gro.mcp1.mdc.crude.model <- glm(W_C ~ gro.mcp1.mdc
                           + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                       family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

gro.mcp1.mdc.crude.model.df <- data.frame(gro.mcp1.mdc.crude.model$y, gro.mcp1.mdc.crude.model$linear.predictors)
gro.mcp1.mdc.crude.model.roc <- roc(gro.mcp1.mdc.crude.model.df$gro.mcp1.mdc.crude.model.y, gro.mcp1.mdc.crude.model.df$gro.mcp1.mdc.crude.model.linear.predictors)
auc(gro.mcp1.mdc.crude.model.roc)
ci(gro.mcp1.mdc.crude.model.roc, of = "auc")
summary(gro.mcp1.mdc.crude.model)

##Odds ratios for above model
exp(gro.mcp1.mdc.crude.model$coefficients)
exp(confint(gro.mcp1.mdc.crude.model, level = 0.95))

##SBP cut point
sbp.roc <- roc(dat.wc$W_C, dat.wc$BP_SYS)
sbp.cp <- coords(sbp.roc, x = "best", best.method = "youden")[[1]]
dat.wc$BP_SYS.grtr.youden.cp <- dat.wc$BP_SYS > sbp.cp
sbp.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$BP_SYS.grtr.youden.cp != (dat.wc$BP_SYS > 127.5))

##sbp cutpoint crude model
sbp.cp.crude.model <- glm(W_C ~ BP_SYS.grtr.youden.cp
                          + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                       family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

sbp.cp.crude.model.df <- data.frame(sbp.cp.crude.model$y, sbp.cp.crude.model$linear.predictors)
sbp.cp.crude.model.roc <- roc(sbp.cp.crude.model.df$sbp.cp.crude.model.y, sbp.cp.crude.model.df$sbp.cp.crude.model.linear.predictors)
auc(sbp.cp.crude.model.roc)
ci(sbp.cp.crude.model.roc, of = "auc")
summary(sbp.cp.crude.model)

##Odds ratios for above model
exp(sbp.cp.crude.model$coefficients)
exp(confint(sbp.cp.crude.model, level = 0.95))

##Apo AII cut point
apo.roc <- roc(dat.wc$W_C, dat.wc$Apo_AII)
apo.cp <- coords(apo.roc, x = "best", best.method = "youden")[[1]]
dat.wc$Apo_AII.grtr.youden.cp <- dat.wc$Apo_AII > apo.cp
apo.cp

#check if rounded cutpoint is equivalent #extra division by 1000 in rounded cp is to make sure it agrees w/ units used in table 2 of paper)
sum(dat.wc$Apo_AII.grtr.youden.cp != (dat.wc$Apo_AII/1000 > 1794.22))

##Apo AII cutpoint crude model
apo.cp.crude.model <- glm(W_C ~ Apo_AII.grtr.youden.cp
                          + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                       family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

apo.cp.crude.model.df <- data.frame(apo.cp.crude.model$y, apo.cp.crude.model$linear.predictors)
apo.cp.crude.model.roc <- roc(apo.cp.crude.model.df$apo.cp.crude.model.y, apo.cp.crude.model.df$apo.cp.crude.model.linear.predictors)
auc(apo.cp.crude.model.roc)
ci(apo.cp.crude.model.roc, of = "auc")
summary(apo.cp.crude.model)

##Odds ratios for above model
exp(apo.cp.crude.model$coefficients)
exp(confint(apo.cp.crude.model, level = 0.95))


##mmp1 cut point
mmp1.roc <- roc(dat.wc$W_C, dat.wc$MMP1)
mmp1.cp <- coords(mmp1.roc, x = "best", best.method = "youden")[[1]]
dat.wc$MMP1.grtr.youden.cp <- dat.wc$MMP1 > mmp1.cp
mmp1.cp

#check if rounded cutpoint is equivalent
sum(dat.wc$MMP1.grtr.youden.cp != (dat.wc$MMP1 > 832.48))

##MMP-1 cutpoint crude model
mmp1.cp.crude.model <- glm(W_C ~ MMP1.grtr.youden.cp
                          + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                       family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

mmp1.cp.crude.model.df <- data.frame(mmp1.cp.crude.model$y, mmp1.cp.crude.model$linear.predictors)
mmp1.cp.crude.model.roc <- roc(mmp1.cp.crude.model.df$mmp1.cp.crude.model.y, mmp1.cp.crude.model.df$mmp1.cp.crude.model.linear.predictors)
auc(mmp1.cp.crude.model.roc)
ci(mmp1.cp.crude.model.roc, of = "auc")
summary(mmp1.cp.crude.model)

##Odds ratios for above model
exp(mmp1.cp.crude.model$coefficients)
exp(confint(mmp1.cp.crude.model, level = 0.95))


#now fit the final model
final.model <- glm(W_C ~ pedf.mip4 + MMP1.grtr.youden.cp
                   + Apo_AII.grtr.youden.cp + gro.mcp1.mdc + BP_SYS.grtr.youden.cp
                   + age_on_911 + Exposure_group_123.factor + BMI_SPE + pre_911_fev1_pct_pred,
                      family = "binomial", data = dat.wc, na.action = na.exclude, maxit = 100, trace = FALSE)

final.model.df <- data.frame(final.model$y, final.model$linear.predictors)
final.model.roc <- roc(final.model.df$final.model.y, final.model.df$final.model.linear.predictors)

summary(final.model)

##Odds ratio for above model
exp(final.model$coefficients)
exp(confint(final.model, level = 0.95))

##5-fold cross validation
cv.final.model <- cv.glm(dat.wc, final.model, K = 5)
cv.final.model$delta

##get ROC coordinates for plotting curves for all variables in the final model
m <- length(roc.crude.models)
cols <- list("pedf_mip4_cp_1_spec" = pedf.mip4.crude.model.roc$specificities, "pedf_mip4_cp_sense" = pedf.mip4.crude.model.roc$sensitivities,
                                "mmp1_1-spec" = mmp1.cp.crude.model.roc$specificities, "mmp1_sens" = mmp1.cp.crude.model.roc$sensitivities,
                                "apoaII_1-spec" = apo.cp.crude.model.roc$specificities, "apoaII_sens" = apo.cp.crude.model.roc$sensitivities,
                                "gro_mcp1_mdc-spec" = gro.mcp1.mdc.crude.model.roc$specificities, "gro_mcp1_mdc_sens" = gro.mcp1.mdc.crude.model.roc$sensitivities,
                                "sbp_cp_1-spec" = sbp.cp.crude.model.roc$specificities, "sbp_cp_sens" = sbp.cp.crude.model.roc$sensitivities,
                                "final_1-spec" = final.model.roc$specificities, "final_sens" = final.model.roc$sensitivities
                                )


pad.na <- function(x, max.col.length) {
  ##pads x with NA until length of x is max.col.length
  if(length(x) < max.col.length) {
    x[(length(x) + 1):max.col.length] <- NA
  }
  return(x)
}

max.col.length <- max(unlist(lapply(cols, length)))

##pad columns in cols with NA
for (i in 1:length(cols)){
  cols[[i]] <- pad.na(cols[[i]], max.col.length)
}

roc.coords.df <- data.frame(cols)

roc.coords.df[,seq(from = 1, to = ncol(roc.coords.df), by = 2)] <- 1 - roc.coords.df[,seq(from = 1, to = ncol(roc.coords.df), by = 2)]
write.csv(roc.coords.df, file = "roc_coords.csv")

plot(final.model.roc)

auc(final.model.roc)
ci(final.model.roc, of = "auc")
