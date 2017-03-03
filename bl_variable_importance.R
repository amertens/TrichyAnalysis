#Baseline variable importance


library("devtools")
library("varImpact")
library(doSNOW)
library(RhpcBLASctl)

load("prescreened_blcov.Rdata")

dim(df)
Y<-df$diar7d
X<-W



# Customize Q and g libraries for TMLE estimation.
Q_lib <- c("SL.gam","SL.glmnet", "SL.stepAIC", "SL.randomForest", "SL.rpartPrune", "SL.bayesglm")
g_lib <- c("SL.stepAIC", "SL.glmnet")

# doSNOW parallel
library(doSNOW)
library(RhpcBLASctl)
# Detect the number of physical cores on this computer using RhpcBLASctl.
cluster <- makeCluster(get_num_cores())

registerDoSNOW(cluster)
vim <- varImpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib)
vim
stopCluster(cluster)
