# Schizophrenia-linked genetic influences on cognition
# Authors: AMC & ABZ, May 2017

# control analyses for revision
# residuals 

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "pROC", "permute", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "psych", "e1071", "dplyr", "nlme", 
          "psych","stringr", "r2glmm")
# libs <- c("caret", "permute", "randomForest", "doMC")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/ML_genetics/Schizophrenia-Zheutlin")
workDir <- getwd()
dataDir <- paste0(workDir, "/dfs/") 

# # Load custom functions from library
# source("/data/ML_genetics/functions/functions.R")

### Read in pre-processed data 
sweden_snp  <- read.table(paste0(dataDir,"GWASinput2-swe.txt"), header=TRUE, as.is=TRUE)
cnp         <- read.table(paste0(dataDir,"GWASinput2-cnp.txt"), header=TRUE, as.is=TRUE)

# exclude symbol / spatial span
sweden_snp <- sweden_snp[,c(1:11,13:92)]
cnp <- cnp[,c(1:8,10:89)]

# Names of snps and outcomes
snps       <- names(cnp)[12:88] 
outcomes   <- names(cnp)[5:11]

# Rename swedish outcomes to match cnp
names(sweden_snp)[8:14] <- outcomes 
# Some proxy snps used, need to rename
names(sweden_snp)[71:73] <- c("rs2535627", "rs4388249", "rs2905426")

# Remove subjects with missing outcomes in CNP
cnp.df     <- cnp[(complete.cases(cnp$CRT_TIME1)),]
swe.df     <- sweden_snp 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(snps)) %>% as.matrix() 

cnp.pc    <- read.table("/data/swe_gwas/ABZ/ML/cnp_MLsubs_pca.eigenvec", stringsAsFactors = F)
cnp.df.pc <- merge(cnp.df, cnp.pc, by.x = "ptid", by.y = "V1")

residuals <- matrix(nrow=739, ncol=length(outcomes))
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~ age + gender + V3 + V4 + V5 + V6 + V7"))
  m1 <- lm(f1, data = cnp.df.pc, na.action = na.omit)
  residuals[,i] <- residuals(m1)
}

residuals        <- as.data.frame(residuals)
outcomes.res     <- NULL
for (i in 1:7){outcomes.res[i] <- paste0(outcomes[i],".res")}  
names(residuals) <- outcomes.res

targets          <- residuals %>% apply(2,as.numeric) %>% as.data.frame()

#  apply scaling
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)
targets  <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame




### Machine Learning preparation
## Set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE,
                           selectionFunction = "oneSE")

## Hyperparameter grid search
rfGrid     <- expand.grid(.mtry = c(3,4,10)) 


##### ##### ##### ##### ##### 
##### Model training
##### ##### ##### ##### ##### 

glmPipeline <- function(response, Xmat = predictors, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "glm",
                 metric = "Rsquared",
                 trControl = cvpar)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}

glm.loop     <- lapply(targets, function(i) glmPipeline(i))
glm.models   <- lapply(glm.loop, function(i) i[[1]])
glm.perfs    <- lapply(glm.loop, function(i) i[[2]])


rfPipeline <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance = TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop     <- lapply(targets, function(i) rfPipeline(i))
rf.models   <- lapply(rf.loop, function(i) i[[1]])
rf.perfs    <- lapply(rf.loop, function(i) i[[2]])



#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),] #rm NAs

rep.outcomes   <- raw_rep %>% dplyr::select(one_of(outcomes))

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW   <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW    <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398

PCs    <- read.table("/data/swe_gwas/ABZ/ML/swe_MLsubs2.mds", header=T, stringsAsFactors = F)
swe.pc <- merge(raw_rep, PCs, by = "IID")

swe.pc.use <- swe.pc[complete.cases(swe.pc[,6:14]),]

resid.swe <- matrix(nrow=339, ncol=length(outcomes))
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~ age + sex + C1 + C2 + C3 + C4 + C5"))
  m1 <- lm(f1, data = swe.pc.use, na.action = na.omit)
  resid.swe[,i] <- residuals(m1)
}

resid.swe        <- as.data.frame(resid.swe)
names(resid.swe) <- outcomes.res

rep.outcomes.resid <- resid.swe %>% apply(2,as.numeric) %>% as.data.frame()


# Apply scaling transformations from CNP
rep.outcomes.s   <- (rep.outcomes.resid - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))

rep.predictors   <- swe.pc.use %>% dplyr::select(one_of(snps)) %>% as.matrix

rf.pred.tr       <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr      <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))


##### ##### ##### ##### ##### 
##### Mixed effects models controlling for family ID
##### ##### ##### ##### ##### 

rf.pred  <- rf.pred.tr %>% as.data.frame
glm.pred <- glm.pred.tr %>% as.data.frame
mix      <- cbind(rep.outcomes.s, rf.pred, glm.pred, swe.pc.use[,"FID.x"]) %>% as.data.frame

outcomes.rf <- NULL  # start empty list
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7]   <- outcomes.res            # renaming using the list of names
names(mix)[8:14]  <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22]    <- "FID"

outstats <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes.res[i],"~",outcomes[i],".rf"))
  f2 <- formula(paste0(outcomes.res[i],"~",outcomes[i],".glm"))
  m1 <- lme(fixed=f1, random=~1|FID,data=mix,na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|FID,data=mix,na.action=na.omit)
  outstats[i,1] <- outcomes[i]
  outstats[i,2] <- summary(m1)$tTable[2,4]
  outstats[i,3] <- summary(m1)$tTable[2,5]
  outstats[i,4] <- summary(m2)$tTable[2,4]
  outstats[i,5] <- summary(m2)$tTable[2,5]
}

outstats <- as.data.frame(outstats)
names(outstats) <- c("var","tval.rf","pval.rf","tval.glm","pval.glm")

write.table(outstats, "revision1/originalresults_pcs.txt", col.names=T, row.names=F, sep="\t", quote=F)

save.image("../results5_amended.RData")
