# Schizophrenia-linked genetic influences on cognition
# Authors: AMC & ABZ, May 2017

# control analyses for revision
# add covariates in training 

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "pROC", "permute", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "psych", "e1071", "dplyr", "nlme", 
          "psych","stringr", "r2glmm")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/ML_genetics/Schizophrenia-Zheutlin")
workDir <- getwd()
dataDir <- paste0(workDir, "/dfs/") 

### Read in pre-processed data 
sweden_snp  <- read.table(paste0(dataDir,"GWASinput2-swe.txt"), header=TRUE, as.is=TRUE)
cnp         <- read.table(paste0(dataDir,"GWASinput2-cnp.txt"), header=TRUE, as.is=TRUE)

# exclude symbol / spatial span
sweden_snp <- sweden_snp[,c(1:11,13:92)]
cnp <- cnp[,c(1:8,10:89)]

# add PCs
cnp.pc    <- read.table("/data/swe_gwas/ABZ/ML/cnp_MLsubs_pca.eigenvec", stringsAsFactors = F)
cnp.df.pc <- merge(cnp, cnp.pc, by.x = "ptid", by.y = "V1")

pc.names <- NULL
for (i in 1:20){pc.names[i] <- paste0("PC",i)}  
names(cnp.df.pc)[89:109] <- c("ptid2",pc.names)

swe.pc    <- read.table("/data/swe_gwas/ABZ/ML/swe_MLsubs2.mds", header=T, stringsAsFactors = F)
swe.df.pc <- merge(sweden_snp, swe.pc, by = "IID")

# Names of snps and outcomes
snps       <- names(cnp.df.pc)[12:88]
covar      <- names(cnp.df.pc)[c(2,3,90:94)]
outcomes   <- names(cnp.df.pc)[5:11]
pred_all   <- c(snps,covar)

# Rename swedish outcomes to match cnp
names(swe.df.pc)[8:14] <- outcomes 
# Some proxy snps used, need to rename
names(swe.df.pc)[71:73] <- c("rs2535627", "rs4388249", "rs2905426")
names(swe.df.pc)[c(6,7,94:98)] <- covar

# Remove subjects with missing outcomes in CNP
cnp.df     <- cnp.df.pc[(complete.cases(cnp.df.pc$CRT_TIME1)),]
swe.df     <- swe.df.pc 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(pred_all)) %>% as.matrix() 
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

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
raw_rep1 <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),] #rm NAs
raw_rep  <- raw_rep[complete.cases(raw_rep1[,c(covar)]),] # N = 359

rep.predictors <- raw_rep %>% dplyr::select(one_of(pred_all)) %>% as.matrix
rep.outcomes   <- raw_rep %>% dplyr::select(one_of(outcomes))

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW   <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW    <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398



# Apply scaling transformations from CNP
rep.outcomes.s   <- (rep.outcomes - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))


rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))


##### ##### ##### ##### ##### 
##### Mixed effects models controlling for family ID
##### ##### ##### ##### ##### 

rf.pred  <- rf.pred.tr %>% as.data.frame
glm.pred <- glm.pred.tr %>% as.data.frame
mix <- cbind(rep.outcomes.s, rf.pred, glm.pred, raw_rep[,"FID.x"]) %>% as.data.frame

outcomes.rf <- NULL  # start empty list
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7]   <- outcomes            # renaming using the list of names
names(mix)[8:14]  <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22]    <- "fid"

outstats <- matrix(nrow=7,ncol=7)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf"))
  f2 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm"))
  m1 <- lme(fixed=f1, random=~1|fid, data=mix, na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|fid, data=mix, na.action=na.omit)
  r2m1 = r2beta(m1,method='sgv')
  r2m2 = r2beta(m2,method='sgv')
  outstats[i,1] <- outcomes[i]
  outstats[i,2] <- summary(m1)$tTable[2,4]
  outstats[i,3] <- summary(m1)$tTable[2,5]
  outstats[i,4] <- r2m1[r2m1$Effect==paste0(outcomes[i],".rf"), 6]
  outstats[i,5] <- summary(m2)$tTable[2,4]
  outstats[i,6] <- summary(m2)$tTable[2,5]
  outstats[i,7] <- r2m2[r2m2$Effect==paste0(outcomes[i],".glm"), 6]
}

outstats        <- as.data.frame(outstats, stringsAsFactors = F)
outstats[,2:7]  <- lapply(outstats[,2:7], as.numeric)
names(outstats) <- c("var","tval.rf","pval.rf","r2.rf", "tval.glm","pval.glm", "r2.glm")

write.table(outstats, "revision1/covar-in-training.txt", col.names=T, row.names=F, sep="\t", quote=F)

save.image("../covar-in-training.RData")

# CNP correlations and p-values
glm.r2     <- unlist(lapply(glm.perfs, function(x) x$TrainRsquared))
rf.r2      <- unlist(lapply(rf.perfs, function(x) x$TrainRsquared))

r2   <- rbind(glm.r2,unname(rf.r2))
rval <- apply(r2,2,sqrt)
tval <- apply(r2,2,function(x) sqrt((737*x)/(1-x)))
pval <- apply(tval,2,function(x) 2*(1-abs(pt(x,df=737))))

