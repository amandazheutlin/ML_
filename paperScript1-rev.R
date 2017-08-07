# Schizophrenia-linked genetic influences on cognition

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "pROC", "permute", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "e1071", "dplyr", "nlme", 
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
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

#  apply scaling
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)
targets  <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame

# recode genotypes to dosage of minor allele (continuous)
# 1/1 = 1, 1/2 = 2, 2/2 = 3
# predictors <- apply(predictors,2,function(x) ifelse(x==11,1,ifelse(x==12,2,3)))

## Scale variables manually (so that we can apply these params to replication set)
#  scale based on controls
#  record m/sd
# targets.c <- cnp.df[cnp.df$Group=="control",] %>% dplyr::select(one_of(outcomes)) %>% 
#   apply(2,as.numeric) %>% as.data.frame()
# means     <- apply(targets.c,2,mean)
# sds       <- apply(targets.c,2,sd)
# 
# #  apply scaling
# targets   <- (targets - means)/sds

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


# check performance in training with covariates
# rf.pred.cnp  <- lapply(rf.models, function(i) predict(i, newdata=predictors))
# glm.pred.cnp <- lapply(glm.models, function(i) predict(i, newdata=predictors))
# 
# mix.cnp      <- cbind(targets, rf.pred.cnp, glm.pred.cnp, cnp.df[,c("ptid", "age", "gender")]) %>% as.data.frame
# outcomes.rf  <- NULL  # start empty list
# outcomes.glm <- NULL
# for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
# for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}
# 
# names(mix.cnp)[1:7]   <- outcomes            # renaming using the list of names
# names(mix.cnp)[8:14]  <- outcomes.rf
# names(mix.cnp)[15:21] <- outcomes.glm
# 
# cnp.pc     <- read.table("/data/swe_gwas/ABZ/ML/cnp_MLsubs_pca.eigenvec", stringsAsFactors = F)
# mix.cnp.pc <- merge(mix.cnp, cnp.pc, by.x = "ptid", by.y = "V1")
# 
# outstats.cnp <- matrix(nrow=7,ncol=5)
# for (i in 1:length(outcomes)){
#   f3 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + gender + V3 + V4 + V5 + V6 + V7"))
#   f4 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + gender + V3 + V4 + V5 + V6 + V7"))
#   m1 <- lm(f3, data=mix.cnp.pc, na.action=na.omit)
#   m2 <- lm(f4, data=mix.cnp.pc, na.action=na.omit)
#   rf.tval  <- summary(m1)$coefficients[2,3]
#   glm.tval <- summary(m2)$coefficients[2,3]
#   outstats.cnp[i,1] <- outcomes[i]
#   outstats.cnp[i,2] <- rf.tval
#   outstats.cnp[i,3] <- summary(m1)$coefficients[2,4]
#   outstats.cnp[i,4] <- sqrt(rf.tval^2/(rf.tval^2+summary(m1)$df[2]))
#   outstats.cnp[i,5] <- glm.tval
#   outstats.cnp[i,6] <- summary(m2)$coefficients[2,4]
#   outstats.cnp[i,7] <- sqrt(glm.tval^2/(glm.tval^2+summary(m2)$df[2]))
# }
# 
# 
# outstats.cnp        <- as.data.frame(outstats.cnp, stringsAsFactors = F)
# outstats.cnp[,2:7]  <- lapply(outstats.cnp[,2:7], as.numeric)
# names(outstats.cnp) <- c("var","tval.rf","pval.rf","r2.rf", "tval.glm","pval.glm", "r2.rf")


#### Permutation testing for the RF models

rfPipelinePerm <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
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

bootMe <- function(permutations, target, xmat = predictors,...){
  out_perm <- list()
  for (perm in 1:permutations){
    if(perm %% 5 == 0){print(paste((perm - 1), "perms done so far"))}
    set.seed(perm)
    permutedTarget <- target[shuffle(length(target))]
    mod <- rfPipelinePerm(permutedTarget, Xmat = xmat)[[1]]
    names(mod) <- paste0(mod, perm)
    out_perm[[perm]] <- mod
  }
  return(out_perm)
}



# set.seed(1)
# test <- bootMe(2, targets$CVLT_TOTCOR, predictors)

# 
# set.seed(1)
# perm.loop     <- lapply(select(targets, -CRT_TIME2), function(i) bootMe(2, i, predictors))

### loop over bootstraps and phenotypes
set.seed(1)
perm.loop     <- lapply(select(targets, CRT_TIME1, VR2DR_TOTALRAW), function(i) bootMe(100, i, predictors))
CRT1_TIME1.models     <- perm.loop[[1]]
VR2DR_TOTALRAW.models <- perm.loop[[2]]

### get training performance
crt1.trainmod  <- lapply(CRT1_TIME1.models, getTrainPerf)
crt1.trainperf <- list() 
for (i in 1:100){ crt1.trainperf[[i]] <- crt1.trainmod[[i]][[2]]}
crt1.trainperf <- unlist(crt1.trainperf)
summary(crt1.trainperf)

vr2d.trainmod  <- lapply(VR2DR_TOTALRAW.models, getTrainPerf)
vr2d.trainperf <- list() 
for (i in 1:100){ vr2d.trainperf[[i]] <- vr2d.trainmod[[i]][[2]]}
vr2d.trainperf <- unlist(vr2d.trainperf)
summary(vr2d.trainperf)


### get predictions and test in validation sample
# trails
p.crt1.predictions  <- lapply(CRT1_TIME1.models, function(x) predict(x, newdata = rep.predictors))
names(p.crt1.predictions) <- paste0("CRT_TIME1.", 1:100)

p.mix.crt  <- cbind(rep.outcomes.s, p.crt1.predictions, raw_rep[,c("IID", "FID", "age", "sex")]) %>% as.data.frame
p.mix.crt1 <- merge(p.mix.crt, PCs, by = "IID")

permstats.p.crt1 <- matrix(nrow=100,ncol=3)
for (perm in 1:100){
  f1.crt1 <- formula(paste0("CRT_TIME1 ~ CRT_TIME1.", perm, "  + age + sex + C1 + C2 + C3 + C4 + C5" ) )
  m1.crt1 <- nlme::lme(fixed=f1.crt1, random = ~1 | FID.x, data = p.mix.crt1, na.action = na.omit)
  r2m1 = r2beta(m1.crt1,method='sgv')
  permstats.p.crt1[perm,1] <- summary(m1.crt1)$tTable[2,4]
  permstats.p.crt1[perm,2] <- summary(m1.crt1)$tTable[2,5]
  permstats.p.crt1[perm,3] <- r2m1[r2m1$Effect==paste0("CRT_TIME1.", perm), 6]
}

permstats.p.crt1 <- as.data.frame(permstats.p.crt1) # 7 p < .05; 7/101 = .0693
names(permstats.p.crt1) <- c("tval","pval", "R2")   # 1 t > 2.75; 1/101 = .0099
permstats.p.crt1$R      <- sqrt(permstats.p.crt1$R2)
permstats.p.crt1 <- permstats.p.crt1[order(permstats.p.crt1$pval),]

# vr2d
p.vr2d.predictions  <- lapply(VR2DR_TOTALRAW.models, function(x) predict(x, newdata = rep.predictors))
names(p.vr2d.predictions) <- paste0("VR2DR_TOTALRAW.", 1:100)

p.mix.vr2d.1 <- cbind(rep.outcomes.s, p.vr2d.predictions, raw_rep[,c("IID", "FID", "age", "sex")]) %>% as.data.frame
p.mix.vr2d   <- merge(p.mix.vr2d.1, PCs, by = "IID")

permstats.p.vr2d <- matrix(nrow=100,ncol=2)
for (perm in 1:100){
  f1.vr2d <- formula(paste0("VR2DR_TOTALRAW ~ VR2DR_TOTALRAW.", perm, "  + age + sex + C1 + C2 + C3 + C4 + C5" ) )
  m1.vr2d <- nlme::lme(fixed=f1.vr2d, random = ~1 | FID.x, control = lmeControl(opt = "optim"), 
                       data = p.mix.vr2d, na.action = na.omit)
  permstats.p.vr2d[perm,1] <- summary(m1.vr2d)$tTable[2,4]
  permstats.p.vr2d[perm,2] <- summary(m1.vr2d)$tTable[2,5]
}

permstats.p.vr2d <- as.data.frame(permstats.p.vr2d) # 8 p < .05; 8/101 = .0792
names(permstats.p.vr2d) <- c("tval","pval")         # 10 t > 1.73; 10/101 = .0990
permstats.p.vr2d <- permstats.p.vr2d[order(permstats.p.vr2d$pval),]
  
save.image(file = "permutations_for_amanda.RData")

### graph permutations

ggplot(permstats.p.crt1,aes(x=tval)) + 
  geom_histogram(binwidth=.6,aes(fill=..count..)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_segment(x=2.75,y=0,xend=2.75,yend=25,linetype="dashed") +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("Prediction Strength of Permutations")

ggplot(permstats.p.vr2d,aes(x=tval)) + 
  geom_histogram(binwidth=.55,aes(fill=..count..)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_segment(x=1.73,y=0,xend=1.73,yend=25,linetype="dashed") +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("Prediction Strength of Permutations")


#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),] #rm NAs

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
# rep.predictors <- apply(rep.predictors,2,function(x) ifelse(x==11,1,ifelse(x==12,2,3)))
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

# out.cors <- list()
# for (i in 1:7){
#   glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
#   rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
#   out.cors[[i]] <- c(glm.cor, rf.cor) 
# }
# names(out.cors) <- names(rep.outcomes)
# 
# out.ps <- list()
# for (i in 1:7){
#   glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
#   rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
#   out.ps[[i]] <- c(glm.cor, rf.cor)
# }
# names(out.ps) <- names(rep.outcomes)
# 

##### ##### ##### ##### ##### 
##### Mixed effects models controlling for family ID
##### ##### ##### ##### ##### 

mix <- cbind(rep.outcomes.s, rf.pred.tr, glm.pred.tr, raw_rep[,c("FID", "IID", "age", "sex", "dx")]) %>% as.data.frame

outcomes.rf <- NULL  # start empty list
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7]   <- outcomes            # renaming using the list of names
names(mix)[8:14]  <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22:26] <- c("fid","iid", "age","sex", "dx")

outstats <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + sex"))
  f2 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + sex"))
  m1 <- lme(fixed=f1, random=~1|fid,data=mix,na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|fid,data=mix,na.action=na.omit)
  outstats[i,1] <- outcomes[i]
  outstats[i,2] <- summary(m1)$tTable[2,4]
  outstats[i,3] <- summary(m1)$tTable[2,5]
  outstats[i,4] <- summary(m2)$tTable[2,4]
  outstats[i,5] <- summary(m2)$tTable[2,5]
}

outstats <- as.data.frame(outstats)
names(outstats) <- c("var","tval.rf","pval.rf","tval.glm","pval.glm")

# add PCs
PCs    <- read.table("/data/swe_gwas/ABZ/ML/swe_MLsubs2.mds", header=T, stringsAsFactors = F)
mix.pc <- merge(mix, PCs, by.x = "iid", by.y = "IID")

outstats.pc <- matrix(nrow=7,ncol=11)
for (i in 1:length(outcomes)){
  f3 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + sex + C1 + C2 + C3 + C4 + C5"))
  f4 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + sex + C1 + C2 + C3 + C4 + C5"))
  m1 <- lme(fixed=f3, random=~1|fid,data=mix.pc,na.action=na.omit)
  m2 <- lme(fixed=f4, random=~1|fid,data=mix.pc,na.action=na.omit)
  r2m1 = r2beta(m1,method='sgv')
  r2m2 = r2beta(m2,method='sgv')
  outstats.pc[i,1] <- outcomes[i]
  outstats.pc[i,2] <- summary(m1)$tTable[2,4]
  outstats.pc[i,3] <- summary(m1)$tTable[2,5]
  outstats.pc[i,4] <- r2m1[r2m1$Effect==paste0(outcomes[i],".rf"), 6]
  outstats.pc[i,5] <- sqrt(mean(m1$residuals^2))
  outstats.pc[i,6] <- summary(m2)$tTable[2,4]
  outstats.pc[i,7] <- summary(m2)$tTable[2,5]
  outstats.pc[i,8] <- r2m2[r2m2$Effect==paste0(outcomes[i],".glm"), 6]
  outstats.pc[i,9] <- sqrt(mean(m2$residuals^2))
  outstats.pc[i,10] <- mean(m1$residuals)
  outstats.pc[i,11] <- mean(m2$residuals)
}

outstats.pc        <- as.data.frame(outstats.pc, stringsAsFactors = F)
outstats.pc[,2:11]  <- lapply(outstats.pc[,2:11], as.numeric)
names(outstats.pc) <- c("var","tval.rf","pval.rf","r2.rf", "rmse.rf", 
                        "tval.glm","pval.glm", "r2.glm", "rmse.glm",
                        "mae.rf", "mae.glm")

write.table(outstats.pc, "revision1/originalresults_pcs.txt", col.names=T, row.names=F, sep="\t", quote=F)

outstats.rmse <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f3 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf"))
  f4 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm"))
  m1 <- lme(fixed=f3, random=~1|fid, data=mix.pc, control = lmeControl(opt = "optim"), na.action=na.omit)
  m2 <- lme(fixed=f4, random=~1|fid, data=mix.pc, control = lmeControl(opt = "optim"), na.action=na.omit)
  outstats.rmse[i,1] <- outcomes[i]
  outstats.rmse[i,2] <- sqrt(mean(m1$residuals^2))
  outstats.rmse[i,3] <- sqrt(mean(m2$residuals^2))
  outstats.rmse[i,4] <- mean(m1$residuals)
  outstats.rmse[i,5] <- mean(m2$residuals)
}

outstats.rmse        <- as.data.frame(outstats.rmse, stringsAsFactors = F)
names(outstats.rmse) <- c("var", "rmse.rf", "rmse.glm", "mae.rf", "mae.glm")


outstats.rmse <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f3 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf"))
  f4 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm"))
  m1 <- lm(f3, data=mix.pc, na.action=na.omit)
  m2 <- lm(f4, data=mix.pc, na.action=na.omit)
  outstats.rmse[i,1] <- outcomes[i]
  outstats.rmse[i,2] <- sqrt(mean(m1$residuals^2))
  outstats.rmse[i,3] <- sqrt(mean(m2$residuals^2))
  outstats.rmse[i,4] <- mean(m1$residuals)
  outstats.rmse[i,5] <- mean(m2$residuals)
}

outstats.rmse        <- as.data.frame(outstats.rmse, stringsAsFactors = F)
names(outstats.rmse) <- c("var", "rmse.rf", "rmse.glm", "mae.rf", "mae.glm")


save.image("../results5_amended.RData")

## graph predictions and observed values 
# trails A
ggplot(mix, aes(x=CRT_TIME1,y=CRT_TIME1.rf)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_point() +
  geom_smooth(method=lm) + 
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("Predicted Scores (Random Forest)") +
  xlab("Observed Scores (Scaled)")

ggplot(mix, aes(x=VR2DR_TOTALRAW,y=VR2DR_TOTALRAW.rf)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_point() +
  geom_smooth(method=lm) + 
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("Predicted Scores (Random Forest)") +
  xlab("Observed Scores (Scaled)")

# in each patient group
mix.sz <- mix[mix$dx==1,]
mix.bp <- mix[mix$dx==3,]

outstats.sz <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + sex"))
  f2 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + sex"))
  m1 <- lme(fixed=f1, random=~1|fid,data=mix.sz,na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|fid,data=mix.sz,na.action=na.omit)
  outstats.sz[i,1] <- outcomes[i]
  outstats.sz[i,2] <- summary(m1)$tTable[2,4]
  outstats.sz[i,3] <- summary(m1)$tTable[2,5]
  outstats.sz[i,4] <- summary(m2)$tTable[2,4]
  outstats.sz[i,5] <- summary(m2)$tTable[2,5]
}

outstats.sz <- as.data.frame(outstats.sz)
names(outstats.sz) <- c("var","tval.rf","pval.rf","tval.glm","pval.glm")

outstats.bp <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + sex"))
  f2 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + sex"))
  m1 <- lme(fixed=f1, random=~1|fid,data=mix.bp,na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|fid,data=mix.bp,na.action=na.omit)
  outstats.bp[i,1] <- outcomes[i]
  outstats.bp[i,2] <- summary(m1)$tTable[2,4]
  outstats.bp[i,3] <- summary(m1)$tTable[2,5]
  outstats.bp[i,4] <- summary(m2)$tTable[2,4]
  outstats.bp[i,5] <- summary(m2)$tTable[2,5]
}

outstats.bp <- as.data.frame(outstats.bp)
names(outstats.bp) <- c("var","tval.rf","pval.rf","tval.glm","pval.glm")

cor.sz <- cor(mix.sz,use="complete.obs") %>% as.data.frame

# CNP correlations and p-values
glm.r2     <- unlist(lapply(glm.perfs, function(x) x$TrainRsquared))
rf.r2      <- unlist(lapply(rf.perfs, function(x) x$TrainRsquared))

r2   <- rbind(glm.r2,unname(rf.r2))
rval <- apply(r2,2,sqrt)
tval <- apply(r2,2,function(x) sqrt((737*x)/(1-x)))
pval <- apply(tval,2,function(x) 2*(1-abs(pt(x,df=737))))

# general ability check 
# calculated from all tasks excepting model-specific task
# scaled, winsorized variables 

# all 6 tasks
g_all <- principal(mix.pc[,c(2:4,6:8)], nfactors=1, rotation="none",scores=T)
g_all # fit = .97, proportion var = .61
mix.pc$g_all <- as.numeric(g_all$scores)

# leave out trails1
g_tr1 <- principal(mix.pc[,c(2:3,6:8)], nfactors=1, rotation="none",scores=T)
g_tr1 # fit = .98, proportion var = .72
mix.pc$g_tr1 <- as.numeric(g_tr1$scores)

# leave out vr2
g_vr <- principal(mix.pc[,c(2:4,6:7)], nfactors=1, rotation="none",scores=T)
g_vr # fit = .96, proportion var = .59
mix.pc$g_vr <- as.numeric(g_vr$scores)

# mixed effect models with g vars
mx.tr1.g <- lme(g_tr1 ~ CRT_TIME1.rf + age + sex + C1 + C2 + C3 + C4 + C5, random = ~1 | fid, data = mix.pc, na.action = na.omit)
mx.vr.g <- lme(g_vr ~ VR2DR_TOTALRAW.rf + age + sex + C1 + C2 + C3 + C4 + C5, random = ~1 | fid, data = mix.pc, na.action = na.omit)

# to load later
save.image(file=paste0(workDir,"results5_amended.Rdata"))

#################### figure
load(paste0(workDir,"rev_results.Rdata"))
library(ggplot2); library(stringr); library(dplyr)
     
# snp-gene table
code       <- read.csv("scz2.anneal.108.genes.csv")
proxy      <- read.table("index-proxy-table_allsnps.txt",header=T)
proxy$gene <- code[match(proxy$index,code$bestsnp),8]

# trails1
fig_trails1        <- caret::varImp(rf.models$CRT_TIME1, scale=TRUE)[[1]] %>% as.data.frame
fig_trails1        <- cbind(rownames(fig_trails1), fig_trails1) %>% arrange(-Overall)
names(fig_trails1) <- c("SNP", "Importance")
fig_trails1$Gene   <- proxy[ match(fig_trails1$SNP, proxy$proxy), 3]
fig_trails1$SNP    <- factor(fig_trails1$SNP,
                             levels = fig_trails1$SNP[order(fig_trails1$Importance)])
fig_trails1$Gene2 <- fig_trails1$Gene %>% as.character()
fig_trails1[1,5] <- "C1orf132, CD46*"
fig_trails1[3,5] <- "GLT8D1, GNL3*"
fig_trails1[11,5] <- "AC005609.1, CD14*"
fig_trails1[13,5] <- "C2orf82, EFHD1*"
fig_trails1[17,5] <- "ACD, C16orf86*"
fig_trails1[20,5] <- "CILP2, GATAD2A*"

# visual reproduction - delayed recall
fig_vrdr           <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr           <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr)    <- c("SNP", "Importance")
fig_vrdr$Gene      <- proxy[match(fig_vrdr$SNP, proxy$proxy), 3]
fig_vrdr$SNP       <- factor(fig_vrdr$SNP,
                             levels = fig_vrdr$SNP[order(fig_vrdr$Importance)])
fig_vrdr$Gene2 <- fig_vrdr$Gene %>% as.character()
fig_vrdr[1,5] <- "C4orf27, CLCN3*"
fig_vrdr[7,5] <- "AC005609.1, CD14*"
fig_vrdr[8,5] <- "GLT8D1, GNL3*"
fig_vrdr[14,5] <- "ADAMTSL3, GOLGA6L4*"
fig_vrdr[15,5] <- "SGSM2, SMG6*"
fig_vrdr[20,5] <- "C2orf82, EFHD1*"
  
# code by shared vs. task-specific genes
fig_trails1$shared         <- ifelse(fig_trails1$Gene %in% fig_vrdr[1:20,3],1,2) %>% as.factor()
fig_vrdr$shared            <- ifelse(fig_vrdr$Gene %in% fig_trails1[1:20,3],1,2) %>% as.factor()
levels(fig_trails1$shared) <- c("both tasks","task-specific")
levels(fig_vrdr$shared)    <- c("both tasks","task-specific")

# trails 1 figure
ggplot(fig_trails1[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.85) +
  geom_text(stat="identity",y=46, hjust=0, aes(label=Gene2), colour="white", size=4) +
  scale_fill_manual(values=c("#B21212","#1485CC")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.06,vjust=1.8,size=18,face="bold")) +
  ggtitle("A) Trails 1/A") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(45,100))

# vis reprod figure
ggplot(fig_vrdr[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.85) +
  geom_text(stat="identity",y=46, hjust=0, aes(label=Gene2), colour="white", size=4) +
  scale_fill_manual(values=c("#B21212","#1485CC")) +  
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.06,vjust=1.8,size=18,face="bold")) +
  ggtitle("B) Visual Reproduction II") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(45,100))


### population stratification within CNP 
cnp_pcs <- read.table("cnp_MLsubs_pca.eigenvec",header=F)
cnp2    <- merge(cnp.df, cnp_pcs[,c(1,3:22)], by.x="ptid",by.y="V1")
cormat  <- cor(cnp2[,c(5:7,9:11,89:108)]) %>% as.data.frame
cormat_edit <- cormat[7:26,1:6]
summary(cormat_edit)