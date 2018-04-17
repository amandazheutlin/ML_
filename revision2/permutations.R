# Schizophrenia-linked genetic influences on cognition
# Authors: AMC & ABZ, September 2017

# load libraries and main analyses
mldir <- "/data/ML_genetics/Schizophrenia-Zheutlin"
setwd(mldir)

libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "psych", "dplyr", "nlme", "psych","stringr", "r2glmm")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

load("revision2/covar-resid-targets-10PCs-g.Rdata")

#### Permutation testing for the RF models

rfPipelinePerm <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1000)
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
    set.seed(perm+1000)
    permutedTarget <- sample(target)
    mod <- rfPipelinePerm(permutedTarget, Xmat = xmat)[[1]]
#     names(mod) <- paste0(mod, perm)
    out_perm[[perm]] <- mod
  }
  return(out_perm)
}

# testing the formulas
# perm.test     <- lapply(select(targets, CRT_TIME1, VR2DR_TOTALRAW), function(i) bootMe(2, i, predictors))
# CRT1_TIME1.models     <- perm.test[[1]]
# VR2DR_TOTALRAW.models <- perm.test[[2]]

### loop over bootstraps and phenotypes
set.seed(1000)
perm.loop <- lapply(select(targets, CRT_TIME1, VR2DR_TOTALRAW), function(i) bootMe(250, i, predictors))
CRT1_TIME1.models     <- perm.loop[[1]]
VR2DR_TOTALRAW.models <- perm.loop[[2]]

### get training performance
crt1.trainmod  <- lapply(CRT1_TIME1.models, function(i) getTrainPerf(i))
crt1.trainperf <- list() 
for (i in 1:250){ crt1.trainperf[[i]] <- crt1.trainmod[[i]][[2]]}
crt1.trainperf <- unlist(crt1.trainperf)
summary(crt1.trainperf)

vr2d.trainmod  <- lapply(VR2DR_TOTALRAW.models, getTrainPerf)
vr2d.trainperf <- list() 
for (i in 1:250){ vr2d.trainperf[[i]] <- vr2d.trainmod[[i]][[2]]}
vr2d.trainperf <- unlist(vr2d.trainperf)
summary(vr2d.trainperf)


### get predictions and test in validation sample
# trails
p.crt1.predictions        <- lapply(CRT1_TIME1.models, function(x) predict(x, newdata = rep.predictors))
names(p.crt1.predictions) <- paste0("CRT_TIME1.", 1:250)

p.mix.crt1 <- cbind(rep.targets, p.crt1.predictions, raw_rep[,"FID.x"]) %>% as.data.frame
names(p.mix.crt1)[length(names(p.mix.crt1))] <- "fid"

permstats.p.crt1 <- matrix(nrow=250,ncol=4)
for (perm in 1:250){
  f1.crt1 <- formula(paste0("CRT_TIME1 ~ CRT_TIME1.", perm) )
  m1.crt1 <- nlme::lme(fixed=f1.crt1, random = ~1 | fid, data = p.mix.crt1)
  r2m1 = r2beta(m1.crt1,method='sgv')
  permstats.p.crt1[perm,1] <- summary(m1.crt1)$tTable[2,4]
  permstats.p.crt1[perm,2] <- summary(m1.crt1)$tTable[2,5]
  permstats.p.crt1[perm,3] <- r2m1[r2m1$Effect==paste0("CRT_TIME1.", perm), 6]
  permstats.p.crt1[perm,4] <- sqrt(mean(m1.crt1$residuals^2))
}

permstats.p.crt1        <- as.data.frame(permstats.p.crt1)  # 4 p < .05; 4/101 = .0396
names(permstats.p.crt1) <- c("tval","pval", "R2", "RMSE")   # round 1: 3 R2 > .013; 3/101 = .0297 (all pos R)
                                                            # round 2: 4 (pos R) R2 > .013; 7/201 = .035
                                                            # round 3: 8 (pos R) R2 > .013; 15/501 = .030
                                                            # round 4: 7 (pos R) R2 > .013; 22/751 = .029
                                                            # round 5: 5 (pos R) R2 > .013; 27/1001 = .027
permstats.p.crt1        <- permstats.p.crt1[order(-permstats.p.crt1$R2),]
permstats.p.crt1$R      <- sqrt(permstats.p.crt1$R2)
permstats.p.crt1$R      <- ifelse(permstats.p.crt1$tval < 1, permstats.p.crt1$R * -1, permstats.p.crt1$R)

# vr2d
p.vr2d.predictions        <- lapply(VR2DR_TOTALRAW.models, function(x) predict(x, newdata = rep.predictors))
names(p.vr2d.predictions) <- paste0("VR2DR_TOTALRAW.", 1:250)

p.mix.vr2d <- cbind(rep.targets, p.vr2d.predictions, raw_rep[,"FID.x"]) %>% as.data.frame
names(p.mix.vr2d)[length(names(p.mix.vr2d))] <- "fid"

permstats.p.vr2d <- matrix(nrow=250,ncol=4)
for (perm in 1:250){
  f1.vr2d <- formula(paste0("VR2DR_TOTALRAW ~ VR2DR_TOTALRAW.", perm) )
  m1.vr2d <- nlme::lme(fixed=f1.vr2d, random = ~1 | fid, control = lmeControl(opt = "optim"), data = p.mix.vr2d)
  r2m1    <- r2beta(m1.vr2d,method='sgv')
  permstats.p.vr2d[perm,1] <- summary(m1.vr2d)$tTable[2,4]
  permstats.p.vr2d[perm,2] <- summary(m1.vr2d)$tTable[2,5]
  permstats.p.vr2d[perm,3] <- r2m1[r2m1$Effect==paste0("VR2DR_TOTALRAW.", perm), 6]
  permstats.p.vr2d[perm,4] <- sqrt(mean(m1.vr2d$residuals^2))
}

permstats.p.vr2d        <- as.data.frame(permstats.p.vr2d)   # 9 p < .05; 9/101 = .0891
names(permstats.p.vr2d) <- c("tval","pval","R2","RMSE")      # 9 (pos R) R2 > .010; 9/101 = .0891
                                                             # round 2: 6 (pos R) R2 > .010; 15/201 = .075
                                                             # round 3: 19 (pos R) R2 > .010; 34/501 = .068
                                                             # round 4: 10 (pos R) R2 > .010; 44/751 = .059
                                                             # round 5: 12 (pos R) R2 > .010; 56/1001 = .056
permstats.p.vr2d        <- permstats.p.vr2d[order(-permstats.p.vr2d$R2),]
permstats.p.vr2d$R      <- sqrt(permstats.p.vr2d$R2)
permstats.p.vr2d$R      <- ifelse(permstats.p.vr2d$tval < 1, permstats.p.vr2d$R * -1, permstats.p.vr2d$R)

# save results
write.table(permstats.p.crt1,"revision2/crt1_perm_results1250.txt")
write.table(permstats.p.vr2d,"revision2/vr_perm_results1250.txt")

### training performance
crt.train <- read.table("revision2/crt1_perm_train.txt", header=T, stringsAsFactors = F)
vr.train  <- read.table("revision2/vr_perm_train.txt", header=T, stringsAsFactors = F)

weighted.mean(crt.train$Mean, crt.train$PermN)
weighted.mean(vr.train$Mean, vr.train$PermN)

### graph permutations
# concatenate results
perm_files_crt <- list.files(path = "revision2/.", pattern="crt1_perm_results*")
perm_crt       <- lapply(perm_files_crt, function(x) read.table(paste0("revision2/", x),header=T))

perm_crt_all <- NULL
for (i in 1: length(perm_crt)){
  perm_crt_all <- rbind(perm_crt_all, perm_crt[[i]])
}

perm_files_vr <- list.files(path = "revision2/.", pattern="vr_perm_results*")
perm_vr       <- lapply(perm_files_vr, function(x) read.table(paste0("revision2/", x),header=T))

perm_vr_all <- NULL
for (i in 1: length(perm_vr)){
  perm_vr_all <- rbind(perm_vr_all, perm_vr[[i]])
}

# graphs
crt1.r <- sqrt(.013)

ggplot(perm_crt_all,aes(x=tval)) +                   # 6in x 5in
  geom_histogram(binwidth=.35,aes(fill=..count..)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_segment(x=2.09,y=0,xend=2.09,yend=150,linetype="dashed") +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("Permuted Model Prediction Strength (t)")

vr2d.r <- sqrt(.010)

ggplot(perm_vr_all,aes(x=tval)) +                    # 6in x 5in
  geom_histogram(binwidth=.35,aes(fill=..count..)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_segment(x=1.81,y=0,xend=1.81,yend=150,linetype="dashed") +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("Permuted Model Prediction Strength (t)")


save.image(file = "revision2/rf.permutations_covar.RData")
save.image(file = "revision2/rf.permutations_covar_r2.RData")
save.image(file = "revision2/rf.permutations_covar_r3.RData")
save.image(file = "revision2/rf.permutations_covar_r4.RData")
save.image(file = "revision2/rf.permutations_covar_r5.RData")
save.image(file = "revision2/rf.permutations_covar_r6.RData")


