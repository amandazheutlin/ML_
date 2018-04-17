# Schizophrenia-linked genetic influences on cognition
# Author: AMC & ABZ, Nov 2015

#### Housekeeping
# Data manipulation and plotting
library("dplyr"); library("ggplot2"); library("RColorBrewer")

# Statistical packages
library("glmnet") # https://urldefense.proofpoint.com/v2/url?u=http-3A__web.stanford.edu_-7Ehastie_glmnet_glmnet-5Falpha.html&d=AwIGAg&c=-dg2m7zWuuDZ0MUcV7Sdqw&r=1eBIqVAhmnkZRWlCX5iU07Sqp55OtRSRubwrZkpK9R4&m=x0sYKnmPdCd6tz6RllAyYboYJJPTi-LQXgCBMesayGo&s=W3iRrgwrQgyU6pVDZrYbPv8yYV80ZJeVa7JpwoHDfqM&e= 
library("e1071"); library("caret"); library("pROC")
library("permute"); library("gbm"); library("klaR"); library("randomForest")
library("dplyr")

# Parallel Computing 
library("doMC"); registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/ML_genetics/Schizophrenia-Zheutlin")
workDir <- getwd()
dataDir <- paste0(workDir, "/dfs/") 

# Load custom functions from adam's library
source("/data/ML_genetics/functions/functions.R")


### Read in pre-processed data 
# sweden_eQTL <- read.table(paste0(dataDir,"eQTL-snps-genes.txt"), header=TRUE, as.is=TRUE)
# sweden_snp  <- read.table(paste0(dataDir,"GWAS-snps-swe-input.txt"), header=TRUE, as.is=TRUE)
# cnp         <- read.table(paste0(dataDir,"GWAS-snps-cnp-input.txt"), header=TRUE, as.is=TRUE)
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
# na_count <- sapply(cnp, function(y) sum(length(which(is.na(y)))))
# na_count <- data.frame(na_count)

cnp.df     <- cnp[(complete.cases(cnp$CRT_TIME1)),]
swe.df     <- sweden_snp 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(snps)) %>% as.matrix() 
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()
# targets$VR2DR_TOTALRAW <- targets$VR2DR_TOTALRAW +1 #someone scored 0 and you can't log 0 and it isnt normal so we added 1 to all scores (range 0-43)

## Scale variables manually (so that we can apply these params to replication set)
#  record m/sd
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)

#  apply scaling
#  log before scaling
# log.targets   <- apply(targets, 2, log) %>% as.data.frame # alt:log transform 
targets       <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame
# s.log.targets <- apply(log.targets, 2, scale) %>% as.data.frame

# get residuals
# targets <- cbind(targets, cnp.df[,2:3])
# resid <- lapply(outcomes, function(x) {
#   lm(eval(substitute(i ~ age + gender, list(i = as.name(x)))),data = cnp.df)
# })
# 
# resid_vals <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
# names(resid_vals) <- outcomes
# 
# targets <- resid_vals

### Machine Learning preparation
## Set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE,
                           selectionFunction = "oneSE")

## Hyperparameter grid search
gbmGrid    <- expand.grid(.interaction.depth = c(1,2),
                          .n.trees = seq(10,4000,by=25),
                          .shrinkage = c(0.01, 0.05),
                          .n.minobsinnode = 5)
rfGrid     <- expand.grid(.mtry = c(3,4,10)) 

# Train models

### just checking which models are now available for regression/dual use
t <- getModelInfo()
m <- list();
for (i in names(t)){
  if (t[[i]]$type != "Classification"){
    m <- c(m, t[i])
  }
}


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
# glm.loop     <- lapply(targets, function(i) glmPipeline(i))
# glm.loop.log <- lapply(log.targets, function(i) glmPipeline(i))
# glm.loop.log.s <- lapply(s.log.targets, function(i) glmPipeline(i))
glm.models   <- lapply(glm.loop, function(i) i[[1]])
glm.perfs    <- lapply(glm.loop, function(i) i[[2]])
# glm.models.log  <- lapply(glm.loop.log, function(i) i[[1]])
# glm.perfs.log   <- lapply(glm.loop.log, function(i) i[[2]])
# glm.models.log.s  <- lapply(glm.loop.log.s, function(i) i[[1]])
# glm.perfs.log.s   <- lapply(glm.loop.log.s, function(i) i[[2]])


gbmPipeline <- function(response, Xmat = predictors, grid = gbmGrid, cvpar = fitControl,...){
    set.seed(1)
    model <- train(x = Xmat, y = response,
                   method = "gbm",
                   metric = "Rsquared",
                   tuneGrid = grid,
                   trControl = cvpar)
    perf <- getTrainPerf(model)
    return(list(model,perf))
}

gbm.loop     <- lapply(targets, function(i) gbmPipeline(i))
# gbm.loop.log <- lapply(log.targets, function(i) gbmPipeline(i))
gbm.models   <- lapply(gbm.loop, function(i) i[[1]])
gbm.perfs    <- lapply(gbm.loop, function(i) i[[2]])
# gbm.models.log  <- lapply(gbm.loop.log, function(i) i[[1]])
# gbm.perfs.log   <- lapply(gbm.loop.log, function(i) i[[2]])

rfPipeline <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance=TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop     <- lapply(targets, function(i) rfPipeline(i))
# rf.loop     <- lapply(targets, function(i) rfPipeline(i))
# rf.loop.log <- lapply(log.targets, function(i) rfPipeline(i))
# rf.loop.log.s <- lapply(s.log.targets, function(i) rfPipeline(i))
rf.models   <- lapply(rf.loop, function(i) i[[1]])
rf.perfs    <- lapply(rf.loop, function(i) i[[2]])
# rf.models.log  <- lapply(rf.loop.log, function(i) i[[1]])
# rf.perfs.log   <- lapply(rf.loop.log, function(i) i[[2]])
# rf.models.log.s  <- lapply(rf.loop.log.s, function(i) i[[1]])
# rf.perfs.log.s   <- lapply(rf.loop.log.s, function(i) i[[2]])

plot(varImp(rf.models$CRT_TIME1, scale=TRUE), top=15)


# crt1.out  <- promising[[1]]
# crt1.vars <- varImp(crt1.out[[1]], scale = FALSE)[[1]]
# crt1.vars <- cbind(crt1.vars, rownames(crt1.vars)) %>% as.data.frame
# crt2.out  <- promising[[2]]
# crt2.vars <- varImp(crt2.out[[1]], scale = FALSE)[[1]]
# crt2.vars <- cbind(crt2.vars, rownames(crt2.vars)) %>% as.data.frame
# 
# c1 <- crt1.vars %>% dplyr::filter(Overall > 0)
# c2 <- crt2.vars %>% dplyr::filter(Overall > 0)

# write.csv(c1,file="crt1.csv")
# write.csv(c2,file="crt2.csv")
# 
# dplyr::setdiff(c1$`rownames(crt1.vars)`, c2$`rownames(crt1.vars)`)






#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),]

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
rep.outcomes   <- raw_rep %>% dplyr::select(one_of(outcomes))
# rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW +1

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398

# twin1 <- raw_rep %>% dplyr::filter(tvab == 1) #twin specific
# twin2 <- raw_rep %>% dplyr::filter(tvab == 2)
# twin1.pred <- twin1 %>% dplyr::select(one_of(snps)) %>% as.matrix
# twin2.pred <- twin2 %>% dplyr::select(one_of(snps)) %>% as.matrix
# twin1.out.log  <- twin1 %>% dplyr::select(CRT_TIME1, CRT_TIME2) %>%
#   apply(2, log) %>% as.data.frame
# twin2.out.log  <- twin2 %>% dplyr::select(CRT_TIME1, CRT_TIME2) %>%
#   apply(2, log) %>% as.data.frame

# Apply scaling transformations from CNP
# rep.outcomes.log <- apply(rep.outcomes, 2, log)
rep.outcomes.s   <- (rep.outcomes - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))
# rep.outcomes.log.s <- (rep.outcomes.log - means)/sds

# scale within sample
# means_swe    <- apply(rep.outcomes,2,mean,na.rm=T)
# sds_swe      <- apply(rep.outcomes,2,sd,na.rm=T)
# 
# rep.outcomes.s2   <- (rep.outcomes - means_swe)/sds_swe
# rep.outcomes.s2   <- apply(rep.outcomes.s2, 2, function(i) psych::winsor(i, trim=0.01))
# rep.outcomes.log.s2 <- (rep.outcomes.log - means_swe)/sds_swe

# rep.outcomes <- apply(rep.outcomes, 2, function(i) psych::winsor(i, trim=0.03))

gbm.pred.tr     <- lapply(gbm.models, function(i) predict(i, newdata=rep.predictors))
# gbm.pred.log.tr <- lapply(gbm.models.log, function(i) predict(i, newdata=rep.predictors))
rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
# rf.pred.log.tr    <- lapply(rf.models.log, function(i) predict(i, newdata=rep.predictors))
# rf.pred.log.s.tr  <- lapply(rf.models.log.s, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))
# glm.pred.log.tr   <- lapply(glm.models.log, function(i) predict(i, newdata=rep.predictors))
# glm.pred.log.s.tr <- lapply(glm.models.log.s, function(i) predict(i, newdata=rep.predictors))

out.cors <- list()
for (i in 1:7){
  glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  gbm.cor <- cor(gbm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  # glm.log.cor <- cor(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
  # rf.log.cor  <- cor(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
  #   glm.log.s.cor <- cor(glm.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")
  #   rf.log.s.cor  <- cor(rf.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")
  #   glm.cor2 <- cor(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
  #   rf.cor2  <- cor(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
  #   glm.log.s.cor2 <- cor(glm.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")
  #   rf.log.s.cor2  <- cor(rf.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")
  #   out.cors[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.log.s.cor, 
  #                      rf.log.s.cor, glm.cor2, rf.cor2, glm.log.s.cor2, rf.log.s.cor2)
  out.cors[[i]] <- c(glm.cor, gbm.cor, rf.cor) #, glm.log.cor, rf.log.cor)
}
names(out.cors) <- names(rep.outcomes)

out.ps <- list()
for (i in 1:7){
  glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  gbm.cor <- cor.test(gbm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  # glm.log.cor <- cor.test(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
  # rf.log.cor  <- cor.test(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
  #   glm.log.s.cor <- cor.test(glm.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")$p.value
  #   rf.log.s.cor  <- cor.test(rf.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")$p.value
  #   glm.cor2 <- cor.test(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
  #   rf.cor2  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
  #   glm.log.s.cor2 <- cor.test(glm.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")$p.value
  #   rf.log.s.cor2  <- cor.test(rf.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")$p.value
  #   out.ps[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.log.s.cor, 
  #                    rf.log.s.cor, glm.cor2, rf.cor2, glm.log.s.cor2, rf.log.s.cor2)
  out.ps[[i]] <- c(glm.cor, gbm.cor, rf.cor) #, glm.log.cor, rf.log.cor)
}
names(out.ps) <- names(rep.outcomes)

# correlations controlling for FID for models that work
library(nlme)

mix <- cbind(rep.outcomes.s, rf.pred.tr, glm.pred.tr, raw_rep[,c(2,6,7)]) %>% as.data.frame

outcomes.rf <- NULL
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7] <- outcomes
names(mix)[8:14] <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22:24] <- c("fid","age","sex")

outstats <- matrix(,nrow=7,ncol=5)
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

# CNP correlations and p-values
glm.r2     <- unlist(lapply(glm.perfs, function(x) x$TrainRsquared))
rf.r2      <- unlist(lapply(rf.perfs, function(x) x$TrainRsquared))
# glm.log.r2 <- unlist(lapply(glm.perfs.log, function(x) x$TrainRsquared))
# rf.log.r2  <- unlist(lapply(rf.perfs.log, function(x) x$TrainRsquared))

# r2 <- rbind(glm.r2,unname(rf.r2),unname(glm.log.r2),unname(rf.log.r2))
r2 <- rbind(glm.r2,unname(rf.r2))
rval <- apply(r2,2,sqrt)
tval <- apply(r2,2,function(x) sqrt((737*x)/(1-x)))
pval <- apply(tval,2,function(x) 2*(1-abs(pt(x,df=737))))

# to load later
save.image(file=paste0(workDir,"amanda.results4.Rdata"))

#################### figure
library(ggplot2); library(stringr)

# snp-gene table
code <- read.csv("scz2.anneal.108.genes.csv")
proxy <- read.table("index-proxy-table_allsnps.txt",header=T)
proxy$gene <- code[match(proxy$index,code$bestsnp),8]

# trails1
fig_trails1 <- caret::varImp(rf.models$CRT_TIME1, scale=TRUE)[[1]] %>% as.data.frame
fig_trails1 <- cbind(rownames(fig_trails1), fig_trails1) %>% arrange(-Overall)
names(fig_trails1) <- c("SNP", "Importance")
fig_trails1$Gene <- proxy[match(fig_trails1$SNP,proxy$proxy),3]
fig_trails1$SNP <- factor(fig_trails1$SNP,levels=fig_trails1$SNP[order(fig_trails1$Importance)])

# visual reproduction - delayed recall
fig_vrdr <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr) <- c("SNP", "Importance")
fig_vrdr$Gene <- proxy[match(fig_vrdr$SNP,proxy$proxy),3]
fig_vrdr$SNP <- factor(fig_vrdr$SNP,levels=fig_vrdr$SNP[order(fig_vrdr$Importance)])

# code by shared vs. task-specific genes
fig_trails1$shared <- ifelse(fig_trails1$Gene %in% fig_vrdr[1:20,3],1,2) %>% as.factor()
fig_vrdr$shared <- ifelse(fig_vrdr$Gene %in% fig_trails1[1:20,3],1,2) %>% as.factor()
levels(fig_trails1$shared) <- c("both tasks","task-specific")
levels(fig_vrdr$shared) <- c("both tasks","task-specific")

ggplot(fig_trails1[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_trails1[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=13, face="bold"),
        axis.text = element_text(size=9,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(50,100))

ggplot(fig_vrdr[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_vrdr[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=13, face="bold"),
        axis.text = element_text(size=9,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(50,100))

# cor(gbm.pred.tr[[1]], rep.outcomes[,1])
# cor(rf.pred.tr[[1]], rep.outcomes[,1])
# cor(gbm.pred.log.tr[[1]], rep.outcomes.log[,1])
# cor(rf.pred.log.tr[[1]], rep.outcomes.log[,1])
# 
# cor(gbm.pred.tr[[2]], rep.outcomes[,2])
# cor(rf.pred.tr[[2]], rep.outcomes[,2])
# cor(gbm.pred.log.tr[[2]], rep.outcomes.log[,2])
# cor(rf.pred.log.tr[[2]], rep.outcomes.log[,2])


# Twin specific
# t1.gbm.pred.log.tr <- lapply(list(gbm.models.log[[5]],gbm.models.log[[6]]),
#                              function(i) predict(i, newdata=twin1.pred))
# t1.rf.pred.log.tr  <- lapply(list(rf.models.log[[5]],rf.models.log[[6]]),
#                              function(i) predict(i, newdata=twin1.pred))
# t2.gbm.pred.log.tr <- lapply(list(gbm.models.log[[5]],gbm.models.log[[6]]),
#                              function(i) predict(i, newdata=twin2.pred))
# t2.rf.pred.log.tr  <- lapply(list(rf.models.log[[5]],rf.models.log[[6]]),
#                              function(i) predict(i, newdata=twin2.pred))
# cor(t1.gbm.pred.log.tr[[1]], twin1.out.log[,1])
# cor(t1.rf.pred.log.tr[[1]],  twin1.out.log[,1])
# cor(t2.gbm.pred.log.tr[[1]], twin2.out.log[,1])
# cor(t2.rf.pred.log.tr[[1]],  twin2.out.log[,1])
# 
# 
# 
# rf.varimp.log.t1  <- lapply(rf.models.log, varImp)
# gbm.varimp.log.t1 <- lapply(gbm.models.log, varImp)
# plot(rf.varimp.log.t1[[5]])
# plot(gbm.varimp.log.t1[[5]])
# 
