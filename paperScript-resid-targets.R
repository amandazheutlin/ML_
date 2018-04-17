# Schizophrenia-linked genetic influences on cognition
# Authors: AMC & ABZ, September 2017

# control analyses for revision II
# regress covariates out of predictors and outcomes 

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "psych", "dplyr", "nlme", "psych","stringr", "r2glmm")
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
covar      <- names(cnp.df.pc)[c(2,3,90:99)]
outcomes   <- names(cnp.df.pc)[5:11]
pred_all   <- c(snps,covar)

# Rename swedish outcomes to match cnp
names(swe.df.pc)[8:14] <- outcomes 
# Some proxy snps used, need to rename
names(swe.df.pc)[71:73] <- c("rs2535627", "rs4388249", "rs2905426")
names(swe.df.pc)[c(6,7,94:103)] <- covar

# Remove subjects with missing outcomes in CNP
cnp.df     <- cnp.df.pc[(complete.cases(cnp.df.pc$CRT_TIME1)),]
swe.df     <- swe.df.pc 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(snps)) %>% as.matrix() 
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

#  apply scaling
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)
targets  <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame


# create residual variables for all predictors and outcomes
# pred.tar    <- c(outcomes, snps) # things to get residuals from
resid.input <- cbind(targets,cnp.df[c(snps,covar)])
resid       <- lapply(outcomes, function(x) {
  lm(eval(substitute(i ~ age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                     list(i = as.name(x)))), data = resid.input)
})

resid_vals        <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
names(resid_vals) <- outcomes

#reset targets and predictions to residuals
# predictors <- resid_vals %>% dplyr::select(one_of(snps)) %>% as.matrix() 
targets    <- resid_vals %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()


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
  set.seed(10)
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
# raw_rep  <- raw_rep1[complete.cases(raw_rep1[,c(covar)]),] # N = 359
# raw_rep1 <- raw_rep1[,c(pred_all,outcomes,"FID.x")] # N = 363
raw_rep  <- raw_rep1[complete.cases(raw_rep1),] # N = 336

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
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

# create residual variables for all predictors and outcomes
resid_swe.input <- cbind(rep.predictors,rep.outcomes.s,raw_rep[,covar]) %>% as.data.frame
resid_swe <- lapply(outcomes, function(x) {
  lm(eval(substitute(i ~ age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                     list(i = as.name(x)))), data = resid_swe.input)
})

resid_swe_vals        <- lapply(resid_swe, function(x) residuals(x)) %>% as.data.frame()
names(resid_swe_vals) <- outcomes

#reset targets and predictions to residuals
# rep.predictors <- resid_swe_vals %>% dplyr::select(one_of(snps)) %>% as.matrix() 
rep.targets    <- resid_swe_vals %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

# apply models to swedish dataset
rf.pred.tr     <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr    <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))


##### ##### ##### ##### ##### 
##### Mixed effects models controlling for family ID
##### ##### ##### ##### ##### 

rf.pred  <- rf.pred.tr %>% as.data.frame
glm.pred <- glm.pred.tr %>% as.data.frame
mix <- cbind(rep.targets, rf.pred, glm.pred, raw_rep[,"FID.x"]) %>% as.data.frame

outcomes.rf <- NULL  # start empty list
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7]   <- outcomes            # renaming using the list of names
names(mix)[8:14]  <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22]    <- "fid"

outstats <- matrix(nrow=7,ncol=13)
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
  outstats[i,5] <- sqrt(mean(m1$residuals^2))
  outstats[i,6] <- mean(m1$residuals^2)
  outstats[i,7] <- summary(m2)$tTable[2,4]
  outstats[i,8] <- summary(m2)$tTable[2,5]
  outstats[i,9] <- r2m2[r2m2$Effect==paste0(outcomes[i],".glm"), 6]
  outstats[i,10] <- sqrt(mean(m2$residuals^2))
  outstats[i,11] <- mean(m2$residuals^2)
  outstats[i,12] <- summary(m1)$tTable[2,5]/2
  outstats[i,13] <- summary(m2)$tTable[2,5]/2
}

outstats        <- as.data.frame(outstats, stringsAsFactors = F)
outstats[,2:13]  <- lapply(outstats[,2:13], as.numeric)
names(outstats) <- c("var", "tval.rf", "pval.rf", "r2.rf", "rmse.rf", "mse.rf",
                     "tval.glm", "pval.glm", "r2.glm", "rmse.glm", "msf.rf", 
                     "pval.rf.onetail", "pval.glm.onetail")

# write.table(outstats, "revision2/covar-resid.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(outstats, "revision2/covar-resid-targets-10PCs.txt", col.names=T, row.names=F, sep="\t", quote=F)

save.image("revision2/covar-resid-targets-10PCs.RData")

# # CNP correlations and p-values
# glm.r2     <- unlist(lapply(glm.perfs, function(x) x$TrainRsquared))
# rf.r2      <- unlist(lapply(rf.perfs, function(x) x$TrainRsquared))
# 
# r2   <- rbind(glm.r2,unname(rf.r2))
# rval <- apply(r2,2,sqrt)
# tval <- apply(r2,2,function(x) sqrt((737*x)/(1-x)))
# pval <- apply(tval,2,function(x) 2*(1-abs(pt(x,df=737))))
# 


# general ability check 
# calculated from all tasks excepting model-specific task
# scaled, winsorized variables 

# all 6 tasks
g_all <- principal(rep.targets[,c(1:3,5:7)], nfactors=1, rotation="none",scores=T)
g_all # fit = .97, proportion var = .61
mix$g_all <- as.numeric(g_all$scores)

# leave out trails1
g_tr1 <- principal(rep.targets[,c(1:2,5:7)], nfactors=1, rotation="none",scores=T)
g_tr1 # fit = .98, proportion var = .72
mix$g_tr1 <- as.numeric(g_tr1$scores)

# leave out vr2
g_vr <- principal(rep.targets[,c(1:3,5:6)], nfactors=1, rotation="none",scores=T)
g_vr # fit = .96, proportion var = .60
mix$g_vr <- as.numeric(g_vr$scores)

# mixed effect models with g vars for rf
mx.tr1.g <- lme(g_tr1 ~ CRT_TIME1.rf, random = ~1 | fid, data = mix, na.action = na.omit)
mx.vr.g  <- lme(g_vr ~ VR2DR_TOTALRAW.rf, random = ~1 | fid, data = mix, na.action = na.omit)

# to load later
save.image(file=paste0(workDir,"/revision2/covar-resid-targets-10PCs-g.Rdata"))

#################### figures
load(paste0(workDir,"/revision2/covar-resid-targets-10PCs-g.Rdata"))
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
fig_trails1[4,4] <- "CILP2, GATAD2A, HAPLN4*"
fig_trails1[7,4] <- "NOSIP, PRR12, PRRG2*"
fig_trails1[13,4] <- "ANP32E, APH1A, C1orf51*"
fig_trails1[16,4] <- "MSL2, NCK1, PCCB*"

# visual reproduction - delayed recall
fig_vrdr           <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr           <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr)    <- c("SNP", "Importance")
fig_vrdr$Gene      <- proxy[match(fig_vrdr$SNP, proxy$proxy), 3]
fig_vrdr$SNP       <- factor(fig_vrdr$SNP,
                             levels = fig_vrdr$SNP[order(fig_vrdr$Importance)])
fig_vrdr$Gene2 <- fig_vrdr$Gene %>% as.character()
fig_vrdr[1,4] <- "CDC25C, CTNNA1, EGR1*"
fig_vrdr[3,4] <- "ABCB9, ARL6IP4, C12orf65*"
fig_vrdr[17,4] <- "AL049840.1, APOPT1, BAG5*"
fig_vrdr[20,4] <- "ANP32E, APH1A, C1orf51*"

# code by shared vs. task-specific genes
fig_trails1$shared         <- ifelse(fig_trails1$Gene %in% fig_vrdr[1:25,3],1,2) %>% as.factor()
fig_vrdr$shared            <- ifelse(fig_vrdr$Gene %in% fig_trails1[1:25,3],1,2) %>% as.factor()
levels(fig_trails1$shared) <- c("both tasks","task-specific")
levels(fig_vrdr$shared)    <- c("both tasks","task-specific")

# trails 1 figure
ggplot(fig_trails1[1:10,], aes(x=SNP,y=Importance)) + 
  geom_bar(fill="#CCCCCC",stat="identity",width=.85) +
  geom_text(stat="identity",y=49, hjust=0, aes(label=SNP), colour="black", size=5) +
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
  coord_flip(ylim = c(50,100))

# vis reprod figure
ggplot(fig_vrdr[1:10,], aes(x=SNP,y=Importance)) + 
  geom_bar(fill="#CCCCCC", stat="identity",width=.85) +
  geom_text(stat="identity",y=49, hjust=0, aes(label=SNP), colour="black", size=5) +
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
  coord_flip(ylim = c(50,100))



# supp figures
# correlations between predicted and observed values
ggplot(mix, aes(x=CRT_TIME1,y=CRT_TIME1.rf)) +    # 6in x 5in
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=18,face="bold")) +
  ggtitle("S5A) Processing Speed Random Forest") +
  ylab("\nTrails 1/A Predicted Values") +
  xlab("Residual Trails 1/A Performance (Scaled)")

ggplot(mix, aes(x=VR2DR_TOTALRAW,y=VR2DR_TOTALRAW.rf)) +    # 6in x 5in
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=18,face="bold")) +
  ggtitle("S5B) Visual Memory Random Forest") +
  ylab("\nVR II Predicted Values") +
  xlab("Residual VR II Performance (Scaled)")

