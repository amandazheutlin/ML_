# POLYGENIC SCORES
# CNP data
# use all the caucasians
# lots of outcome variables
# 108 PGC SNPs

mldir <- "/data/ML_genetics/Schizophrenia-Zheutlin"
setwd(mldir)

libs <- c("ggplot2", "dplyr", "nlme", "stringr")
invisible(lapply(libs, require, character.only = TRUE))

# load data
cnp <- read.table("dfs/GWASinput2-cnp.txt",header=T)
swe <- read.table("dfs/GWASinput2-swe.txt",header=T)

# add PCs
cnp.pc    <- read.table("/data/swe_gwas/ABZ/ML/cnp_MLsubs_pca.eigenvec", stringsAsFactors = F)
cnp.df.pc <- merge(cnp, cnp.pc, by.x = "ptid", by.y = "V1")

pc.names <- NULL
for (i in 1:20){pc.names[i] <- paste0("PC",i)}  
names(cnp.df.pc)[89:109] <- c("ptid2",pc.names)

swe.pc    <- read.table("/data/swe_gwas/ABZ/ML/swe_MLsubs2.mds", header=T, stringsAsFactors = F)
swe.df.pc <- merge(swe, swe.pc, by = "IID")


########################################### COG-BASED PRS
# pheno files of residual vars
outcomes <- names(cnp.df.pc)[c(5:7,10:12)]
cnp.df   <- cnp.df.pc[(complete.cases(cnp.df.pc$CRT_TIME1)),]
resid    <- lapply(outcomes, function(x) {
  lm(eval(substitute(i ~ age + gender + PC1 + PC2 + PC3 + PC4 + PC5, list(i = as.name(x)))),data = cnp.df)
})

resid_vals        <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
resid_vals        <- cbind(cnp.df[,c(1,1)],resid_vals)
names(resid_vals) <- c("FID","IID",outcomes)

for (i in 1:6){
  write.table(resid_vals[,c(1,2,i+2)],
              paste0("cogPRS/",names(resid_vals)[i+2],"-resid-PCs-cnp.txt"),
              col.names=T,row.names=F,quote=F,sep="\t")
}

#### run GWAS-cogvars-PCs.sh

# input files for score generation
assoc_files <- list.files(path = "cogPRS/.", pattern="77SNPassoc-PCs*.*assoc.linear")
assoc <- lapply(assoc_files,function(x) read.table(paste0("cogPRS/", x),header=T))
names(assoc) <- c("CRT","CVLT","DS","VOC","VR1","VR2")

for (i in 1:length(assoc)){
  write.table(assoc[[i]][,c(2,4,7)],
              paste0("cogPRS/77SNP-PCs-",names(assoc)[i],".txt"),
              col.names=T,row.names=F,quote=F,sep="\t")
}

#### run 77SNP-cog-scores-PCs.sh

# test scores in swedish sample
score_files_swe   <- list.files(path = "cogPRS/.", pattern="77SNP-PCs-swe.*\\.profile")
scores_swe        <- lapply(score_files_swe,function(x) read.table(paste0("cogPRS/", x),header=T))
names(scores_swe) <- c("CRT","CVLT","DS","VOC","VR1","VR2")

names(swe.df.pc)[c(8:10,13:15)] <- outcomes
swe.df      <- swe.df.pc[!apply(swe.df.pc[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),]
swe.df$CRT  <- scores_swe[["CRT"]][match(swe.df$IID,scores_swe[["CRT"]][,2]),6]
swe.df$CVLT <- scores_swe[["CVLT"]][match(swe.df$IID,scores_swe[["CVLT"]][,2]),6]
swe.df$DS   <- scores_swe[["DS"]][match(swe.df$IID,scores_swe[["DS"]][,2]),6]
swe.df$VOC  <- scores_swe[["VOC"]][match(swe.df$IID,scores_swe[["VOC"]][,2]),6]
swe.df$VR1  <- scores_swe[["VR1"]][match(swe.df$IID,scores_swe[["VR1"]][,2]),6]
swe.df$VR2  <- scores_swe[["VR2"]][match(swe.df$IID,scores_swe[["VR2"]][,2]),6]

swe.lm1 <- lme(CRT_TIME1 ~ CRT + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)
swe.lm2 <- lme(CVLT_TOTCOR ~ CVLT + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)
swe.lm3 <- lme(DS_TOTALRAW ~ DS + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)
swe.lm4 <- lme(VOC_TOTALRAW ~ VOC + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)
swe.lm5 <- lme(VR1IR_TOTALRAW ~ VR1 + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)
swe.lm6 <- lme(VR2DR_TOTALRAW ~ VR2 + age + sex + C1 + C2 + C3 + C4 + C5,
               random=~1|FID.x, data = swe.df, na.action=na.omit)


########################################### 77 SZ RISK SCORES
# input file for scoring 
pgc          <- read.table("scz2.rep.128.txt",header=T, stringsAsFactors = F)
snp_tbl      <- read.table("index-proxy-table_allsnps.txt",header=T, stringsAsFactors = F)
input        <- pgc[pgc$snpid %in% snp_tbl$index, c("snpid","a1a2","or")] # only 77 SNPs
input$A1     <- sapply(input[,2], function(x) unlist(strsplit(x,""))[1]) # A1
input        <- input[order(input$snpid),]
input[1:2,4] <- c("I2","I12") # manually fix random A1 nonsense
input$lnOR   <- log(input$or) # natual log of OR

# proxy snps
input$prx    <- snp_tbl[match(input$snpid,snp_tbl$index),2] # proxy snps
input$prxswe <- input$prx # proxy snps
snp_tbl[snp_tbl$prxswe == "rs2535627",3] <- "rs4481150"
snp_tbl[snp_tbl$prxswe == "rs4388249",3] <- "rs12656763"
snp_tbl[snp_tbl$prxswe == "rs2905426",3] <- "rs2905425"

write.table(input[,c("prx","A1","lnOR")],"prs_77snp_input.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

write.table(input[,c("prxswe","A1","lnOR")],"prs_77snp_input_swe.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# test scores
cnp_prs <- read.table("SZ-PRS-77SNP-cnp.profile",header=T)
swe_prs <- read.table("SZ-PRS-77SNP-swe.profile",header=T)

cnp$prs <- cnp_prs[match(cnp$ptid,cnp_prs$IID),6]
swe$prs <- swe_prs[match(swe$IID,swe_prs$IID),6]

vars             <- colnames(cnp)[5:12]
names(swe)[8:15] <- vars
names(cnp)[3]    <- "sex"

cnp.pc    <- read.table("/data/swe_gwas/ABZ/ML/cnp_MLsubs_pca.eigenvec", stringsAsFactors = F)
cnp.df.pc <- merge(cnp, cnp.pc, by.x = "ptid", by.y = "V1")

pc.names <- NULL
for (i in 1:20){pc.names[i] <- paste0("PC",i)}  
names(cnp.df.pc)[91:111] <- c("ptid2",pc.names)

swe.pc    <- read.table("/data/swe_gwas/ABZ/ML/swe_MLsubs2.mds", header=T, stringsAsFactors = F)
swe.df.pc <- merge(swe, swe.pc, by = "IID")

covar     <- names(cnp.df.pc)[c(2,3,92:96)]
names(swe.df.pc)[c(6,7,96:100)] <- covar

swe.df1 <- swe.df.pc[!apply(swe.df.pc[,c(vars)], 1, function(x) (mean(is.na(x))==1)),]
swe.df  <- swe.df1[complete.cases(swe.df1[,c(covar)]),] # N = 359

prs_corrs <- matrix(nrow=8,ncol=7)
for (i in 1:length(vars)){
  f <- formula(paste(vars[i],"prs + age + sex + PC1 + PC2 + PC3 + PC4 + PC5",sep="~"))
  m1 <- lm(f, data=cnp.df.pc)
  m2 <- lme(fixed=f, random = ~1|FID.x, data = swe.df, na.action=na.omit)
  r2m1 = r2beta(m1,method='sgv')
  r2m2 = r2beta(m2,method='sgv')
  prs_corrs[i,1] <- vars[i]
  prs_corrs[i,2] <- summary(m1)$coefficients[2,3]
  prs_corrs[i,3] <- r2m1[r2m1$Effect=="prs", 6]
  prs_corrs[i,4] <- summary(m1)$coefficients[2,4] / 2 # one-tailed p-values
  prs_corrs[i,5] <- summary(m2)$tTable[2,4]
  prs_corrs[i,6] <- r2m1[r2m1$Effect=="prs", 6]
  prs_corrs[i,7] <- summary(m2)$tTable[2,5] / 2 # one-tailed p-values
}

prs_corrs         <- data.frame(prs_corrs, stringsAsFactors=F) 
prs_corrs[,2:7]   <- lapply(prs_corrs[,2:7], as.numeric) %>% as.data.frame
names(prs_corrs)  <- c("cogvar","cnp_t","cnp_r2","cnp_p_onetail","swe_t","swe_r2","swe_p_onetail")

prs_corrs_use         <- prs_corrs[c(1:3,6:8),]
prs_corrs_use$cnp_FDR <- p.adjust(prs_corrs_use$cnp_p,method="fdr")
prs_corrs_use$swe_FDR <- p.adjust(prs_corrs_use$swe_p,method="fdr")

write.table(prs_corrs_use,"prs77_corrs_PCs.txt",col.names=T,
            row.names=F,quote=F,sep="\t")

########################################### GWAS SZ RISK SCORES
scoredir <- "~/Google Drive/ABZ/Yale/CNP_Genetics/CVLTscores_imputed/scores/"
cnpsz <- read.table(paste0(scoredir,"SZrisk-pT5-cnpsubs-callgenos.profile"),header=T)
swesz <- read.table(paste0(scoredir,"SZrisk-pT5-omnisubs-callgenos.profile"),header=T)

cnp$SZprs <- cnpsz[match(cnp$ptid,cnpsz$IID),6]
swe$SZprs <- swesz[match(swe$IID,swesz$IID),6]

vars <- colnames(cnp)[5:12]
names(swe)[8:15] <- vars
names(cnp)[3] <- "sex"

prs_corrs <- matrix(nrow=8,ncol=7)
for (i in 1:length(vars)){
  f <- formula(paste(vars[i],"SZprs + age + sex",sep="~"))
  m1 <- lm(f, data=cnp)
  m2 <- lme(fixed=f, random = ~1|FID, data = swe, na.action=na.omit)
  cnptval <- summary(m1)$coefficients[2,3]
  swetval <- summary(m2)$tTable[2,4]
  prs_corrs[i,1] <- vars[i]
  prs_corrs[i,2] <- cnptval
  prs_corrs[i,3] <- sqrt(cnptval^2/(cnptval^2+summary(m1)$df[2]))
  prs_corrs[i,4] <- summary(m1)$coefficients[2,4]
  prs_corrs[i,5] <- swetval
  prs_corrs[i,6] <- sqrt(swetval^2/(swetval^2+summary(m2)$tTable[2,3]))
  prs_corrs[i,7] <- summary(m2)$tTable[2,5]
}

prs_corrs <- as.data.frame(prs_corrs)
names(prs_corrs) <- c("cogvar","cnp_t","cnp_corr","cnp_p","swe_t","swe_corr","swe_p")

write.table(prs_corrs,"SZprs_pT5_corrs.txt",col.names=T,
            row.names=F,quote=F,sep="\t")


