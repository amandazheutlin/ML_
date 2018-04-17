# POLYGENIC SCORES
# CNP data
# use all the caucasians
# lots of outcome variables
# 108 PGC SNPs

mldir <- "~/Google Drive/ABZ/Yale/ML_PGC-SZ/"
pgcdir <- "~/Google Drive/ABZ/Yale/PGC-Results/2014/"
setwd(mldir)

# load data
cnp <- read.table("GWASinput2-cnp.txt",header=T)
swe <- read.table("GWASinput2-swe.txt",header=T)

# GWAS hit - proxy SNP table (same for CNP and swedish)
# proxy <- read.table("index-proxy-table.txt",header=T)
# proxy[,1:2] <- data.frame(lapply(proxy[,1:2],as.character),stringsAsFactors=F)
# 
# # SNPs (77) - CNP
# cnp_snps <- names(cnp)[!(names(cnp)%in%proxy$PROXY)] # index snps
# cnp_snps_idx <- c(cnp_snps[13:71],proxy[,1]) # index snps + proxy index snps
# cnp_snps_prx <- c(cnp_snps[13:71],proxy[,2]) # index snps + proxy snps
# snp_tbl <- as.data.frame(cnp_snps_idx)
# snp_tbl[,2] <- cnp_snps_prx
# names(snp_tbl) <- c("index","proxy")
# 
# write.table(snp_tbl,"index-proxy-table_allsnps.txt",quote=F,
#             col.names=T,row.names=F,sep="\t")

# SNPs (77) - swedish
# same SNPs except 3
# index_snp  proxy_snp
# rs2535627  rs4481150
# rs4388249	 rs12656763
# rs2905426	 rs2905425
# snp_tbl$proxy_swe <- cnp_snps_prx
# snp_tbl[snp_tbl$proxy_swe=="rs2535627",3] <- "rs4481150"
# snp_tbl[snp_tbl$proxy_swe=="rs4388249",3] <- "rs12656763"
# snp_tbl[snp_tbl$proxy_swe=="rs2905426",3] <- "rs2905425"

########################################### COG-BASED PRS
# pheno files of residual vars
library(dplyr)
outcomes <- names(cnp)[c(5:7,10:12)]
cnp.df   <- cnp[(complete.cases(cnp$CRT_TIME1)),]
resid    <- lapply(outcomes, function(x) {
  lm(eval(substitute(i ~ age + gender, list(i = as.name(x)))),data = cnp.df)
})

resid_vals <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
resid_vals <- cbind(cnp.df[,c(1,1)],resid_vals)
names(resid_vals) <- c("FID","IID",outcomes)

for (i in 1:6){
  write.table(resid_vals[,c(1,2,i+2)],
              paste0(names(resid_vals)[i+2],"-resid-cnp.txt"),
              col.names=T,row.names=F,quote=F,sep="\t")
}

# input files for score generation
assoc_files <- list.files(pattern="*.assoc.linear")
assoc <- lapply(assoc_files,function(x) read.table(x,header=T))
names(assoc) <- c("CRT","CVLT","DS","VOC","VR1","VR2")

for (i in 1:length(assoc)){
  write.table(assoc[[i]][,c(2,4,7)],
              paste0("77SNP-",names(assoc)[i],".txt"),
              col.names=T,row.names=F,quote=F,sep="\t")
}

# test scores in swedish sample
score_files_swe <- list.files(pattern="77SNP-swe.*\\.profile")
scores_swe <- lapply(score_files_swe,function(x) read.table(x,header=T))
names(scores_swe) <- c("CRT","CVLT","DS","VOC","VR1","VR2")

names(swe)[c(8:10,13:15)] <- outcomes
swe.df <- swe[!apply(swe[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),]
swe.df$CRT <- scores_swe[["CRT"]][match(swe.df$IID,scores_swe[["CRT"]][,2]),6]
swe.df$CVLT <- scores_swe[["CVLT"]][match(swe.df$IID,scores_swe[["CVLT"]][,2]),6]
swe.df$DS <- scores_swe[["DS"]][match(swe.df$IID,scores_swe[["DS"]][,2]),6]
swe.df$VOC <- scores_swe[["VOC"]][match(swe.df$IID,scores_swe[["VOC"]][,2]),6]
swe.df$VR1 <- scores_swe[["VR1"]][match(swe.df$IID,scores_swe[["VR1"]][,2]),6]
swe.df$VR2 <- scores_swe[["VR2"]][match(swe.df$IID,scores_swe[["VR2"]][,2]),6]

library(nlme)
swe.lm1 <- lme(CRT_TIME1 ~ CRT + age + sex,
      random=~1|FID, data = swe.df, na.action=na.omit)
swe.lm2 <- lme(CVLT_TOTCOR ~ CVLT + age + sex,
               random=~1|FID, data = swe.df, na.action=na.omit)
swe.lm3 <- lme(DS_TOTALRAW ~ DS + age + sex,
               random=~1|FID, data = swe.df, na.action=na.omit)
swe.lm4 <- lme(VOC_TOTALRAW ~ VOC + age + sex,
               random=~1|FID, data = swe.df, na.action=na.omit)
swe.lm5 <- lme(VR1IR_TOTALRAW ~ VR1 + age + sex,
               random=~1|FID, data = swe.df, na.action=na.omit)
swe.lm6 <- lme(VR2DR_TOTALRAW ~ VR2 + age + sex,
               random=~1|FID, data = swe.df, na.action=na.omit)


########################################### 77 SZ RISK SCORES
# input file for scoring 
pgc <- read.table(paste0(pgcdir,"scz2.rep.128.txt"),header=T)
input <- pgc[pgc$snpid%in%cnp_snps_idx,c("snpid","a1a2","or")] # only 77 SNPs
input[,1:2] <- data.frame(lapply(input[,1:2],as.character),stringsAsFactors=F)
input$A1 <- sapply(input[,2], function(x) unlist(strsplit(x,""))[1]) # A1
input <- input[order(input$snpid),]
input[1:2,4] <- c("I2","I12") # manually fix random A1 nonsense
input$proxyid <- snp_tbl[match(input$snpid,snp_tbl$index),2] # proxy snps
input$proxyid_swe <- snp_tbl[match(input$snpid,snp_tbl$index),3] # proxy snps

write.table(input[,c("proxyid","A1","or")],"prs_77snp_input.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

write.table(input[,c("proxyid_swe","A1","or")],"prs_77snp_input_swe.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# test scores
cnp_prs <- read.table("cnp_77prs.profile",header=T)
swe_prs <- read.table("swe_77prs.profile",header=T)

cnp$prs <- cnp_prs[match(cnp$ptid,cnp_prs$IID),6]
swe$prs <- swe_prs[match(swe$IID,swe_prs$IID),6]

vars <- colnames(cnp)[5:12]
names(swe)[8:15] <- vars
names(cnp)[3] <- "sex"

swe.df <- swe[!apply(swe[,c(vars)], 1, function(x) (mean(is.na(x))==1)),]

library(nlme)
prs_corrs <- matrix(,nrow=8,ncol=7)
for (i in 1:length(vars)){
  f <- formula(paste(vars[i],"prs + age + sex",sep="~"))
  m1 <- lm(f, data=cnp)
  m2 <- lme(fixed=f, random = ~1|FID, data = swe.df, na.action=na.omit)
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

write.table(prs_corrs,"prs77_corrs.txt",col.names=T,
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

library(nlme)
prs_corrs <- matrix(,nrow=8,ncol=7)
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

