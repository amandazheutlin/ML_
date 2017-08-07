# MACHINE LEARNING
# CNP data
# use all the caucasians
# lots of outcome variables
# 108 PGC SNPs

setwd("~/Google Drive/ABZ/Yale/ML_PGC-SZ")
############################################ CNP
# load data
data <- read.table("~/Google Drive/ABZ/Yale/CNP_Genetics/CNP_NP/mm_data.txt",header=T)
cauc <- subset(data,RACE_MAIN==5 & ethnicity==2,select=c("ptid","age","gender","Group",
                                                         "CVLT_TOTCOR","DS_TOTALRAW",
                                                         "MR_TOTALRAW","VOC_TOTALRAW",
                                                         "CRT_TIME1","CRT_TIME2",
                                                         "rk_dprime","scap_dprime"))

# missing data
# too many missing from rk_dprime and scap_dprime
na_count <- sapply(cauc, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
cauc <- cauc[,1:10]

# add genotypes
# 60 out of 108 SNPs (56%)
snps <- read.table("108-PGC-genos-cnp.txt",header=T)
cauc_snps <- merge(cauc,snps,by.x="ptid",by.y="IID")
cauc_snps <- cauc_snps[,c(1:10,13:72)]

# write df
write.table(cauc_snps,"GWAS-snps-cnp-input.txt",col.names=T,row.names=F,
            quote=F,sep="\t")

### new df with proxy snps and diff cog variables
# only use cog vars in swedish sample (8)
# CVLT, vocab, trails 1(A) and 2(B)
# spatial span, digit span
# immediate visual reproduction, delayed VR
cauc2 <- subset(data,RACE_MAIN==5 & ethnicity==2,
                select=c("ptid","age","gender","Group","CVLT_TOTCOR",
                         "VOC_TOTALRAW","CRT_TIME1","CRT_TIME2","SSP_TOTALRAW",
                         "DS_TOTALRAW","VR1IR_TOTALRAW","VR2DR_TOTALRAW"))

# missing data
# one missing from spatial span and trails1
na_count <- sapply(cnp, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

# add genotypes
# 60 out of 108 SNPs (56%) 
# + 18 proxy SNPs 
# - 1 SNP not in swedish data (rs35518360)
# --> 77 SNPs out of 108 (71%)
snps <- read.table("108-PGC-genos-cnp.txt",header=T)
proxy <- read.table("108-proxysnps-cnp.txt",header=T)
cauc2_snps <- merge(cauc2,snps,by.x="ptid",by.y="IID")
cauc2_snps <- merge(cauc2_snps,proxy,by.x="ptid",by.y="IID")
cauc2_snps <- cauc2_snps[,c(1:12,15:27,29:74,76:93)] #exclude outlier

# write df2
write.table(cauc2_snps,"GWASinput2-cnp.txt",col.names=T,row.names=F,
            quote=F,sep="\t")

############################################ Swedish
# Trails A & B (Trail_A_test, Trail_B_test) and 60 snps
# need 3 proxy snps for CNP overlap: rs2535627, rs4388249, rs2905426
# index_snp  proxy_snp
# rs2535627	rs4481150
# rs4388249	rs12656763
# rs2905426	rs2905425

# np
data <- read.table("~/Google Drive/ABZ/Yale/RNA_SwedishTwins/swedishNP_PCA.txt",header=T)
data <- data[,c(1:2,7:9,14:15,24:25,20,17)]

# snps
snps_swe <- read.table("108-PGC-genos.txt",header=T)
data_snps <- merge(data,snps_swe,by="IID")
data_snps <- data_snps[,c(1:11,14:81)]

# df
df <- data_snps[,names(data_snps)%in%names(snps)]
df <- cbind(data_snps[,1:11],df[,2:57])
names(df)[2:3] <- c("FID","dx")
names(df)[7:9] <- c("sex","trailsA","trailsB")

# proxy snps
proxy <- read.table("proxysnps.txt",header=T)
df <- merge(df,proxy[,1:6],by="IID")
df <- df[,c(1:67,69:72)]

# missing data
na_count <- sapply(df, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

# write df
write.table(df,"GWAS-snps-swe-input_morecog.txt",col.names=T,row.names=F,
            quote=F,sep="\t")

### new df with proxy snps and diff cog variables
# cog vars from CNP (8)
# CVLT, vocab, trails 1(A) and 2(B)
# spatial span (wms_ss_t), digit span (waisiii3)
# immediate visual reproduction, delayed VR
data <- data[,c(1:2,7:9,14:15,20,17,24:25,21:22,19)]
names(data)[c(3,10:14)] <- c("dx","trailsA","trailsB","sspan","dspan","VRIR")
np <- read.xls("~/Google Drive/ABZ/Yale/SwedishData/SwedishData_All/swedish_Form_NP.xls",sheet=1,header=T)
data$VRDR <- np[match(data$IID,np$PTID),29]

# missing data
# ~45 people missing from all the cog data (41-54)
na_count <- sapply(data, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

# snps (56) + proxy snps (3) + more proxy snps (18) == 77 SNPs
snps_swe <- read.table("108-PGC-genos.txt",header=T)
data_snps <- merge(data,snps_swe,by="IID")
data_snps <- data_snps[,c(1:15,18:85)]

# overlapping CNP snps
df <- data_snps[,names(data_snps)%in%names(snps)]
df <- cbind(data_snps[,1:15],df[,2:57])
names(df)[c(2,7)] <- c("FID","sex")

# proxy swe snps
proxy <- read.table("proxysnps.txt",header=T)
df <- merge(df,proxy[,2:5],by="IID")

# proxy cnp/swe snps
proxy2 <- read.table("108-proxysnps-swe.txt",header=T)
df <- merge(df,proxy2[,2:20],by="IID")

# write df2
write.table(df,"GWASinput2-swe.txt",col.names=T,row.names=F,
            quote=F,sep="\t")

####### general ability variable
library(psych)
load("~/Downloads/amanda.results_2.Rdata")

# calculated from all 8 tasks
# use scaled variables 
# (they are also winsorized, trim = .01)
rep.outcomes.s <- as.data.frame(rep.outcomes.s)
g_all <- principal(rep.outcomes.s, nfactors=1, rotation="none",scores=T)
g_all # fit = .65, proportion var = .39
mix$g_all <- as.numeric(g_all$scores)

# leave out vocab OR trails1 OR VRV2
g_voc <- principal(rep.outcomes.s[,c(1,3:8)], nfactors=1, rotation="none",scores=T)
g_voc # fit = .55, proportion var = .36
mix$g_voc <- as.numeric(g_voc$scores)

g_tr1 <- principal(rep.outcomes.s[,c(1:2,4:8)], nfactors=1, rotation="none",scores=T)
g_tr1 # fit = .66, proportion var = .42
mix$g_tr1 <- as.numeric(g_tr1$scores)

g_vr <- principal(rep.outcomes.s[,c(1:7)], nfactors=1, rotation="none",scores=T)
g_vr # fit = .67, proportion var = .43
mix$g_vr <- as.numeric(g_vr$scores)

# mixed effect models with g vars
mx.tr1.g <- lme(g_tr1 ~ pred.tr1, random = ~1 | fid, data = mix, na.action = na.omit)
mx.voc.g <- lme(g_voc ~ pred.voc, random = ~1 |fid, data = mix, na.action=na.omit)
mx.vr.g <- lme(g_vr ~ pred.vr, random = ~1 | fid, data = mix, na.action = na.omit)

# mixed effect models testing outcome x prediction
mx.tr1 <- lme(out.tr1.log ~ pred.tr1, random = ~1 | fid, data = mix, na.action = na.omit)
mx.voc <- lme(outvoc ~ pred.voc, random = ~1 |fid, data = mix, na.action=na.omit)
mx.vr <- lme(outvr ~ pred.vr, random = ~1 | fid, data = mix, na.action = na.omit)



## Figure 1
# important variables annotated by 
# gene name and function

mldir <- "~/Google Drive/ABZ/Yale/ML_PGC-SZ/"
pgcdir <- "~/Google Drive/ABZ/Yale/PGC-Results/2014/"
setwd(mldir)

# snp-gene table
code <- read.csv(paste0(pgcdir,"scz2.anneal.108.genes.csv"))
proxy <- read.table(paste0(mldir,"index-proxy-table_allsnps.txt"),header=T)
proxy$gene <- code[match(proxy$index,code$bestsnp),8]

# vocabulary
fig_vocab <- caret::varImp(rf.models$VOC_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vocab <- cbind(rownames(fig_vocab), fig_vocab) %>% arrange(-Overall)
names(fig_vocab) <- c("SNP", "Importance")
fig_vocab$Gene <- proxy[match(fig_vocab$SNP,proxy$proxy),3]
fig_vocab$SNP <- factor(fig_vocab$SNP,levels=fig_vocab$SNP[order(fig_vocab$Importance)])

# really, really sketchy interpretation of genes
# not the best
vochits <- as.character(fig_vocab[1:20,3])
vochits <- strsplit(vochits,", ") %>% unlist() %>% as.data.frame
vochits[,2] <- vochits[,1] %in% ions[,1]
fig_vocab$ion <- c(2,2,1,1,2,1,1,1,1,1,2,1,1,2,2,2,2,2,1,2,rep(2,57)) %>% as.factor() # only care about top 20
levels(fig_vocab$ion) <- c("ion binding", "not ion binding")

library(ggplot2)
ggplot(fig_vocab[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=ion),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_vocab[1:20,3]),width=30)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=12,color="black"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip()

#ggsave("figure_GWAShits.eps", dpi=600)

# trails1
fig_trails1 <- caret::varImp(rf.models$CRT_TIME1, scale=TRUE)[[1]] %>% as.data.frame
fig_trails1 <- cbind(rownames(fig_trails1), fig_trails1) %>% arrange(-Overall)
names(fig_trails1) <- c("SNP", "Importance")
fig_trails1$Gene <- proxy[match(fig_trails1$SNP,proxy$proxy),3]
fig_trails1$SNP <- factor(fig_trails1$SNP,levels=fig_trails1$SNP[order(fig_trails1$Importance)])

# sketch interpretation
tr1hits <- as.character(fig_trails1[1:20,3])
tr1hits <- strsplit(tr1hits,", ") %>% unlist() %>% as.data.frame
tr1hits[,2] <- tr1hits[,1] %in% ions[,1]
tr1hits
fig_trails1$ion <- c(2,1,2,1,2,2,1,1,2,1,2,2,1,1,1,1,2,1,2,1,rep(2,57)) %>% as.factor() # only care about top 20
levels(fig_trails1$ion) <- c("ion binding", "not ion binding")

ggplot(fig_trails1[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=ion),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_trails1[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=12,color="black"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip()

# visual reproduction - delayed recall
fig_vrdr <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr) <- c("SNP", "Importance")
fig_vrdr$Gene <- proxy[match(fig_vrdr$SNP,proxy$proxy),3]
fig_vrdr$SNP <- factor(fig_vrdr$SNP,levels=fig_vrdr$SNP[order(fig_vrdr$Importance)])

# sketch interpretation
vrdrhits <- as.character(fig_vrdr[1:20,3])
vrdrhits <- strsplit(vrdrhits,", ") %>% unlist() %>% as.data.frame
vrdrhits[,2] <- vrdrhits[,1] %in% ions[,1]
vrdrhits
fig_vrdr$ion <- c(2,1,1,1,1,2,1,1,2,2,2,1,1,1,2,1,2,1,2,1,rep(2,57)) %>% as.factor() # only care about top 20
levels(fig_vrdr$ion) <- c("ion binding", "not ion binding")

ggplot(fig_vrdr[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=ion),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_vrdr[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=10,color="black"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip()



#figData %>%
#   dplyr::arrange(-Importance) %>%
#   dplyr::slice(1:20) %>%
#   ggplot() +
#   geom_bar(aes(x= gene, y=Importance), stat = "identity") +
#   coord_flip()
