# demo stats for ML paper
# take out age, sex, and diagnosis

mldir <- "/data/ML_genetics/Schizophrenia-Zheutlin"
setwd(mldir)

library(psych); library(dplyr); library(nlme)

# load data
swe <- read.table("dfs/GWASinput2-swe.txt", header=TRUE, as.is=TRUE)
cnp <- read.table("dfs/GWASinput2-cnp.txt", header=TRUE, as.is=TRUE)

# exclude symbol / spatial span
swe <- swe[,c(1:11,13:92)]
cnp <- cnp[,c(1:8,10:89)]

# only the people we used 
# CNP (N = 739)
# swedish (N = 336)
cnp.subs <- read.table("revision2/cnp.subs.txt",header=F,stringsAsFactors = F) 
swe.subs <- read.table("revision2/swe.subs.txt",header=F,stringsAsFactors = F)

outcomes <- names(cnp)[c(5:7,9:11)]
names(swe)[c(8:10,12:14)] <- outcomes

cnp.df <- cnp[cnp$ptid %in% cnp.subs$V1,]
swe.df <- swe[swe$IID %in% swe.subs$V1,]

# scale the swedish variables to the CNP ranges
swe.df$VOC_TOTALRAW <- swe.df$VOC_TOTALRAW * 1.25
swe.df$DS_TOTALRAW <- swe.df$DS_TOTALRAW * 1.6
swe.df$VR1IR_TOTALRAW <- swe.df$VR1IR_TOTALRAW * .398
swe.df$VR2DR_TOTALRAW <- swe.df$VR2DR_TOTALRAW * .398

# CNP
describe(cnp.df$age)
table(cnp.df$gender)
table(cnp.df$Group)
describe(cnp.df[,c(outcomes)])

# swedish (mixed effect)
describe(swe.df$age)
table(swe.df$sex)
table(swe.df$dx)
describe(swe.df[,c(outcomes)])

# merge samples
cnp.merge             <- cnp.df[,c("ptid","ptid","age","gender","Group",outcomes)]
names(cnp.merge)[1:5] <- c("FID","IID","age","sex","dx")
cnp.merge$dx          <- as.factor(cnp.merge$dx)
levels(cnp.merge$dx)  <- c(15,3,7,1) # match levels to swedish (adhd is new = 15)
cnp.merge$sample <- 1

swe.merge        <- swe.df[,c("FID","IID","age","sex","dx",outcomes)]
swe.merge$dx     <- as.factor(swe.merge$dx)
swe.merge$sample <- 2

merge        <- rbind(cnp.merge,swe.merge)
merge$sex    <- as.factor(merge$sex)
merge$sample <- as.factor(merge$sample)

require(car)
merge$patient <- recode(merge$dx,"1=1; 2=2; 3=1; 4:7=2; 9=2; 10=2; 15=1")

# test for demo differences
chisq.test(merge$sex,merge$sample)
chisq.test(merge$patient,merge$sample)
age <- lme(age ~ sample, random = ~1|FID, data = merge)
summary(age)

# test for differences in cog accounting for age and dx
# (not sex bc no sample differences)
CVLT <- lme(CVLT_TOTCOR ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)
VOC  <- lme(VOC_TOTALRAW ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)
CRT  <- lme(CRT_TIME1 ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)
DS   <- lme(DS_TOTALRAW ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)
VR1  <- lme(VR1IR_TOTALRAW ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)
VR2  <- lme(VR2DR_TOTALRAW ~ sample + age + patient, random = ~1|FID, data = merge, na.action = na.omit)

summary(CVLT)
summary(VOC)
summary(CRT)
summary(DS)
summary(VR1)
summary(VR2)
