# Schizophrenia-linked genetic influences on cognition
# Authors: AMC & ABZ, September 2017

#load things
mldir <- "/data/ML_genetics/Schizophrenia-Zheutlin"
setwd(mldir)

libs <- c("ggplot2", "RColorBrewer", "psych", "dplyr", "stringr")
invisible(lapply(libs, require, character.only = TRUE))


# figure of importance variables
trails <- read.table("revision2/trails_imp.txt", header=T, stringsAsFactors = F)
vr     <- read.table("revision2/vr_imp.txt", header=T, stringsAsFactors = F)
cvlt   <- read.table("cogPRS/77SNP-10PCs-CVLT.txt", header=T, stringsAsFactors = F)
ds     <- read.table("cogPRS/77SNP-10PCs-DS.txt", header=T, stringsAsFactors = F)

# align 'ranking' of SNP based on trails
trails$trailsimp <- seq(1,77,by=1) %>% as.factor()
vr$trailsimp     <- trails[match(vr$SNP, trails$SNP), "trailsimp"] %>% as.factor()
cvlt$trailsimp   <- trails[match(cvlt$SNP, trails$SNP), "trailsimp"] %>% as.factor()
ds$trailsimp     <- trails[match(ds$SNP, trails$SNP), "trailsimp"] %>% as.factor()

# make betas positive
cvlt$BETA.pos    <- abs(cvlt$BETA)
ds$BETA.pos      <- abs(ds$BETA)

# cosmetics
trails$SNP       <- factor(trails$SNP, levels = trails$SNP[order(trails$Importance)])
vr$SNP           <- factor(vr$SNP, levels = vr$SNP[order(vr$Importance)])
cvlt$SNP         <- factor(cvlt$SNP, levels = cvlt$SNP[order(cvlt$BETA.pos)])
ds$SNP           <- factor(ds$SNP, levels = ds$SNP[order(ds$BETA.pos)])


# trails 1 figure
trails.plot <- ggplot(trails, aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=trailsimp),stat="identity",width=.85) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(77)) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text.x = element_text(size=13,color="black"),
        axis.text.y = element_text(size=8,color="black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=16,face="bold")) +
  ggtitle("\n") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(0,100))

# vis reprod figure
vr.plot <- ggplot(vr, aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=trailsimp), stat="identity",width=.85) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(77)) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text.x = element_text(size=13,color="black"),
        axis.text.y = element_text(size=8,color="black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=16,face="bold")) +
  ggtitle("\n") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(0,100))

# cvlt figure
cvlt.plot <- ggplot(cvlt, aes(x=SNP,y=BETA.pos)) + 
  geom_bar(aes(fill=trailsimp), stat="identity",width=.85) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(77)) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text.x = element_text(size=13,color="black"),
        axis.text.y = element_text(size=8,color="black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=16,face="bold")) +
  ggtitle("\n") +
  ylab("\nBeta Weight") +
  xlab("") +
  coord_flip(ylim = c(0,.15))

# ds figure
ds.plot <- ggplot(ds, aes(x=SNP,y=BETA.pos)) + 
  geom_bar(aes(fill=trailsimp), stat="identity",width=.85) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(77)) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text.x = element_text(size=13,color="black"),
        axis.text.y = element_text(size=8,color="black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.4,vjust=1.8,size=16,face="bold")) +
  ggtitle("\n") +
  ylab("\nBeta Weight") +
  xlab("") +
  coord_flip(ylim = c(0,.18))


multiplot(trails.plot, vr.plot, cvlt.plot, ds.plot, cols=4) # 14in x 10in








# graph as facets
# names(ds)[5]   <- "Importance"
# names(cvlt)[5] <- "Importance"
# trails$task    <- "TRAILS 1/A"
# vr$task        <- "VR II"
# cvlt$task      <- "CVLT"
# ds$task        <- "Digit Symbol"
# 
# 
# combo <- rbind(trails,vr,cvlt[,c("SNP","Importance","trailsimp","task")],
#                ds[,c("SNP","Importance","trailsimp","task")])
