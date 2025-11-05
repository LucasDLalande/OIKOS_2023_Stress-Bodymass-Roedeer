############################# ANALYSE STRESS-MASS ##############################
#################### EARLY-LATE RELATION STRESS X BODY MASS ####################
################# A./ PREDICTIONS POIDS MOYENS-----------------#################
################# B./ FIGURES----------------------------------#################
################################# 20/06/2022 ###################################
rm(list=ls())

library(lme4)
library(MuMIn)
library(ggplot2)
library(rptR)
library(lubridate)
library(dplyr)
library(ggthemes)
library(VGAM) #pour AICcmodavg
library(AICcmodavg) #pour les predicts
library(gameofthrones)
library(jtools) # pour la repr?sentation graphique 
library(interactions) # pour la repr?sentation graphique
library(visreg) # pour visualiser une relation (graph)
library(grid)
library(scales)
library(interplot) # visualiser l'interaction entre 2 variables quantitatives
library(gridExtra)
library(ggsci)

{
  data <- read.csv("dataset/dataset-20220516.csv", h=T, sep=";", dec=".")
  data <- data[!is.na(data$ageannee),]
  data$age <- data$ageannee
  data <- data[!is.na(data$masse),]
  data <- data[!is.na(data$sexe),]
  data$pop <- as.factor(data$pop)
  data$numind <- as.factor(data$numind)
  data$sexe <- as.factor(data$sexe)
  data$cohorte <- as.factor(data$cohorte)
  data <- data[!is.na(data$cohorte),]
  data <- data[!is.na(data$qualite_cohorte),]
  data$idkit <- as.factor(data$idkit)
  data$id_bioch <- as.factor(data$id_bioch)
  data$date <- dmy(data$date)
  
  # FGM of Priority = 3 --> NA
  data <- data %>%
    mutate((replace(FGMngg, FGM_Priority == "3", NA)))
  data <- data[,-18]
  names(data)[19] <- "FGMngg"
  
  rownames(data) <- 1:nrow(data)
  str(data)
  
  
  # Keep at least 3 obs per age class (to calculate averaged body mass per age/pop/sex)
  tapply(data$ageannee, data$pop:data$sexe, table)
  
  dataMCH <- data[data$pop=="C" & data$sexe=="M",]
  dataMCH <- dataMCH[dataMCH$ageannee<=10,]
  dataFCH <- data[data$pop=="C" & data$sexe=="F",]
  dataFCH <- dataFCH[dataFCH$ageannee<=15,]
  dataMTF <- data[data$pop=="3F" & data$sexe=="M",]
  dataMTF <- dataMTF[dataMTF$ageannee<=12,]
  dataFTF <- data[data$pop=="3F" & data$sexe=="F",]
  dataFTF <- dataFTF[dataFTF$ageannee<=13,]
  
  data <- rbind(dataMCH, dataFCH, dataMTF, dataFTF)
  
  
  # Calculate fawn mass corrected for the date of capture (brought back to the 27th of january)
  # Based on Douhard et al., 2017
  dataCH <- data[data$pop=="C",]
  dataCH$mass_jul <- ifelse(dataCH$ageannee==1, dataCH$masse+0.012*(57-dataCH$datejulienne), dataCH$masse)
  dataTF <- data[data$pop=="3F",]
  dataTF$mass_jul <- ifelse(dataTF$ageannee==1, dataTF$masse+0.024*(57-dataTF$datejulienne), dataTF$masse)
  
  data <- rbind(dataCH, dataTF)
  
  # Number of body mass obs per age class/pop/sex
  tapply(data$ageannee, data$pop:data$sexe, table)
  
  # Full sample size
  tapply(data$sexe, data$pop, table)
  
  
  # Sample size FGM
  datagc <- data[!is.na(data$FGMngg),]
  tapply(datagc$sexe, datagc$pop, table)
  tapply(datagc$ageannee, datagc$pop:datagc$sexe, table)
  
  
  # A./ Calculate averaged body mass according to age/pop/sex ----
  data$age_factor <- as.factor(data$ageannee)
  
  dataC <- data[data$pop=="C",]
  
  MCH <- subset(dataC, sexe=="M") 
  unique(factor(MCH$numind)) 
  MCH$predict <- ave(MCH$mass_jul, MCH$age_factor, FUN=mean)
  MCH$dev <- MCH$mass_jul - MCH$predict
  # 405 obs / 217 ind
  
  FCH <- subset(dataC, sexe=="F") 
  unique(factor(FCH$numind))
  FCH$predict <- ave(FCH$mass_jul, FCH$age_factor, FUN=mean)
  FCH$dev <- FCH$mass_jul - FCH$predict
  # 512 obs / 224 ind
  
  dataTF <- data[data$pop=="3F",]
  
  MTF <- subset(dataTF, sexe=="M") 
  unique(factor(MTF$numind)) 
  MTF$predict <- ave(MTF$mass_jul, MTF$age_factor, FUN=mean)
  MTF$dev <- MTF$mass_jul - MTF$predict
  # 456 obs / 253 ind
  
  FTF <- subset(dataTF, sexe=="F") 
  unique(factor(FTF$numind))
  FTF$predict <- ave(FTF$mass_jul, FTF$age_factor, FUN=mean)
  FTF$dev <- FTF$mass_jul - FTF$predict
  # 508 obs / 254 ind
  
  data <- rbind(MCH, FCH, MTF, FTF) 
  unique(factor(data$numind)) #1881 obs / 947 ind
  
  
  # Calculate a relative environmental quality per population
  data$annee <- as.factor(data$annee)
  data$mean_cohort_quality <- ave(data$qualite_cohorte, data$pop, FUN=function(x) mean(x, na.rm=T)) # mean cohort quality in CH and TF
  data$rel_cohort_quality <- data$qualite_cohorte - data$mean_cohort_quality # relative cohort quality
}


# B./ Early effect of stress on late body mass ----
# B1./ FGMs ----
datagc <- data


# a) MCHgc: Males, Chiz? ----
{MCHgc1 <- subset(datagc, pop=="C" & sexe=="M" & ageannee >= 2) # 252 obs / 135 ind
unique(factor(MCHgc1$numind))

MCHgc2 <- subset(datagc, pop=="C" & sexe=="M" & ageannee == 1) # 153 obs.
MCHgc2 <- MCHgc2[!is.na(MCHgc2$FGMngg),]
MCHgc2 <- MCHgc2[!is.na(MCHgc2$numind),]
MCHgc2 <- MCHgc2[!MCHgc2$idkit=="2014-79-031",] # removing an obs for an individual with two captures the same year
unique(factor(MCHgc2$numind)) # 96 obs / 96 ind

MCHgc2 <- MCHgc2[,c(2,19)]
colnames(MCHgc2)[2] <- "fgm1" # FGM level during the 1st year

MCHgc <- merge(MCHgc1, MCHgc2,  by="numind", all=F)
unique(factor(MCHgc$numind)) # 87 obs / 44 ind


# b) FCHgc: Females, Chiz? ----
FCHgc1 <- subset(datagc, pop=="C" & sexe=="F" & ageannee >= 2) # 384 obs / 149 ind
unique(factor(FCHgc1$numind))

FCHgc2 <- subset(datagc, pop=="C" & sexe=="F" & ageannee == 1)
unique(factor(FCHgc2$numind)) 
FCHgc2 <- FCHgc2[!is.na(FCHgc2$FGMngg),]
FCHgc2 <- FCHgc2[!is.na(FCHgc2$numind),]
unique(factor(FCHgc2$numind)) # 78 obs / 78 ind

FCHgc2 <- FCHgc2[,c(2,19)]
colnames(FCHgc2)[2] <- "fgm1"

FCHgc <- merge(FCHgc1, FCHgc2,  by="numind", all=F)
unique(factor(FCHgc$numind)) # 78 obs / 34 ind

# c) MTFgc: Males, Trois-Fontaines ----
MTFgc1 <- subset(datagc, pop=="3F" & sexe=="M" & ageannee >= 2) # 278 obs / 128 ind
unique(factor(MTFgc1$numind))

MTFgc2 <- subset(datagc, pop=="3F" & sexe=="M" & ageannee == 1) 
unique(factor(MTFgc2$numind)) 
MTFgc2 <- MTFgc2[!is.na(MTFgc2$FGMngg),]
MTFgc2 <- MTFgc2[!is.na(MTFgc2$numind),]
unique(factor(MTFgc2$numind)) # 103 obs / 103 ind

MTFgc2 <- MTFgc2[,c(2,19)]
colnames(MTFgc2)[2] <- "fgm1"

MTFgc <- merge(MTFgc1, MTFgc2,  by="numind", all=F)
unique(factor(MTFgc$numind)) # 81 obs / 39 ind


# d) FTFgc: Females, Trois-Fontaines ----
FTFgc1 <- subset(datagc, pop=="3F" & sexe=="F" & ageannee >= 2) 
FTFgc1 <- FTFgc1[!FTFgc1$idkit=="2010-51-043",] # removing an obs for an individual with two captures the same year
unique(factor(FTFgc1$numind)) # 344 obs / 150 ind

FTFgc2 <- subset(datagc, pop=="3F" & sexe=="F" & ageannee == 1) 
unique(factor(FTFgc2$numind)) 
FTFgc2 <- FTFgc2[!is.na(FTFgc2$FGMngg),]
FTFgc2 <- FTFgc2[!is.na(FTFgc2$numind),]
unique(factor(FTFgc2$numind)) # 92 obs / 92 ind

FTFgc2 <- FTFgc2[,c(2,19)]
colnames(FTFgc2)[2] <- "fgm1"

FTFgc <- merge(FTFgc1, FTFgc2,  by="numind", all=F)
unique(factor(FTFgc$numind)) # 102 obs / 43 ind
}

# FIGURES ----
# Fig3A - Males

MCHgc$logfgm1 <- log(MCHgc$fgm1)
moda1 <- lmer(dev ~ scale(logfgm1, scale=F) + (1|numind), data=MCHgc)
summary(moda1)
r.squaredGLMM(moda1)

newdataa1=expand.grid(logfgm1=seq(min(log(MCHgc$fgm1)),
                                max(log(MCHgc$fgm1)), 0.005))

preda1=predictSE(moda1,newdata=newdataa1,se.fit=TRUE)
newdataa1=cbind(newdataa1,preda1)
newdataa1$low=newdataa1$fit-1.44*newdataa1$se.fit
newdataa1$upp=newdataa1$fit+1.44*newdataa1$se.fit

MTFgc$logfgm1 <- log(MTFgc$fgm1)
moda2 <- lmer(dev ~ scale(logfgm1, scale=F) + (1|numind), data=MTFgc)
summary(moda2)
r.squaredGLMM(moda2)

newdataa2=expand.grid(logfgm1=seq(min(log(MTFgc$fgm1)),
                                max(log(MTFgc$fgm1)), 0.005))

preda2=predictSE(moda2,newdata=newdataa2,se.fit=TRUE)
newdataa2=cbind(newdataa2,preda2)
newdataa2$low=newdataa2$fit-1.44*newdataa2$se.fit
newdataa2$upp=newdataa2$fit+1.44*newdataa2$se.fit

groba <- grobTree(textGrob("Males", x=0.05,  y=0.97, hjust=0,
                           gp=gpar(col="black", fontsize=14, face="bold")))
fig3a <- ggplot(data = MCHgc, aes(y=dev, x=logfgm1)) + 
  geom_point(data = MCHgc, aes(y=dev, x=logfgm1, color="#008EA0FF")) +
  geom_point(data = MTFgc, aes(y=dev, x=logfgm1, color="#FF6F00FF")) +
  scale_y_continuous(name = "Relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("Juvenile FGM levels (log-transformed)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  annotation_custom(groba) +
  labs(tag="A") +
  theme(plot.margin=margin(0,1,0,0, "cm"),
        plot.title = element_text(face = "bold", hjust=0.5, size=18), # plot title
        axis.title.x = element_text(size=12), #change size x axis title
        axis.text.x = element_text(size=11.5), #change size text x axis
        axis.title.y = element_text(size=12), #change size y axis title
        axis.text.y = element_text(size=11.5), #change size text y axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # no grid in the background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=10),#change legend text font size
        legend.position=c(0.2, 0.1),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#008EA0FF", "#FF6F00FF"), 
                        name="Population",
                        labels=c("Chizé", "Trois-Fontaines"))
fig3a


# Fig3B - Females

FCHgc$logfgm1 <- log(FCHgc$fgm1)
modb1 <- lmer(dev ~ scale(logfgm1, scale=F) + scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc)
summary(modb1)
r.squaredGLMM(modb1)

newdatab1=expand.grid(logfgm1=seq(min(log(FCHgc$fgm1)),
                                  max(log(FCHgc$fgm1)), 0.005),
                      qualite_cohorte=mean(FCHgc$qualite_cohorte))

predb1=predictSE(modb1,newdata=newdatab1,se.fit=TRUE)
newdatab1=cbind(newdatab1,predb1)
newdatab1$low=newdatab1$fit-1.44*newdatab1$se.fit
newdatab1$upp=newdatab1$fit+1.44*newdatab1$se.fit

FTFgc <- FTFgc[FTFgc$fgm1<5000,] 
FTFgc$logfgm1 <- log(FTFgc$fgm1)
modb2 <- lmer(dev ~ scale(logfgm1, scale=F) + (1|numind), data=FTFgc)
summary(modb2)
r.squaredGLMM(modb2)

newdatab2=expand.grid(logfgm1=seq(min(log(FTFgc$fgm1)),
                                  max(log(FTFgc$fgm1)), 0.005))

predb2=predictSE(modb2,newdata=newdatab2,se.fit=TRUE)
newdatab2=cbind(newdatab2,predb2)
newdatab2$low=newdatab2$fit-1.44*newdatab2$se.fit
newdatab2$upp=newdatab2$fit+1.44*newdatab2$se.fit

grobb <- grobTree(textGrob("Females", x=0.05,  y=0.97, hjust=0,
                           gp=gpar(col="black", fontsize=14, face="bold")))
fig3b <- ggplot(data = FCHgc, aes(y=dev, x=logfgm1)) + 
  geom_point(data = FCHgc, aes(y=dev, x=logfgm1, color="#008EA0FF")) +
  geom_point(data = FTFgc, aes(y=dev, x=logfgm1, color="#FF6FOOFF")) +
  scale_y_continuous(name = "Relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("Juvenile FGM levels (log-transformed)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  annotation_custom(grobb) +
  labs(tag="B") +
  theme(plot.margin=margin(0,1,0,0, "cm"),
        plot.title = element_text(face = "bold", hjust=0.5, size=18), # plot title
        axis.title.x = element_text(size=12), #change size x axis title
        axis.text.x = element_text(size=11.5), #change size text x axis
        axis.title.y = element_text(size=12), #change size y axis title
        axis.text.y = element_text(size=11.5), #change size text y axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # no grid in the background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_text(size=11), #change legend title font size
        legend.text = element_text(size=10),#change legend text font size
        legend.position=c(0.2, 0.1),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#008EA0FF", "#FF6F00FF"), 
                     name="Population",
                     labels=c("Chizé", "Trois-Fontaines"))
fig3b

f3 <- arrangeGrob(fig3a, fig3b, ncol=2, nrow=1)
grid.newpage()
grid.draw(f3)

# Save this figure
ggsave(filename = "results/figures/fig3.png", 
       f3, width = 30, height = 15, dpi = 300, units = "cm", device='png')

