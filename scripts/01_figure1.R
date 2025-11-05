############################# ANALYSE STRESS-MASS ##############################
#################### SHORT-TERM RELATION STRESS X BODY MASS ####################
################# A./ PREDICTIONS POIDS MOYENS-----------------#################
################# B./ FIGURES----------------------------------#################
################################# 16/06/2022 ###################################
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
  
  # B./ Relation between body mass and Stress ----
  # B1./ FGM ----
  datagc <- data[!is.na(data$FGMngg),] # 982
  
  # a) MCHgc: Males, Chiz? ----
  MCHgc <- subset(datagc, pop=="C" & sexe=="M")
  MCHgc <- MCHgc[!is.na(MCHgc$numind),]
  unique(factor(MCHgc$numind)) 
  # 233 obs / 149 ind
  
  par(mfrow=c(2,2))
  plot(MCHgc$age, MCHgc$predict)
  plot(MCHgc$age, MCHgc$dev)
  plot(log(MCHgc$FGMngg) ~ MCHgc$age)
  hist(MCHgc$dev)
  
  # b) FCHgc: Females, Chiz? ----
  FCHgc <- subset(datagc, pop=="C" & sexe=="F")
  FCHgc <- FCHgc[!is.na(FCHgc$numind),]
  unique(factor(FCHgc$numind)) 
  # 297 obs / 160 ind (with high FGM value)
  # 296 obs / 159 ind
  
  par(mfrow=c(2,2))
  plot(FCHgc$age, FCHgc$predict)
  plot(FCHgc$age, FCHgc$dev)
  plot(log(FCHgc$FGMngg) ~ FCHgc$age)
  hist(FCHgc$dev)
  
  
  # c) MTFgc: Males, Trois-Fontaines ----
  MTFgc <- subset(datagc, pop=="3F" & sexe=="M")
  MTFgc <- MTFgc[!is.na(MTFgc$numind),]
  unique(factor(MTFgc$numind)) 
  # 240 obs / 159 ind
  
  par(mfrow=c(2,2))
  plot(MTFgc$age, MTFgc$predict)
  plot(MTFgc$age, MTFgc$dev)
  plot(log(MTFgc$FGMngg) ~ MTFgc$age)
  hist(MTFgc$dev)
  
  
  # d) FTFgc: Females, Trois-Fontaines ----
  FTFgc <- subset(datagc, pop=="3F" & sexe=="F")
  FTFgc <- FTFgc[!is.na(FTFgc$numind),]
  unique(factor(FTFgc$numind)) 
  # 256 obs / 161 ind (with high FGM value)
  # 255 obs / 161 ind
  
  par(mfrow=c(2,2))
  plot(FTFgc$age, FTFgc$predict)
  plot(FTFgc$age, FTFgc$dev)
  plot(log(FTFgc$FGMngg) ~ FTFgc$age)
  hist(FTFgc$dev)
  
  
  datagc <- rbind(FCHgc, MCHgc, FTFgc, MTFgc)
  unique(factor(datagc$numind)) 
  # 1026 obs / 629 ind (with high FGM values)
  # 1024 obs / 628 ind
  
  par(mfrow=c(1,1))
  
  # FGM sample size according to age/pop/sex
  tapply(datagc$ageannee, datagc$pop:datagc$sexe, table)
  table(datagc$pop:datagc$sexe)
}


# B./ FIGURES ----
# a) Figure 1A ----
# Males, CHIZE ----
MCHgc1 <- MCHgc[MCHgc$ageannee==1,]
MCHgc1 <- MCHgc1[!MCHgc1$idkit=="2014-79-031",] # removing the only individuals with two captures at 1 yo
unique(factor(MCHgc1$numind)) # 96 ind/obs

# Check how the slope of FGM changes with environmental quality
MCHgc1$logfgm <- scale(log(MCHgc1$FGMngg), scale=F)
MCHgc1$qual_coh <- scale(MCHgc1$qualite_cohorte, scale=F)
mod1 <- lm(mass_jul ~ qual_coh*logfgm, data=MCHgc1)

summary(mod1)


fig1asupp <- interplot(mod1, var1 = "logfgm", var2 = "qual_coh") +
  xlab("Environmental quality") +
  ylab("FGM (log-transformed)") +
  ggtitle("Estimated coefficient (slope) for FGM on body mass according\nto environmental quality")
fig1asupp
ggsave(filename = "results/figures/fig1a-slope.png", 
       fig1asupp, width = 15, height = 15, dpi = 300, units = "cm", device='png')

MCHgc1$LogFGM <- log(MCHgc1$FGMngg)
hist(MCHgc1$qualite_cohorte)
MCHgc1a <- subset(MCHgc1, qualite_an<13)
MCHgc1b <- subset(MCHgc1, qualite_an>13)

mod1ca <- lm(mass_jul ~ scale(qualite_an, scale=F) + scale(LogFGM, scale=F), data=MCHgc1a)

min(MCHgc1a$LogFGM)
max(MCHgc1a$LogFGM)

newdataMCHa=expand.grid(LogFGM=seq(min(MCHgc1a$LogFGM),
                                max(MCHgc1a$LogFGM), 0.05),
                     qualite_an=mean(MCHgc1a$qualite_an))

predMCHa=predict(mod1ca,newdata=newdataMCHa,se.fit=TRUE)
newdataMCHa=cbind(newdataMCHa,predMCHa)
newdataMCHa$low=newdataMCHa$fit-1.44*newdataMCHa$se.fit
newdataMCHa$upp=newdataMCHa$fit+1.44*newdataMCHa$se.fit

mod1cMCHb <- lm(mass_jul ~ scale(qualite_an, scale=F) + scale(LogFGM, scale=F), data=MCHgc1b)

min(MCHgc1b$LogFGM)
max(MCHgc1b$LogFGM)

newdataMCHb=expand.grid(LogFGM=seq(min(MCHgc1b$LogFGM),
                                max(MCHgc1b$LogFGM), 0.05),
                     qualite_an=mean(MCHgc1b$qualite_an))

predMCHb=predict(mod1cMCHb,newdata=newdataMCHb,se.fit=TRUE)
newdataMCHb=cbind(newdataMCHb,predMCHb)
newdataMCHb$low=newdataMCHb$fit-1.44*newdataMCHb$se.fit
newdataMCHb$upp=newdataMCHb$fit+1.44*newdataMCHb$se.fit





# b) Figure 1B ----
# Males, TF ----
MTFgc1 <- MTFgc[MTFgc$ageannee==1,]
unique(factor(MTFgc1$numind)) # 103 ind/obs

MTFgc1$LogFGM <- log(MTFgc1$FGMngg)

mod1b <- lm(mass_jul ~ scale(qualite_an, scale=F) + scale(LogFGM, scale=F), data=MTFgc1)# modèle 4
min(log(MTFgc1$FGMngg))
max(log(MTFgc1$FGMngg))

newdataMTF=expand.grid(LogFGM=seq(min(log(MTFgc1$FGMngg)),
                               max(log(MTFgc1$FGMngg)),0.05),
                    qualite_an=mean(MTFgc1$qualite_an))

predMTF=predict(mod1b,newdata=newdataMTF,se.fit=TRUE)
newdataMTF=cbind(newdataMTF,predMTF)
newdataMTF$low=newdataMTF$fit-1.44*newdataMTF$se.fit
newdataMTF$upp=newdataMTF$fit+1.44*newdataMTF$se.fit

show_col(pal_futurama("planetexpress", alpha = 1)(12))

grob1a <- grobTree(textGrob("Juvenile males", x=0.05,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=14, face="bold")))
fig1A <- ggplot(data = MCHgc1, aes(y=mass_jul, x=LogFGM)) + 
  geom_point(data = MCHgc1a, aes(y=mass_jul, x=LogFGM, color="#008EA0FF"), shape=16) +
  geom_point(data = MCHgc1b, aes(y=mass_jul, x=LogFGM, color="#84D7E1FF"), shape=16) +
  geom_ribbon(data = newdataMCHa, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMCHa, aes(x=LogFGM, y=fit, color="#008EA0FF"), size=1) +
  geom_ribbon(data = newdataMCHb, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMCHb, aes(x=LogFGM, y=fit, color="#84D7E1FF"), size=1) +
  geom_point(data = MTFgc1, aes(y=mass_jul, x=LogFGM, color="#FF6F00FF"), shape=16) +
  geom_ribbon(data = newdataMTF, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMTF, aes(x=LogFGM, y=fit, color="#FF6F00FF"), size=1) +
  scale_y_continuous(name = "Mass (kg, corrected for the date of capture)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
  labs(tag="A") +
  annotation_custom(grob1a) +
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
        legend.position=c(0.32, 0.1),
        legend.background = element_rect(fill='transparent')) +
  scale_colour_manual(values=c("#008EA0FF", "#84D7E1FF", "#FF6F00FF"), 
                      name="Population",
                      labels=c("Chizé - Low cohort quality (<12.5kg)", "Chizé - High cohort quality (>13kg)", "Trois-Fontaines")) 

fig1A


# c) Figure 1C ----
# Females, CH ----
FCHgc1 <- FCHgc[FCHgc$ageannee==1,]
unique(factor(FCHgc1$numind)) # 78 ind/obs

FCHgc1$LogFGM <- log(FCHgc1$FGMngg)

mod1c <- lm(mass_jul ~ scale(LogFGM, scale=F) + scale(qualite_cohorte, scale=F), data=FCHgc1)# modèle 4
min(log(FCHgc1$FGMngg))
max(log(FCHgc1$FGMngg))

newdataFCH=expand.grid(LogFGM=seq(min(log(FCHgc1$FGMngg)),
                               max(log(FCHgc1$FGMngg)),0.05),
                       qualite_cohorte=mean(FCHgc1$qualite_cohorte))

predFCH=predict(mod1c,newdata=newdataFCH,se.fit=TRUE)
newdataFCH=cbind(newdataFCH,predFCH)
newdataFCH$low=newdataFCH$fit-1.44*newdataFCH$se.fit
newdataFCH$upp=newdataFCH$fit+1.44*newdataFCH$se.fit


# d) Figure 1D ----
# Females, TF ----
FTFgc1 <- FTFgc[FTFgc$ageannee==1,]
FTFgc1 <- FTFgc1[FTFgc1$FGMngg<=6000,] 
unique(factor(FTFgc1$numind)) # 91 ind/obs

FTFgc1$LogFGM <- log(FTFgc1$FGMngg)

mod1d <- lm(mass_jul ~ scale(LogFGM, scale=F) + scale(qualite_cohorte, scale=F), data=FTFgc1)# modèle 4
min(log(FTFgc1$FGMngg))
max(log(FTFgc1$FGMngg))

newdataFTF=expand.grid(LogFGM=seq(min(log(FTFgc1$FGMngg)),
                               max(log(FTFgc1$FGMngg)),0.05),
                       qualite_cohorte=mean(FTFgc1$qualite_cohorte))

predFTF=predict(mod1d,newdata=newdataFTF,se.fit=TRUE)
newdataFTF=cbind(newdataFTF,predFTF)
newdataFTF$low=newdataFTF$fit-1.44*newdataFTF$se.fit
newdataFTF$upp=newdataFTF$fit+1.44*newdataFTF$se.fit

grob1d <- grobTree(textGrob("Juvenile females", x=0.05,  y=0.97, hjust=0,
                            gp=gpar(col="black", fontsize=14, face="bold")))
fig1B <- ggplot(data = FCHgc1, aes(y=mass_jul, x=log(FGMngg))) + 
  geom_point(data = FCHgc1, aes(y=mass_jul, x=log(FGMngg), color="#008EA0FF"), shape=16) +
  geom_point(data = FTFgc1, aes(y=mass_jul, x=log(FGMngg), color="#FF6F00FF"), shape=16) +
  scale_y_continuous(name = "Mass (kg, corrected for the date of capture)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
  labs(tag="B") +
  annotation_custom(grob1d) +
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
        legend.position=c(0.2, 0.83),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#008EA0FF", "#FF6F00FF"), 
                        name="Population",
                        labels=c("Chizé", "Trois-Fontaines"))
fig1B



# e) Figure 1E ----
# Males, Chizé ----
MCHgc2 <- MCHgc[!MCHgc$age==1,]
unique(factor(MCHgc2$numind))
# 136 obs / 83 ind

MCHgc2$LogFGM <- log(MCHgc2$FGMngg)

mod2a <- lmer(dev ~ scale(LogFGM, scale=F) + (1|numind), data=MCHgc2)
min(log(MCHgc2$FGMngg))
max(log(MCHgc2$FGMngg))

newdataMCH2=expand.grid(LogFGM=seq(min(log(MCHgc2$FGMngg)),
                                   max(log(MCHgc2$FGMngg)),0.05))

predMCH2=predictSE(mod2a,newdata=newdataMCH2,se.fit=TRUE)
newdataMCH2=cbind(newdataMCH2,predMCH2)
newdataMCH2$low=newdataMCH2$fit-1.44*newdataMCH2$se.fit
newdataMCH2$upp=newdataMCH2$fit+1.44*newdataMCH2$se.fit

# Males, TF ----
MTFgc2 <- MTFgc[!MTFgc$age==1,]
unique(factor(MTFgc2$numind))
# 137 obs / 84 ind

MTFgc2$LogFGM <- log(MTFgc2$FGMngg)

# Check how the slope of FGM changes with environmental quality
MTFgc2$logfgm <- scale(log(MTFgc2$FGMngg), scale=F)
MTFgc2$qual_coh <- scale(MTFgc2$qualite_cohorte, scale=F)
mod1 <- lm(dev ~ qual_coh*logfgm, data=MTFgc2)

summary(mod1)


fig1csupp <- interplot(mod1, var1 = "logfgm", var2 = "qual_coh") +
  xlab("Environmental quality") +
  ylab("FGM (log-transformed)") +
  ggtitle("Estimated coefficient (slope) for FGM on body mass according\nto environmental quality")
fig1csupp
ggsave(filename = "results/figures/fig1c-slope.png", 
       fig1csupp, width = 15, height = 15, dpi = 300, units = "cm", device='png')

hist(MTFgc2$qualite_cohorte)
MTFgc2a <- subset(MTFgc2, qualite_cohorte<16)
MTFgc2b <- subset(MTFgc2, qualite_cohorte>16)

mod2ba <- lm(dev ~ scale(qualite_cohorte, scale=F) + scale(LogFGM, scale=F), data=MTFgc2a)

min(MTFgc2a$LogFGM)
max(MTFgc2a$LogFGM)

newdataMTF2a=expand.grid(LogFGM=seq(min(MTFgc2a$LogFGM),
                                   max(MTFgc2a$LogFGM), 0.05),
                        qualite_cohorte=mean(MTFgc2a$qualite_cohorte))

predMTF2a=predict(mod2ba,newdata=newdataMTF2a,se.fit=TRUE)
newdataMTF2a=cbind(newdataMTF2a,predMTF2a)
newdataMTF2a$low=newdataMTF2a$fit-1.44*newdataMTF2a$se.fit
newdataMTF2a$upp=newdataMTF2a$fit+1.44*newdataMTF2a$se.fit

mod2bb <- lm(dev ~ scale(qualite_cohorte, scale=F) + scale(LogFGM, scale=F), data=MTFgc2b)

min(MTFgc2b$LogFGM)
max(MTFgc2b$LogFGM)

newdataMTF2b=expand.grid(LogFGM=seq(min(MTFgc2b$LogFGM),
                                   max(MTFgc2b$LogFGM), 0.05),
                        qualite_cohorte=mean(MTFgc2b$qualite_cohorte))

predMTF2b=predict(mod2bb,newdata=newdataMTF2b,se.fit=TRUE)
newdataMTF2b=cbind(newdataMTF2b,predMTF2b)
newdataMTF2b$low=newdataMTF2b$fit-1.44*newdataMTF2b$se.fit
newdataMTF2b$upp=newdataMTF2b$fit+1.44*newdataMTF2b$se.fit

show_col(pal_futurama("planetexpress", alpha = 0.3)(12))

grob1c <- grobTree(textGrob("Adult males", x=0.05,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=14, face="bold")))
fig1C <- ggplot(data = MTFgc2, aes(y=dev, x=LogFGM)) + 
  geom_point(data = MTFgc2a, aes(y=dev, x=LogFGM, color="#FF6F004C"), shape=16) +
  geom_point(data = MTFgc2b, aes(y=dev, x=LogFGM, color="#FF6F00FF"), shape=16) +
  geom_point(data = MCHgc2, aes(y=dev, x=LogFGM, color="#008EA0FF"), shape=16) +
  geom_ribbon(data = newdataMCH2, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMCH2, aes(x=LogFGM, y=fit, color="#008EA0FF"), size=1) +
  geom_ribbon(data = newdataMTF2a, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMTF2a, aes(x=LogFGM, y=fit, color="#FF6F004C"), size=1) +
  geom_ribbon(data = newdataMTF2b, aes(y=fit, ymin=low, ymax=upp, x=LogFGM), alpha = 0.1) +
  geom_line(data = newdataMTF2b, aes(x=LogFGM, y=fit, color="#FF6F00FF"), size=1) +
  scale_y_continuous(name = "Relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
  labs(tag="C") +
  annotation_custom(grob1c) +
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
        legend.position=c(0.37, 0.15),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#008EA0FF", "#FF6F00FF", "#FF6F004C"), 
                      name="Population",
                      labels=c("Chizé", "Trois-Fontaines - Low cohort quality (<16kg)", "Trois-Fontaines - High cohort quality (>16kg)")) 
fig1C

# Females, CH ----
FCHgc2 <- FCHgc[!FCHgc$age==1,]
FCHgc2 <- FCHgc2[FCHgc2$FGMngg<=5000,] 
unique(factor(FCHgc2$numind))
# 218 obs / 104 ind

FCHgc2$LogFGM <- log(FCHgc2$FGMngg)

mod1c <- lmer(dev ~ scale(LogFGM, scale=F) + scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc2)
min(log(FCHgc2$FGMngg))
max(log(FCHgc2$FGMngg))

newdataFCH2=expand.grid(LogFGM=seq(min(log(FCHgc2$FGMngg)),
                                   max(log(FCHgc2$FGMngg)),0.05),
                        qualite_cohorte=mean(FCHgc2$qualite_cohorte))

predFCH2=predictSE(mod1c,newdata=newdataFCH2,se.fit=TRUE)
newdataFCH2=cbind(newdataFCH2,predFCH2)
newdataFCH2$low=newdataFCH2$fit-1.44*newdataFCH2$se.fit
newdataFCH2$upp=newdataFCH2$fit+1.44*newdataFCH2$se.fit

# Females, TF ----
FTFgc2 <- FTFgc[!FTFgc$age==1,]
unique(factor(FTFgc2$numind))
# 152 obs / 98 ind

FTFgc2$LogFGM <- log(FTFgc2$FGMngg)

mod1c <- lmer(dev ~ scale(LogFGM, scale=F) + (1|numind), data=FTFgc2)
min(log(FTFgc2$FGMngg))
max(log(FTFgc2$FGMngg))

newdataFTF2=expand.grid(LogFGM=seq(min(log(FTFgc2$FGMngg)),
                               max(log(FTFgc2$FGMngg)),0.05))

predFTF2=predictSE(mod1c,newdata=newdataFTF2,se.fit=TRUE)
newdataFTF2=cbind(newdataFTF2,predFTF2)
newdataFTF2$low=newdataFTF2$fit-1.44*newdataFTF2$se.fit
newdataFTF2$upp=newdataFTF2$fit+1.44*newdataFTF2$se.fit


grob1d <- grobTree(textGrob("Adult females", x=0.05,  y=0.97, hjust=0,
                            gp=gpar(col="black", fontsize=14, face="bold")))
fig1D <- ggplot(data = FCHgc2, aes(y=dev, x=log(FGMngg))) + 
  geom_point(data = FCHgc2, aes(y=dev, x=log(FGMngg), color="#008EA0FF"), shape=16) +
  geom_point(data = FTFgc2, aes(y=dev, x=log(FGMngg), color="#FF6F00FF"), shape=16) +
  scale_y_continuous(name = "Relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
  labs(tag="D") +
  annotation_custom(grob1d) +
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
        legend.position=c(0.2, 0.83),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#008EA0FF", "#FF6F00FF"), 
                     name="Population",
                     labels=c("Chizé", "Trois-Fontaines"))
fig1D



f1 <- arrangeGrob(fig1A, fig1B, fig1C, fig1D, ncol=2, nrow=2)
grid.newpage()
grid.draw(f1)

# Save this figure
ggsave(filename = "results/figures/fig1.png", 
       f1, width = 30, height = 30, dpi = 300, units = "cm", device='png')


