############################# ANALYSE STRESS-MASS ##############################
#################### EARLY-LATE RELATION STRESS X BODY MASS ####################
################# A./ PREDICTIONS POIDS MOYENS-----------------#################
################# B./ ANALYSES---------------------------------#################
#################   1) FGM-------------------------------------#################
################################# 31/05/2022 ###################################
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

# B2./ Models ----
# 0) Pooling
malesgc2 <- rbind(MCHgc, MTFgc)
unique(factor(malesgc2$numind)) # 134 obs / 71 males
femalesgc2 <- rbind(FCHgc, FTFgc)
unique(factor(femalesgc2$numind)) # 141 obs / 65 females


# pop-sex pooled ----
datagc2a <- rbind(MCHgc, FCHgc, MTFgc, FTFgc) # 348 obs / 160 ind
datagc2 <- datagc2a[datagc2a$fgm1<5000,] 
unique(factor(datagc2$numind)) # 345 obs / 159 ind

mod <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(rel_cohort_quality, scale=F) + scale(log(fgm1), scale=F)*pop + scale(log(fgm1), scale=F)*sexe + (1|numind), data=datagc2)
summary(mod)

qqnorm(resid(mod))
qqline(resid(mod))

options(na.action="na.fail")
s1pool <- dredge(mod, rank="AICc", evaluate=T)
s1pool
write.table(s1pool, "results/model selection tables/longpool.csv", row.names = F, sep=";", dec=".")

mod1 <- lmer(dev ~ scale(rel_cohort_quality, scale=F) + (1|numind), data=malesgc2)
summary(mod1)
-0.3933-1.44*0.2136;-0.3933+1.44*0.2136
0.4187-1.44*0.2545;0.4187+1.44*0.2545
r.squaredGLMM(mod1)

# with high fgm value
mod <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(rel_cohort_quality, scale=F) + scale(log(fgm1), scale=F)*pop + scale(log(fgm1), scale=F)*sexe + (1|numind), data=datagc2a)
summary(mod)

qqnorm(resid(mod))
qqline(resid(mod))

options(na.action="na.fail")
s1poolhigh <- dredge(mod, rank="AICc", evaluate=T)
s1poolhigh
write.table(s1poolhigh, "results/model selection tables/longpoolhigh.csv", row.names = F, sep=";", dec=".")

mod1 <- lmer(dev ~ scale(rel_cohort_quality, scale=F) + (1|numind), data=datagc2a)
summary(mod1)
-0.2916-1.44*0.1592;-0.2916+1.44*0.1592
0.5332-1.44*0.1983;0.5332+1.44*0.1983
r.squaredGLMM(mod1)


# a) MCH ----
# Quality cohort ----
m1b <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=MCHgc)
summary(m1b)

options(na.action="na.fail")
s1mch <- dredge(m1b, rank="AICc", evaluate=T)
s1mch
write.table(s1mch, "results/model selection tables/longMCH.csv", row.names = F, sep=";", dec=".")

dd=subset(s1, delta<2)
dd

nst=nested(dd)
nst
#dd$nested <- sapply(nst, paste, collapse = ",") # dit si nested ou non
d=subset(dd, !nested(.)) # table avec que les mod?les non nich?s
d

mod1 <- lmer(dev ~ (1|numind), data=MCHgc)
summary(mod1)
-0.1443-1.44*0.2492;-0.1443+1.44*0.2492

r.squaredGLMM(mod1)

# b) FCH ----
# Quality cohort ----
m2b <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc)
summary(m2b)

options(na.action="na.fail")
s1fch <- dredge(m2b, rank="AICc", evaluate=T)
s1fch
write.table(s1fch, "results/model selection tables/longFCH.csv", row.names = F, sep=";", dec=".")

dd=subset(s1, delta<2)
dd


# Mod?le averaging #
mod1 <- lmer(dev ~ scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc)
summary(mod1)
-0.4138-1.44*0.2835;-0.4138+1.44*0.2835
0.8510-1.44*0.3453;0.8510+1.44*0.3453

# R2m et R2c
r.squaredGLMM(mod1)

# Normalité du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# Normalit? des effets al?atoires
plot(ranef(mod1))[1]


# c) MTF ----
# Quality cohort ----
m3b <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=MTFgc)
summary(m3b)

options(na.action="na.fail")
s1mtf <- dredge(m3b, rank="AICc", evaluate=T)
s1mtf
write.table(s1mtf, "results/model selection tables/longMTF.csv", row.names = F, sep=";", dec=".")


# Mod?le averaging #
mod1b <- lmer(dev ~ (1|numind), data=MTFgc)
summary(mod1b)
-0.7395-1.44*0.365;-0.7395+1.44*0.365

# R2m et R2c
r.squaredGLMM(mod1b)

# Normalité du meilleur modele
hist(residuals(mod1b))
qqnorm(residuals(mod1b))
qqline(residuals(mod1b))
shapiro.test(residuals(mod1b))
plot(mod1b,ask=TRUE)

# Normalit? des effets al?atoires
plot(ranef(mod1b))[1]


# d) FTF ----
# Quality cohort ----
m4b <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FTFgc)
summary(m4b)

options(na.action="na.fail")
s1ftfhigh <- dredge(m4b, rank="AICc", evaluate=T)
s1ftfhigh
write.table(s1ftfhigh, "results/model selection tables/longFTFhigh.csv", row.names = F, sep=";", dec=".")

mod1 <- lmer(dev ~ (1|numind), data=FTFgc)
summary(mod1)
-0.07136-1.44*0.36125;-0.07136+1.44*0.36125
r.squaredGLMM(mod1)

# without high value
FTFgca <- FTFgc[FTFgc$fgm1<5000,] # 99 obs / 42 ind
unique(factor(FTFgca$numind))
m4ba <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FTFgca)
summary(m4b)

options(na.action="na.fail")
s1ftf <- dredge(m4ba, rank="AICc", evaluate=T)
s1ftf
write.table(s1ftf, "results/model selection tables/longFTF.csv", row.names = F, sep=";", dec=".")

mod1 <- lmer(dev ~ (1|numind), data=FTFgca)
summary(mod1)
-0.1454-1.44*0.3625;-0.1454+1.44*0.3625
r.squaredGLMM(mod1)


################################# REVIEW OIKOS #################################
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
  
  
  # Calculate averaged body mass according to age/pop/sex
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


# Early effect of stress on late body mass
datagc <- data


# a) MCHgc: Males, Chizé
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


# b) FCHgc: Females, Chizé
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

# c) MTFgc: Males, Trois-Fontaines
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


# d) FTFgc: Females, Trois-Fontaines
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

datagc <- rbind(MCHgc, FCHgc, MTFgc, FTFgc)
unique(factor(datagc$numind)) # 348 obs / 160 ind
}

datagcA <- datagc[datagc$fgm1<5000,] 
unique(factor(datagcA$numind)) # 345 obs / 159 ind
datagcA$age_diff <- datagcA$ageannee-1

mod <- lmer(dev ~ scale(log(fgm1), scale=F)*scale(rel_cohort_quality, scale=F)*pop*sexe + (1|numind) + (1|annee), data=datagcA)
summary(mod)

qqnorm(resid(mod))
qqline(resid(mod))

options(na.action="na.fail")
dr <- dredge(mod, rank="AICc", evaluate=T, trace=2)
ddr <- subset(dr, delta<=2)
ddr

d <- subset(ddr, !nested(.))
d

write.table(dr, "results/Review OIKOS/results_long-term.csv", row.names = F, sep=";", dec=".")
