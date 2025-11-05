############################# ANALYSE STRESS-MASS ##############################
#################### SHORT-TERM RELATION STRESS X BODY MASS ####################
################# A./ PREDICTIONS POIDS MOYENS-----------------#################
################# B./ ANALYSES---------------------------------#################
#################   1) FGM-------------------------------------#################
#################     a) 1 year old----------------------------#################
#################     b) 2 year old and more-------------------#################
################################# 30/05/2022 ###################################
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
unique(factor(data$numind)) # 1881 obs / 947 ind


# Calculate a relative environmental quality per population
data$annee <- as.factor(data$annee)
data$mean_cohort_quality <- ave(data$qualite_cohorte, data$pop, FUN=function(x) mean(x, na.rm=T)) # mean cohort quality in CH and TF
data$rel_cohort_quality <- data$qualite_cohorte - data$mean_cohort_quality # relative cohort quality

# B./ Relation between body mass and Stress ----
# B1./ FGM ----
  datagc <- data[!is.na(data$FGMngg),] # 1026
  
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

cohort <- datagc[,c(3,4,5)]
cohort <- cohort[!duplicated(cohort[c(1,2,3)]),]
hist(cohort$cohorte)
hist(cohort$qualite_cohorte)
cor.test(cohort$cohorte, cohort$qualite_cohorte, method="spearman") 
# -0.14 pearson; -0.11 spearman, -0.06 kendall en prenant tout le jeu de données
# -0.26 pearson; -0.24 spearman, -0.16 kendall en prenant qu'une obs par année
plot(cohort$cohorte, cohort$qualite_cohorte)

cohortCH <- cohort[cohort$pop=="C",]
cor.test(cohortCH$cohorte, cohortCH$qualite_cohorte, method="kendall") 
# -0.13 pearson; -0.04 spearman, 0.05 kendall en prenant tout le jeu de données
# -0.15 pearson; -0.05 spearman, 0 kendall en prenant qu'une obs par année
plot(cohortCH$cohorte, cohortCH$qualite_cohorte)

cohortTF <- cohort[cohort$pop=="3F",]
cor.test(cohortTF$cohorte, cohortTF$qualite_cohorte, method="kendall") 
# -0.35 pearson; -0.39 spearman, -0.29 kendall en prenant tout le jeu de données
# -0.31 pearson; -0.27 spearman, -0.19 kendall en prenant qu'une obs par année
plot(cohortTF$cohorte, cohortTF$qualite_cohorte)

data$cohorte <- as.factor(data$cohorte)


  
# B1A./ Only 1 year ----
  # 0) Pooling ----
  datagc1a <- datagc[datagc$ageannee==1,] 
  datagc1a <- datagc1a[!datagc1a$idkit=="2014-79-031",] # removing the only individuals with two captures at 1 yo # 369 obs/ind
  datagc1 <- datagc1a[datagc1a$FGMngg<5000,] 
  unique(factor(datagc1$numind)) # 368 obs/ind
  malesgc1 <- datagc1[datagc1$sexe=="M",]
  unique(factor(malesgc1$numind)) # 194 males
  femalesgc1 <- datagc1[datagc1$sexe=="F",]
  unique(factor(femalesgc1$numind)) # 166 females
  
  mod <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F) + scale(log(FGMngg), scale=F)*pop + scale(log(FGMngg), scale=F)*sexe, data=datagc1)
  summary(mod)
  
  qqnorm(resid(mod))
  qqline(resid(mod))
  
  options(na.action="na.fail")
  s1pool1 <- dredge(mod, rank="AICc", evaluate=T)
  s1pool1
  write.table(s1pool1, "results/model selecton tables/shortpool1.csv", row.names = F, sep=";", dec=".")
  
  dd=subset(s1, delta<2)
  dd # tables of models within 2 delta AICc
  
  nst=nested(dd)
  nst
  d=subset(dd, !nested(.)) # table avec que les mod?les non nich?s
  d
  
  mod1 <- lm(mass_jul ~ pop + scale(log(FGMngg), scale=F) + scale(rel_cohort_quality, scale=F) + sexe, data=datagc1)
  summary(mod1)
  15.2788-1.44*0.1866;15.2788+1.44*0.1866
  -3.1019-1.44*0.2153;-3.1019+1.44*0.2153
  -0.4110-1.44*0.1771;-0.4110+1.44*0.1771
  1.0290-1.44*0.1233;1.0290+1.44*0.1233
  0.5254-1.44*0.2156;0.5254+1.44*0.2156
  r.squaredGLMM(mod1)
  
  # with high fgm value
  modbis <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F) + scale(log(FGMngg), scale=F)*pop + scale(log(FGMngg), scale=F)*sexe, data=datagc1a)
  summary(modbis)
  
  qqnorm(resid(modbis))
  qqline(resid(modbis))
  
  options(na.action="na.fail")
  s1pool1high <- dredge(modbis, rank="AICc", evaluate=T)
  s1pool1high
  write.table(s1pool1high, "results/model selection tables/shortpool1high.csv", row.names = F, sep=";", dec=".")
  
  dd=subset(s1, delta<2)
  dd # tables of models within 2 delta AICc
  
  nst=nested(dd)
  nst
  d=subset(dd, !nested(.)) # table avec que les mod?les non nich?s
  d
  
  ### mod?le averaging avec tous les mod?les gard?s par nesting rule EG 26/01/2021
  avgm <- model.avg(dd,subset=c("80", "16", "48"), fit=T)
  avgm#[["coefficients"]]
  summary(avgm)
  confint(avgm, level=0.85)
  
  mod1 <- lm(mass_jul ~ pop + scale(log(FGMngg), scale=F) + scale(qualite_cohorte, scale=F) + sexe, data=datagc1)
  summary(mod1)
  13.8927-1.44*0.2584;13.8927+1.44*0.2584
  -0.1705-1.44*0.4237;-0.1705+1.44*0.4237
  -0.4110-1.44*0.1771;-0.4110+1.44*0.1771
  1.0290-1.44*0.1233;1.0290+1.44*0.1233
  0.5254-1.44*0.2156;0.5254+1.44*0.2156
r.squaredGLMM(mod1)
  
 
# a) Females, Chiz?----
FCHgc1 <- FCHgc[FCHgc$ageannee==1,] # 78 ind and obs
unique(factor(FCHgc1$numind))

m1gc <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F), data=FCHgc1)
summary(m1gc)

qqnorm(resid(m1gc))
qqline(resid(m1gc))

options(na.action="na.fail")
s1fch1 <- dredge(m1gc, rank="AICc", evaluate=T)
s1fch1
write.table(s1fch1, "results/model selection tables/shortFCH1.csv", row.names = F, sep=";", dec=".")

dd=subset(s1, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

# Retained model
mod1 <- lm(mass_jul ~ scale(qualite_cohorte, scale=F), data=FCHgc1)
summary(mod1)
12.1282-1.44*0.2088;12.1282+1.44*0.2088
0.8251-1.44*0.2386;0.8251+1.44*0.2386

# R2m et R2c
r.squaredGLMM(mod1)

# Normality best model
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))


# b) Males, Chiz? ----
MCHgc1 <- MCHgc[MCHgc$ageannee==1,]
MCHgc1 <- MCHgc1[!MCHgc1$idkit=="2014-79-031",] # removing the only individuals with two captures at 1 yo
unique(factor(MCHgc1$numind)) # 96 ind/obs

m2gc <- lm(mass_jul ~ scale(qualite_cohorte, scale=F)*scale(log(FGMngg), scale=F), data=MCHgc1)
summary(m2gc)

qqnorm(resid(m2gc))
qqline(resid(m2gc))

options(na.action="na.fail")
s1mch1 <- dredge(m2gc, rank="AICc", evaluate=T)
s1mch1
write.table(s1mch1, "results/model selection tables/shortMCH1.csv", row.names = F, sep=";", dec=".")

dd=subset(s1mch1, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

### mod?le averaging avec tous les mod?les gard?s par nesting rule EG 26/01/2021
avgm <- model.avg(dd,subset=c("8", "3", "4"), fit=T)
avgm#[["coefficients"]]
summary(avgm)
confint(avgm, level=0.85)


avgm <- model.avg(dd,subset=c("8", "4"), fit=T)
avgm#[["coefficients"]]
summary(avgm)
confint(avgm, level=0.95)

MCHgc1$logfgm <- c(scale(log(MCHgc1$FGMngg), scale=F))
MCHgc1$qualan <- c(scale(MCHgc1$qualite_an, scale=F))
mod1 <- lm(mass_jul ~ scale(qualite_an, scale=F)*scale(log(FGMngg), scale=F), data=MCHgc1)
summary(mod1)
mod1a <- lm(mass_jul ~ logfgm*qualan, data=MCHgc1)
summary(mod1a)
# fgm coeff < 0: When fgm increases, then bm decreases
# Qcoh coeff > 0: When env qual increases, then bm increases
# interaction coeff > 0: When Qcoh increases, the slope of fgm against bm increases (goes from <0 towards 0 in this case)
par(mfrow=c(1,1))
visreg(mod1a, "logfgm", by="qualan", overlay=TRUE)
visreg(mod1a, "qualan", by="logfgm", overlay=TRUE)
# R2m et R2c
r.squaredGLMM(mod1)


# Normalité du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)


# c) Females, Trois-Fontaines ----
FTFgc1 <- FTFgc[FTFgc$ageannee==1,]
unique(factor(FTFgc1$numind)) # 92 ind/obs

m3gc <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F), data=FTFgc1)
summary(m3gc)

qqnorm(resid(m3gc))
qqline(resid(m3gc))

options(na.action="na.fail")
s1ftf1 <- dredge(m3gc, rank="AICc", evaluate=T)
write.table(s1ftf1, "results/model selection tables/shortFTF1.csv", row.names = F, sep=";", dec=".")

dd=subset(s1ftf1, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

mod1 <- lm(mass_jul ~ scale(qualite_cohorte, scale=F), data=FTFgc1)
summary(mod1)
15.4950-1.44*0.2355;15.4950+1.44*0.2355
1.3952-1.44*0.3270;1.3952+1.44*0.3270

# R2m et R2c
r.squaredGLMM(mod1)

# Normalité du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# Without the observations with the out of range FGM (6153ngg)
FTFgc1a <- FTFgc1[FTFgc1$FGMngg<=6000,] 
unique(factor(FTFgc1a$numind)) # 91 ind/obs
# FTFgc1a <- FTFgc1a[!FTFgc1a$idkit=="2021-51-005",]
m3gc1a <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F), data=FTFgc1a)

options(na.action="na.fail")
s1ftf1high <- dredge(m3gc1a, rank="AICc", evaluate=T)
write.table(s1ftf1high, "results/model selection tables/shortFTF1high.csv", row.names = F, sep=";", dec=".")

mod1 <- lm(mass_jul ~ scale(qualite_cohorte, scale=F), data=FTFgc1a)
summary(mod1)
15.4969-1.44*0.2381;15.4969+1.44*0.2381
1.3961-1.44*0.3289;1.3961+1.44*0.3289

# R2m et R2c
r.squaredGLMM(mod1)

# Normalité du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# d) Males, Trois-Fontaines ----
MTFgc1 <- MTFgc[MTFgc$ageannee==1,] # 103 obs and ind
unique(factor(MTFgc1$numind))

m4gc <- lm(mass_jul ~ scale(log(FGMngg), scale=F)*scale(qualite_an, scale=F), data=MTFgc1)
summary(m4gc)

qqnorm(resid(m4gc))
qqline(resid(m4gc))

options(na.action="na.fail")
s1mtf1 <- dredge(m4gc, rank="AICc", evaluate=T)
write.table(s1mtf1, "results/model selection tables/shortMTF1.csv", row.names = F, sep=";", dec=".")

# Retained model
mod2 <- lm(mass_jul ~ scale(qualite_an, scale=F) + scale(log(FGMngg), scale=F), data=MTFgc1)
summary(mod2)
15.6556-1.44*0.2075;15.6556+1.44*0.2075
1.0735-1.44*0.2020;1.0735+1.44*0.2020
-0.9909-1.44*0.4009;-0.9909+1.44*0.4009


# R2m et R2c
r.squaredGLMM(mod2)

# Normalité du meilleur modele
hist(residuals(mod2))
qqnorm(residuals(mod2))
qqline(residuals(mod2))
shapiro.test(residuals(mod2))
plot(mod2,ask=TRUE)


# B1B./ 2+ year ----
# 0) Pooling
datagc2a <- datagc[datagc$ageannee>=2,]# 656 obs / 378 ind
datagc2 <- datagc2a[datagc2a$FGMngg<5000,] 
unique(factor(datagc2$numind)) # 655 obs / 377 ind

mod <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F) + scale(log(FGMngg), scale=F)*pop + scale(log(FGMngg), scale=F)*sexe + (1|numind), data=datagc2a, REML=T)
summary(mod)

qqnorm(resid(mod))
qqline(resid(mod))

options(na.action="na.fail")
s1pool2high <- dredge(mod, rank="AICc", evaluate=T)
s1pool2high
write.table(s1pool2high, "results/model selection tables/shortpool2high.csv", row.names = F, sep=";", dec=".")

dd=subset(s1pool2, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

### mod?le averaging avec tous les mod?les gard?s par nesting rule EG 26/01/2021
avgm <- model.avg(dd,subset=c("7", "3"), fit=T)
summary(avgm)
confint(avgm, level=0.85)

mod1 <- lmer(dev ~ scale(rel_cohort_quality, scale=F) + scale(log(FGMngg), scale=F) + (1|numind) + (1|cohorte), data=datagc2)
summary(mod1)

r.squaredGLMM(mod1)


# a) Females, Chiz?----
FCHgc2 <- FCHgc[!FCHgc$age==1,]
unique(factor(FCHgc2$numind))
# 219 obs / 105 ind

m1gc <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc2)
summary(m1gc)
qqnorm(resid(m1gc))
qqline(resid(m1gc))

options(na.action="na.fail")
s1fch2high <- dredge(m1gc, rank="AICc", evaluate=T)
s1fch2high
write.table(s1fch2high, "results/model selection tables/shortFCH2high.csv", row.names = F, sep=";", dec=".")

dd=subset(s1, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

mod1 <- lmer(dev ~ scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc2)
summary(mod1)
0.08602 - 1.44*0.1622;0.08602 + 1.44*0.1622
0.59953 - 1.44*0.23924;0.59953 + 1.44*0.23924

# R2m et R2c
r.squaredGLMM(mod1)

# Normalité du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# Normalit? effet al?atoire
plot(ranef(mod1))[1]

# to remove the second high value of FGM (5192 ngg)
FCHgc2a <- FCHgc2[FCHgc2$FGMngg<=5000,] 
unique(factor(FCHgc2a$numind))

m1gca <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc2a)
options(na.action="na.fail")
s1fch2 <- dredge(m1gca, rank="AICc", evaluate=T)
s1fch2
write.table(s1fch2, "results/model selection tables/shortFCH2.csv", row.names = F, sep=";", dec=".")

dd=subset(s1, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

mod1 <- lmer(dev ~ scale(qualite_cohorte, scale=F) + (1|numind), data=FCHgc2a)
summary(mod1)
0.1055 - 1.44*0.1623;0.1055 + 1.44*0.1623
0.5988 - 1.44*0.2385;0.5988 + 1.44*0.2385

# R2m et R2c
r.squaredGLMM(mod1)

# Normalit? du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# Normalit? effet al?atoire
plot(ranef(mod1))[1]

# b) Males, Chiz? ----
MCHgc2 <- MCHgc[!MCHgc$age==1,]
unique(factor(MCHgc2$numind))
# 136 obs / 83 ind

m2gc <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=MCHgc2)
summary(m2gc)
qqnorm(resid(m2gc))
qqline(resid(m2gc))

options(na.action="na.fail")
s1mch2 <- dredge(m2gc, rank="AICc", evaluate=T)
s1mch2
write.table(s1mch2, "results/model selection tables/shortMCH2.csv", row.names = F, sep=";", dec=".")

dd=subset(s1mch2, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

# Mod?le averaging #
mod1 <- lmer(dev ~ scale(log(FGMngg), scale=F) + (1|numind), data=MCHgc2)
summary(mod1)
0.01705-1.44*0.20417;0.01705+1.44*0.20417
-0.52086-1.44*0.17243;-0.52086+1.44*0.17243

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

# c) Females, Trois-Fontaines ----
FTFgc2 <- FTFgc[!FTFgc$age==1,]
unique(factor(FTFgc2$numind))
# 164 obs / 106 ind
# FTFgc2 <- FTFgc2[FTFgc2$FGMngg>10,] to test without the two very low FGM (8 and 9 ng/g)
# results are similar without those

m3gc <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=FTFgc2)
summary(m3gc)
qqnorm(resid(m3gc))
qqline(resid(m3gc))

options(na.action="na.fail")
s1ftf2 <- dredge(m3gc, rank="AICc", evaluate=T)
s1ftf2
write.table(s1ftf2, "results/model selection tables/shortFTF2.csv", row.names = F, sep=";", dec=".")

dd=subset(s1ftf2, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

mod1 <- lmer(dev ~ (1|numind), data=FTFgc2)
summary(mod1)
-0.2076 - 1.44*0.2324;-0.2076 + 1.44*0.2324


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

# d) Males, Trois-Fontaines ----
MTFgc2 <- MTFgc[!MTFgc$age==1,]
unique(factor(MTFgc2$numind))
# 137 obs / 84 ind

m4gc <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=MTFgc2)
summary(m4gc)
qqnorm(resid(m4gc))
qqline(resid(m4gc))

options(na.action="na.fail")
s1mtf2 <- dredge(m4gc, rank="AICc", evaluate=T)
s1mtf2
write.table(s1mtf2, "results/model selection tables/shortMTF2.csv", row.names = F, sep=";", dec=".")

dd=subset(s1mtf2, delta<2)
dd

nst=nested(dd)
nst
d=subset(dd, !nested(.))
d

# Retained model
mod1 <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(qualite_cohorte, scale=F) + (1|numind), data=MTFgc2)
summary(mod1)
0.00521-1.44*0.25283;0.00521+1.44*0.25283
-0.91944-1.44*0.22162;-0.91944+1.44*0.22162
0.07626-1.44*0.2753;0.07626+1.44*0.2753
1.60476-1.44*0.35053;1.60476+1.44*0.35053

# R2m et R2c
r.squaredGLMM(mod1)

# Normalit? du meilleur modele
hist(residuals(mod1))
qqnorm(residuals(mod1))
qqline(residuals(mod1))
shapiro.test(residuals(mod1))
plot(mod1,ask=TRUE)

# Normalit? des effets al?atoires
plot(ranef(mod1))[1]



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
  unique(factor(data$numind)) # 1881 obs / 947 ind
  
  
  # Calculate a relative environmental quality per population
  data$annee <- as.factor(data$annee)
  data$mean_cohort_quality <- ave(data$qualite_cohorte, data$pop, FUN=function(x) mean(x, na.rm=T)) # mean cohort quality in CH and TF
  data$rel_cohort_quality <- data$qualite_cohorte - data$mean_cohort_quality # relative cohort quality
  

  datagc <- data[!is.na(data$FGMngg),] # 1026
  datagc$age_class <- ifelse(datagc$ageannee==1, "juveniles", "adults")
}

# 1) Populations and sexes pooled 
# 2) Juveniles and adults pooled
plot(log(datagc$FGMngg) ~ datagc$age)
plot(datagc$FGMngg ~ datagc$age)
# 3) Add year of capture in random structure
datagc$cohorte <- as.factor(datagc$cohorte)
datagc$annee <- as.factor(datagc$annee)
# 4) Compute the full model and just look at the estimates

datagcA <- datagc[datagc$FGMngg<5000,] 
unique(factor(datagc$numind)) # 1024 obs / 629 ind

# WE CANNOT COMPUTE THE FULL MODEL FGM*COHORT*POP*SEX*AGECLASS
# BECAUSE MUMIN CANNOT PROCESS MODELS WITH MORE THAN 30 PARAMETERS

mod <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F)*pop*sexe + 
                  (1|numind) + (1|annee), data=datagcA, REML=F)
summary(mod)

options(na.action="na.fail")
dr <- dredge(mod, rank="AICc", evaluate=T)

# SO WE SPLIT THE DATASET BY AGE CLASS (JUVENILES V. ADULTS)

# JUVENILES
datagcA1 <- subset(datagcA, ageannee==1) 
datagcA1 <- datagcA1[!datagcA1$idkit=="2014-79-031",] # removing the only individuals with two captures at 1 yo # 369 obs/ind
unique(factor(datagcA1$numind)) # 448 obs / 448 ind

mod1 <- lm(dev ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F)*pop*sexe, data=datagcA1, REML=F)
summary(mod1)

options(na.action="na.fail")
dr1 <- dredge(mod1, rank="AICc", evaluate=T, trace=2)
ddr1 <- subset(dr1, delta<=2)
ddr1
d1 <- subset(ddr1, !nested(.))
d1

write.table(dr1, "results/Review OIKOS/results_short-term_juv.csv", row.names = F, sep=";", dec=".")

datagcA2 <- subset(datagcA, ageannee>=2) 
unique(factor(datagcA2$numind)) # 655 obs / 377 ind

mod2 <- lmer(dev ~ scale(log(FGMngg), scale=F)*scale(rel_cohort_quality, scale=F)*pop*sexe + 
               (1|numind) + (1|annee), data=datagcA2, REML=F)
summary(mod2)

options(na.action="na.fail")
dr2 <- dredge(mod2, rank="AICc", evaluate=T, trace=2)
ddr2 <- subset(dr2, delta<=2)
ddr2
d2 <- subset(ddr2, !nested(.))
d2

write.table(dr2, "results/Review OIKOS/results_short-term_adults.csv", row.names = F, sep=";", dec=".")

