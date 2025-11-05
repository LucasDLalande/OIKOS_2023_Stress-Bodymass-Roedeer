############################# ANALYSE STRESS-MASS ##############################
################################ Repeatability #################################
################# A./ POPULATION-SPECIFIC STRESS---------------#################
#################   1) FGM-------------------------------------#################
#################   2) NLR-------------------------------------#################
################# B./ REPEATABILITY----------------------------#################
#################   1) FGM-------------------------------------#################
#################   2) NLR-------------------------------------#################
################################# 24/05/2022 ###################################
rm(list=ls())

library(lme4)
library(MuMIn)
library(ggplot2)
library(rptR)
library(lubridate)
library(dplyr)


{
data <- read.csv("dataset/dataset-20220516.csv", h=T, sep=";", dec=".")
data <- data[!is.na(data$ageannee),]
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
data$annee <- as.factor(data$annee)
data$date <- dmy(data$date)

# FGM of Priority = 3 --> NA
data <- data %>%
        mutate((replace(FGMngg, FGM_Priority == "3", NA)))
data <- data[,-18]
names(data)[18] <- "FGMngg"

rownames(data) <- 1:nrow(data)
str(data)



dataMCH <- data[data$pop=="C" & data$sexe=="M",]
dataMCH <- dataMCH[dataMCH$ageannee<=10,]
dataFCH <- data[data$pop=="C" & data$sexe=="F",]
dataFCH <- dataFCH[dataFCH$ageannee<=15,]
dataMTF <- data[data$pop=="3F" & data$sexe=="M",]
dataMTF <- dataMTF[dataMTF$ageannee<=12,]
dataFTF <- data[data$pop=="3F" & data$sexe=="F",]
dataFTF <- dataFTF[dataFTF$ageannee<=13,]

data <- rbind(dataMCH, dataFCH, dataMTF, dataFTF)
tapply(data$ageannee, data$pop:data$sexe, table)

dataC <- data[data$pop=="C",] # 917 masses
dataTF <- data[data$pop=="3F",] # 964 masses
}

# A./ Population-specific Stress ----
# FGMs ----
datagc <- data[!is.na(data$FGMngg),] # 1026

datagcMCH <- datagc[datagc$pop=="C" & datagc$sexe=="M",]
#datagcMCH <- datagcMCH[datagcMCH$ageannee<=10,]
datagcFCH <- datagc[datagc$pop=="C" & datagc$sexe=="F",]
#datagcFCH <- datagcFCH[datagcFCH$ageannee<=15,]
datagcMTF <- datagc[datagc$pop=="3F" & datagc$sexe=="M",]
#datagcMTF <- datagcMTF[datagcMTF$ageannee<=12,]
datagcFTF <- datagc[datagc$pop=="3F" & datagc$sexe=="F",]
#datagcFTF <- datagcFTF[datagcFTF$ageannee<=13,]

datagca <- rbind(datagcMCH, datagcFCH, datagcMTF, datagcFTF)
tapply(datagc$ageannee, datagc$pop:datagc$sexe, table)
datagc <- datagca[datagca$FGMngg < 5000,] # 1024


mean(datagc$FGMngg) # 749.2 ng/g
sd(datagc$FGMngg) # 431.1
range(datagc$FGMngg) # 8.0-3427.6
median(datagc$FGMngg) # 689.0

mean(datagc$FGMngg[datagc$pop=="C"]) # 732.9
sd(datagc$FGMngg[datagc$pop=="C"]) # 437.3
range(datagc$FGMngg[datagc$pop=="C"]) # 34.2-3275.3
median(datagc$FGMngg[datagc$pop=="C"]) # 674.0

mean(datagc$FGMngg[datagc$pop=="3F"]) # 766.7
sd(datagc$FGMngg[datagc$pop=="3F"]) # 424.0
range(datagc$FGMngg[datagc$pop=="3F"]) # 8-3427.6
median(datagc$FGMngg[datagc$pop=="3F"]) # 715.0

mean(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"]) # 726.6
sd(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"]) # 466.3
range(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"]) # 64.9-3275.3
median(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"]) # 645.0
shapiro.test(log(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"]))

mean(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"]) # 737.9
sd(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"]) # 413.8
range(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"]) # 34.2-2488.1
median(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"]) # 685.5
shapiro.test(log(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"]))

wilcox.test(datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="F"], datagc$FGMngg[datagc$pop=="C" & datagc$sexe=="M"])

mean(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"]) # 765.9
sd(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"]) # 393.5
range(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"]) # 8.9-3427.6
median(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"]) # 741.7
shapiro.test(log(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"]))

mean(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"]) # 767.4
sd(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"]) # 451.6
range(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"]) # 8-3369.9
median(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"]) # 683.5
shapiro.test(log(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"]))

wilcox.test(datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="F"], datagc$FGMngg[datagc$pop=="3F" & datagc$sexe=="M"])

tfgc <- datagc[datagc$pop=="3F",]
tfgc$numind <- factor(tfgc$numind)

chgc <- datagc[datagc$pop=="C",]
chgc$numind <- factor(chgc$numind)

hist(tfgc$FGMngg)
hist(chgc$FGMngg)
shapiro.test(log(tfgc$FGMngg))
shapiro.test(log(chgc$FGMngg))
qqnorm(log(tfgc$FGMngg))
qqline(log(tfgc$FGMngg))
qqnorm(log(chgc$FGMngg))
qqline(log(chgc$FGMngg))

wilcox.test(tfgc$FGMngg, chgc$FGMngg)
wilcox.test(log(tfgc$FGMngg), log(chgc$FGMngg))
levels(datagc$pop)
datagc$pop <- factor(datagc$pop, levels = c("C", "3F"))

library(wesanderson)
png(file="results/figures/histogram_fgm.png")
ggplot(datagc, aes(x=FGMngg, fill=pop)) + 
        geom_density(size=1, alpha=0.3) +
        #scale_fill_manual(values=wes_palette(n=2, name="Darjeeling1")) +
        scale_fill_discrete(name="Population", labels=c("Chizé", "Trois-Fontaines")) +
        ylab("Density") +
        xlab("Faecal Glucocorticoid Metabolites (FGM)") +
        theme(plot.margin=margin(0,1,0,0, "cm"),
              plot.title = element_text(face = "bold", hjust=0.5, size=18), # plot title
              axis.title.x = element_text(size=14), #change size x axis title
              axis.text.x = element_text(size=11.5), #change size text x axis
              axis.title.y = element_text(size=14), #change size y axis title
              axis.text.y = element_text(size=11.5), #change size text y axis
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # no grid in the background
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=12),#change legend text font size
              legend.position=c(0.8, 0.9),
              legend.background = element_rect(fill='transparent'))
dev.off()

png(file="results/figures/histogram_highfgm.png")
ggplot(datagca, aes(x=FGMngg, fill=pop)) + 
  geom_density(size=1, alpha=0.3) +
  #scale_fill_manual(values=wes_palette(n=2, name="Darjeeling1")) +
  scale_fill_discrete(name="Population", labels=c("Chizé", "Trois-Fontaines")) +
  ylab("Density") +
  xlab("Faecal Glucocorticoid Metabolites (FGM)") +
theme(plot.margin=margin(0,1,0,0, "cm"),
      plot.title = element_text(face = "bold", hjust=0.5, size=18), # plot title
      axis.title.x = element_text(size=14), #change size x axis title
      axis.text.x = element_text(size=11.5), #change size text x axis
      axis.title.y = element_text(size=14), #change size y axis title
      axis.text.y = element_text(size=11.5), #change size text y axis
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # no grid in the background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title = element_text(size=14), #change legend title font size
      legend.text = element_text(size=12),#change legend text font size
      legend.position=c(0.8, 0.9),
      legend.background = element_rect(fill='transparent'))
dev.off()

plot(datagc$ageannee, log(datagc$FGMngg))

var(log(datagc$FGMngg[datagc$sexe=="M"])) 
var(log(datagc$FGMngg[datagc$sexe=="F"]))
shapiro.test(log(datagc$FGMngg[datagc$sexe=="M"]))
shapiro.test(log(datagc$FGMngg[datagc$sexe=="F"]))
var(datagc$FGMngg[datagc$sexe=="M"])
var(datagc$FGMngg[datagc$sexe=="F"])
var.test(log(datagc$FGMngg[datagc$pop=="C"]), log(datagc$FGMngg[datagc$pop=="3F"]))

library(onewaytests)
bf.test(log(FGMngg) ~ sexe, datagc)
bf.test(log(FGMngg) ~ pop, datagc)




# MASS ----
datam <- data[!is.na(data$masse),] # 1772
mean(datam$masse) # 19.4 kg

datam1 <- datam[datam$ageannee==1,]
datam1 <- datam1[!is.na(datam1$ageannee),]
datam2 <- datam[datam$ageannee>=2,]
datam2 <- datam2[!is.na(datam2$ageannee),]

library(wesanderson)
ggplot(datam1, aes(x=masse, fill=pop)) + 
        geom_density(size=1, alpha=0.4) +
        scale_fill_manual(values=wes_palette(n=2, name="Moonrise2"))

ggplot(datam2, aes(x=masse, fill=pop)) + 
        geom_density(size=1, alpha=0.4) +
        scale_fill_manual(values=wes_palette(n=2, name="Moonrise2"))

ggplot(datam, aes(x=masse, fill=pop)) + 
        geom_density(size=1, alpha=0.4) +
        scale_fill_manual(values=wes_palette(n=2, name="Moonrise2"))

mean(datam1$masse)
with(datam1, tapply(masse, pop, mean))

12.58819 - qnorm(0.975)*sd(datam1$masse[datam1$pop=="C"])/sqrt(length(datam1$masse[datam1$pop=="C"]))
12.58819 + qnorm(0.975)*sd(datam1$masse[datam1$pop=="C"])/sqrt(length(datam1$masse[datam1$pop=="C"]))
15.69633 - qnorm(0.975)*sd(datam1$masse[datam1$pop=="3F"])/sqrt(length(datam1$masse[datam1$pop=="3F"]))
15.69633 + qnorm(0.975)*sd(datam1$masse[datam1$pop=="3F"])/sqrt(length(datam1$masse[datam1$pop=="3F"]))

mean(datam2$masse)
with(datam2, tapply(masse, pop, mean))

20.49535 - qnorm(0.975)*sd(datam2$masse[datam2$pop=="C"])/sqrt(length(datam2$masse[datam2$pop=="C"]))
20.49535 + qnorm(0.975)*sd(datam2$masse[datam2$pop=="C"])/sqrt(length(datam2$masse[datam2$pop=="C"]))
23.18170 - qnorm(0.975)*sd(datam2$masse[datam2$pop=="3F"])/sqrt(length(datam2$masse[datam2$pop=="3F"]))
23.18170 + qnorm(0.975)*sd(datam2$masse[datam2$pop=="3F"])/sqrt(length(datam2$masse[datam2$pop=="3F"]))

mean(datam$masse)
with(datam, tapply(masse, pop, mean))

18.15702 - qnorm(0.975)*sd(datam$masse[datam$pop=="C"])/sqrt(length(datam$masse[datam$pop=="C"]))
18.15702 + qnorm(0.975)*sd(datam$masse[datam$pop=="C"])/sqrt(length(datam$masse[datam$pop=="C"]))
20.49576 - qnorm(0.975)*sd(datam$masse[datam$pop=="3F"])/sqrt(length(datam$masse[datam$pop=="3F"]))
20.49576 + qnorm(0.975)*sd(datam$masse[datam$pop=="3F"])/sqrt(length(datam$masse[datam$pop=="3F"]))


# B./ Repeatability ----
# 1) FGMs ----
datagc <- datagc[!is.na(datagc$numind),] # 1711

# Seems better to assess repeatability of log(FGMs) as residuals of the LMM are closer to a normal error distribution 
rpt(log(FGMngg) ~ annee + (1|numind), grname="numind", datatype="Gaussian", npermut=1000, data=datagc)
summary(lmer(log(FGMngg) ~ annee + (1|numind), data=datagc))
0.06402/(0.06402+0.38189)
plot(lmer(log(FGMngg) ~ (1|numind), data=datagc))

tfgc <- datagc[datagc$pop=="3F",]
tfgc$numind <- factor(tfgc$numind)

chgc <- datagc[datagc$pop=="C",]
chgc$numind <- factor(chgc$numind)

rept <- lmer(log(FGMngg) ~ annee + (1|numind), data=tfgc)
summary(rept)
0.06737/(0.06737+0.37335) # 0.153
0.04443/(0.04443+0.34709) # 0.113
plot(rept)

repc <- lmer(log(FGMngg) ~ annee + (1|numind), data=chgc)
summary(repc)
0.06476/(0.06476+0.38456) # 0.144
0.02507/(0.02507+0.30063) # 0.077
plot(repc)

rpt(log(FGMngg) ~ annee + (1|numind), grname="numind", datatype="Gaussian", npermut=1000, data=tfgc)
# 0.113 [0, 0.252]
rpt(log(FGMngg) ~ annee + (1|numind), grname="numind", datatype="Gaussian", npermut=1000, data=chgc)
# 0.077 [0, 0.206]

# Excluding fawns/juveniles
datagc2 <- datagc[datagc$ageannee>=2,]
tfgc2 <- datagc2[datagc2$pop=="3F",]
chgc2 <- datagc2[datagc2$pop=="C",]

rpt(log(FGMngg) ~ (1|numind), grname="numind", datatype="Gaussian", npermut=1000, data=tfgc2)

rpt(log(FGMngg) ~ (1|numind), grname="numind", datatype="Gaussian", npermut=1000, data=chgc2)

# Weak repeatability of the FGMs (log-transformed)