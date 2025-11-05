############################# ANALYSE STRESS-MASS ##############################
################# Mass variation between two consecutive years #################
################# A./ FGM--------------------------------------#################
#################   1) Growth (1-2yo)--------------------------#################
#################     FIGURES----------------------------------#################
#################   2) Stabilised (4-10yo)---------------------#################
#################     FIGURES----------------------------------#################
################################# 01/06/2022 ###################################
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
  names(data)[18] <- "FGMngg"
  
  rownames(data) <- 1:nrow(data)
  str(data)
  
  data$age_factor <- as.factor(data$ageannee)
  
  # Calculate fawn mass corrected for the date of capture (brought back to the 27th of january)
  # Based on Douhard et al., 2017
  dataCH <- data[data$pop=="C",]
  dataCH$masse <- ifelse(dataCH$ageannee==1, dataCH$masse+0.012*(57-dataCH$datejulienne), dataCH$masse)
  dataTF <- data[data$pop=="3F",]
  dataTF$masse <- ifelse(dataTF$ageannee==1, dataTF$masse+0.024*(57-dataTF$datejulienne), dataTF$masse)
  
  data <- rbind(dataCH, dataTF)
}


# A./ FGMs ----
datagc <- data

# 1) Growth (1-2 yo) ----
# a) MCHgc: Males, Chiz? ----
MCH <- datagc[datagc$pop == "C" & datagc$sexe == "M",]
MCHgc <- MCH[MCH$ageannee <= 2,] # only individuals of 1 and 2 yo
MCHgc <- MCHgc[!MCHgc$idkit=="2011-79-050",] # removing an obs for individual MC2237 captured twice at 2yo (keep the first one)
MCHgc <- MCHgc[!MCHgc$idkit=="2014-79-031",] # removing an obs for individual MC2271 captured twice at 1yo (keep the first one)
MCHgc <- MCHgc[!MCHgc$idkit=="2020-79-041",] # removing an obs for individual CC6139 captured twice at 1yo (keep the first one)
unique(factor(MCHgc$numind)) # 196 obs / 163 ind

# Ordering dataset by individuals and measurement date
MCHgc <- MCHgc[with(MCHgc, order(numind, date)),]
rownames(MCHgc) <- 1:nrow(MCHgc)

# Remove last measurement for each individuals in a first dataset
MCHgc1 <- MCHgc %>% group_by(numind) %>% slice(-n())
MCHgc1 <- MCHgc1[,-c(13,14,15,16,17,19)]

# Remove first measurement for each individuals in a second dataset
MCHgc2 <- MCHgc %>% group_by(numind) %>% slice(-1)
MCHgc2 <- MCHgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
MCHgcn <- cbind(MCHgc1, MCHgc2)
MCHgcn <- MCHgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(MCHgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

unique(factor(MCHgcn$id)) # 33 individuals

# Mass change between year 1 and year 2
MCHgcn$mass_change <- MCHgcn$mass2 - MCHgcn$mass1

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
MCHgcn1 <- MCHgcn[!is.na(MCHgcn$FGM1),] # 24 individuals
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
MCHgcn2 <- MCHgcn[!is.na(MCHgcn$FGM2),] # 22 individuals
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
MCHgcmean <- MCHgcn[!is.na(MCHgcn$FGM1),] # 24 individuals
MCHgcmean <- MCHgcmean[!is.na(MCHgcmean$FGM2),] # 16 individuals
MCHgcmean$mean_fgm <- (MCHgcmean$FGM1 + MCHgcmean$FGM2)/2 # mean FGM during two consecutive years

# Mean mass for all individuals of 1y - Mean mass for all individuals of 2y
MCH1 <- MCH[MCH$ageannee==1,]
mean1 <- mean(MCH1$masse)

MCH2 <- MCH[MCH$ageannee==2,]
mean2 <- mean(MCH2$masse)

mass_change_pop <- mean2 - mean1

# Individual mass change between year 1 and 2 - Sample diff in mass between y1 and 2
MCHgcn1$r_mass <- MCHgcn1$mass_change - mass_change_pop # year1
MCHgcn2$r_mass <- MCHgcn2$mass_change - mass_change_pop # year2
MCHgcmean$r_mass <- MCHgcmean$mass_change - mass_change_pop # mean



# b) FCHgc: Females, Chiz? ----
FCH <- datagc[datagc$pop == "C" & datagc$sexe == "F",]
FCHgc <- FCH[FCH$ageannee <= 2,] # only individuals of 1 and 2 yo
FCHgc <- FCHgc[!FCHgc$idkit=="2013-79-077",] # removing an obs for individual FC1202 captured twice at 2yo (keep the first one)
unique(factor(FCHgc$numind)) # 176 obs / 148 ind

# Ordering dataset by individuals and measurement date
FCHgc <- FCHgc[with(FCHgc, order(numind, date)),]
rownames(FCHgc) <- 1:nrow(FCHgc)

# Remove last measurement for each individuals in a first dataset
FCHgc1 <- FCHgc %>% group_by(numind) %>% slice(-n())
FCHgc1 <- FCHgc1[,-c(13,14,15,16,17,19)]

# Remove first measurement for each individuals in a second dataset
FCHgc2 <- FCHgc %>% group_by(numind) %>% slice(-1)
FCHgc2 <- FCHgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
FCHgcn <- cbind(FCHgc1, FCHgc2)
FCHgcn <- FCHgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(FCHgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

unique(factor(FCHgcn$id)) # 28 individuals

# Mass change between year 1 and year 2
FCHgcn$mass_change <- FCHgcn$mass2 - FCHgcn$mass1

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
FCHgcn1 <- FCHgcn[!is.na(FCHgcn$FGM1),] # 19 individuals
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
FCHgcn2 <- FCHgcn[!is.na(FCHgcn$FGM2),] # 18 individuals
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
FCHgcmean <- FCHgcn[!is.na(FCHgcn$FGM1),] # 19 individuals
FCHgcmean <- FCHgcmean[!is.na(FCHgcmean$FGM2),] # 12 individuals
FCHgcmean$mean_fgm <- (FCHgcmean$FGM1 + FCHgcmean$FGM2)/2 # mean FGM during two consecutive years

# Mean mass for all individuals of 2y - Mean mass for all individuals of 2y
FCH1 <- FCH[FCH$ageannee==1,]
mean1 <- mean(FCH1$masse)

FCH2 <- FCH[FCH$ageannee==2,]
mean2 <- mean(FCH2$masse)

mass_change_pop <- mean2 - mean1

# Individual mass change between year 1 and 2 - Sample diff in mass between y1 and 2
FCHgcn1$r_mass <- FCHgcn1$mass_change - mass_change_pop
FCHgcn2$r_mass <- FCHgcn2$mass_change - mass_change_pop
FCHgcmean$r_mass <- FCHgcmean$mass_change - mass_change_pop


# c) MTFgc: Males, Trois-Fontaine ----
MTF <- datagc[datagc$pop == "3F" & datagc$sexe == "M",]
MTFgc <- MTF[MTF$ageannee <= 2,] # only individuals of 1 and 2 yo
unique(factor(MTFgc$numind)) # 243 obs / 206 ind

# Ordering dataset by individuals and measurement date
MTFgc <- MTFgc[with(MTFgc, order(numind, date)),]
rownames(MTFgc) <- 1:nrow(MTFgc)

# Remove last measurement for each individuals in a first dataset
MTFgc1 <- MTFgc %>% group_by(numind) %>% slice(-n())
MTFgc1 <- MTFgc1[,-c(13,14,15,16,17,19)]
unique(factor(MTFgc1$numind))

# Remove first measurement for each individuals in a second dataset
MTFgc2 <- MTFgc %>% group_by(numind) %>% slice(-1)
unique(factor(MTFgc2$numind))
MTFgc2 <- MTFgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
MTFgcn <- cbind(MTFgc1, MTFgc2)
MTFgcn <- MTFgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(MTFgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

unique(factor(MTFgcn$id)) # 37 individuals

# Mass change between year 1 and year 2
MTFgcn$mass_change <- MTFgcn$mass2 - MTFgcn$mass1

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
MTFgcn1 <- MTFgcn[!is.na(MTFgcn$FGM1),] # 28 individuals
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
MTFgcn2 <- MTFgcn[!is.na(MTFgcn$FGM2),] # 23 individuals
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
MTFgcmean <- MTFgcn[!is.na(MTFgcn$FGM1),] # 25 individuals
MTFgcmean <- MTFgcmean[!is.na(MTFgcmean$FGM2),] # 19 individuals
MTFgcmean$mean_fgm <- (MTFgcmean$FGM1 + MTFgcmean$FGM2)/2 # mean FGM during two consecutive years

# Mean mass for all individuals of 2y - Mean mass for all individuals of 2y
MTF1 <- MTF[MTF$ageannee==1,]
mean1 <- mean(MTF1$masse)

MTF2 <- MTF[MTF$ageannee==2,]
mean2 <- mean(MTF2$masse)

mass_change_pop <- mean2 - mean1

# Individual mass change between year 1 and 2 - Sample diff in mass between y1 and 2
MTFgcn1$r_mass <- MTFgcn1$mass_change - mass_change_pop
MTFgcn2$r_mass <- MTFgcn2$mass_change - mass_change_pop
MTFgcmean$r_mass <- MTFgcmean$mass_change - mass_change_pop


# d) FTFgc: Females, Trois-Fontaine ----
FTF <- datagc[datagc$pop == "3F" & datagc$sexe == "F",]
FTFgc <- FTF[FTF$ageannee <= 2,] # only individuals of 1 and 2 yo
FTFgc <- FTFgc[!FTFgc$idkit=="2010-51-043",] # removing an obs for individual FT357 captured twice at 2yo (keep the first one)
unique(factor(FTFgc$numind)) # 239 obs / 198 ind

# Ordering dataset by individuals and measurement date
FTFgc <- FTFgc[with(FTFgc, order(numind, date)),]
rownames(FTFgc) <- 1:nrow(FTFgc)

# Remove last measurement for each individuals in a first dataset
FTFgc1 <- FTFgc %>% group_by(numind) %>% slice(-n())
FTFgc1 <- FTFgc1[,-c(13,14,15,16,17,19)]
unique(factor(FTFgc1$numind))

# Remove first measurement for each individuals in a second dataset
FTFgc2 <- FTFgc %>% group_by(numind) %>% slice(-1)
unique(factor(FTFgc2$numind))
FTFgc2 <- FTFgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
FTFgcn <- cbind(FTFgc1, FTFgc2)
FTFgcn <- FTFgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(FTFgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

unique(factor(FTFgcn$id)) # 41 individuals

# Mass change between year 1 and year 2
FTFgcn$mass_change <- FTFgcn$mass2 - FTFgcn$mass1

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
FTFgcn1 <- FTFgcn[!is.na(FTFgcn$FGM1),] # 29 individuals
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
FTFgcn2 <- FTFgcn[!is.na(FTFgcn$FGM2),] # 26 individuals
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
FTFgcmean <- FTFgcn[!is.na(FTFgcn$FGM1),] # 29 individuals
FTFgcmean <- FTFgcmean[!is.na(FTFgcmean$FGM2),] # 21 individuals
FTFgcmean$mean_fgm <- (FTFgcmean$FGM1 + FTFgcmean$FGM2)/2 # mean FGM during two consecutive years

# Mean mass for all individuals of 2y - Mean mass for all individuals of 2y
FTF1 <- FTF[FTF$ageannee==1,]
mean1 <- mean(FTF1$masse)

FTF2 <- FTF[FTF$ageannee==2,]
mean2 <- mean(FTF2$masse)

mass_change_pop <- mean2 - mean1

# Individual mass change between year 1 and 2 - Sample diff in mass between y1 and 2
FTFgcn1$r_mass <- FTFgcn1$mass_change - mass_change_pop
FTFgcn2$r_mass <- FTFgcn2$mass_change - mass_change_pop
FTFgcmean$r_mass <- FTFgcmean$mass_change - mass_change_pop


# MODELS ----
# Pool the 4 subsets:
datagcn1a <- rbind(FCHgcn1, MCHgcn1, FTFgcn1, MTFgcn1)
unique(factor(datagcn1a$id)) # 100 individuals
datagcn1 <- datagcn1a[datagcn1a$FGM1<5000,]
unique(factor(datagcn1$id)) # 99 individuals
datagcn2 <- rbind(FCHgcn2, MCHgcn2, FTFgcn2, MTFgcn2)
unique(factor(datagcn2$id)) # 89 individuals
datagcmeana <- rbind(FCHgcmean, MCHgcmean, FTFgcmean, MTFgcmean)
unique(factor(datagcmeana$id)) # 68 individuals
datagcmean <- datagcmeana[datagcmeana$FGM1<5000,]
unique(factor(datagcmean$id)) # 67 individuals

plot(r_mass ~ qualite2, data=datagcn1)
plot(r_mass ~ qualite2, data=datagcn2)
plot(r_mass ~ qualite2, data=datagcmean)

# FIGURES ----
# Fig2A - Growing individuals
datagcn1$logFGM1 <- log(datagcn1$FGM1) 
moda1 <- lmer(r_mass ~ scale(logFGM1, scale=F) + (1|cohorte), data = datagcn1)
summary(moda1)

newdataa1=expand.grid(logFGM1=seq(min(log(datagcn1$FGM1)),
                                  max(log(datagcn1$FGM1)), 0.005))

preda1=predictSE(moda1,newdata=newdataa1,se.fit=TRUE)
newdataa1=cbind(newdataa1,preda1)
newdataa1$low=newdataa1$fit-1.44*newdataa1$se.fit
newdataa1$upp=newdataa1$fit+1.44*newdataa1$se.fit

datagcmean$logmeanfgm <- log(datagcmean$mean_fgm) 
moda2 <- lmer(r_mass ~ scale(logmeanfgm, scale=F) + (1|cohorte), data = datagcmean)
summary(moda2)

newdataa2=expand.grid(logmeanfgm=seq(min(log(datagcmean$mean_fgm)),
                                  max(log(datagcmean$mean_fgm)), 0.005))

preda2=predictSE(moda2,newdata=newdataa2,se.fit=TRUE)
newdataa2=cbind(newdataa2,preda2)
newdataa2$low=newdataa2$fit-1.44*newdataa2$se.fit
newdataa2$upp=newdataa2$fit+1.44*newdataa2$se.fit

show_col(pal_futurama("planetexpress", alpha = 1)(12))

groba <- grobTree(textGrob("Growing individuals", x=0.05,  y=0.97, hjust=0,
                           gp=gpar(col="black", fontsize=14, face="bold")))
fig2a <- ggplot(data = datagcn1, aes(y=r_mass, x=logFGM1)) + 
  geom_point(data = datagcn1, aes(y=r_mass, x=logFGM1, color="#5FB55F")) +
  geom_point(data = datagcmean, aes(y=r_mass, x=logmeanfgm, color="#D35FB7")) +
  scale_y_continuous(name = "Change in relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
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
  scale_color_manual(values=c("#D35FB7", "#5FB55F"), 
                        name="FGM measurement",
                        labels=c("1st year", "Mean 1st and 2nd year"))
fig2a


  

# 2) Stabilised (4-10 yo) ----
# a) MCHgc: Males, Chiz? ----
MCH <- datagc[datagc$pop == "C" & datagc$sexe == "M",]
MCHgc <- MCH[MCH$ageannee >= 4 & MCH$ageannee <= 10 ,] # only individuals between 4 and 10 yo
unique(factor(MCHgc$numind)) # 148 obs / 95 ind

# Ordering dataset by individuals and measurement date
MCHgc <- MCHgc[with(MCHgc, order(numind, date)),]
rownames(MCHgc) <- 1:nrow(MCHgc)

# Remove last measurement for each individuals in a first dataset
MCHgc1 <- MCHgc %>% group_by(numind) %>% slice(-n())
MCHgc1 <- MCHgc1[,-c(13,14,15,16,17,19)]

# Remove first measurement for each individuals in a second dataset
MCHgc2 <- MCHgc %>% group_by(numind) %>% slice(-1)
MCHgc2 <- MCHgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
MCHgcn <- cbind(MCHgc1, MCHgc2)
MCHgcn <- MCHgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(MCHgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

# Remove measures with interval >1 year
MCHgcn <- MCHgcn[MCHgcn$age2-MCHgcn$age1==1,]

# Mass change between year 1 and year 2
MCHgcn$mass_change <- MCHgcn$mass2 - MCHgcn$mass1

# Mean mass for all individuals of year1 - Mean mass for all individuals of year2
MCH4 <- MCH[MCH$ageannee==4,]
mean4 <- mean(MCH4$masse)
MCH5 <- MCH[MCH$ageannee==5,]
mean5 <- mean(MCH5$masse)
mean45 <- mean5-mean4
MCH6 <- MCH[MCH$ageannee==6,]
mean6 <- mean(MCH6$masse)
mean56 <- mean6-mean5
MCH7 <- MCH[MCH$ageannee==7,]
mean7 <- mean(MCH7$masse)
mean67 <- mean7-mean6
MCH8 <- MCH[MCH$ageannee==8,]
mean8 <- mean(MCH8$masse)
mean78 <- mean8-mean7
MCH9 <- MCH[MCH$ageannee==9,]
mean9 <- mean(MCH9$masse)
mean89 <- mean9-mean8
MCH10 <- MCH[MCH$ageannee==10,]
mean10 <- mean(MCH10$masse)
mean910 <- mean10-mean9

MCHgcn4 <- MCHgcn[MCHgcn$age1==4,]
MCHgcn4$r_mass <- MCHgcn4$mass_change - mean45
MCHgcn5 <- MCHgcn[MCHgcn$age1==5,]
MCHgcn5$r_mass <- MCHgcn5$mass_change - mean56
MCHgcn6 <- MCHgcn[MCHgcn$age1==6,]
MCHgcn6$r_mass <- MCHgcn6$mass_change - mean67
MCHgcn7 <- MCHgcn[MCHgcn$age1==7,]
MCHgcn7$r_mass <- MCHgcn7$mass_change - mean78
MCHgcn8 <- MCHgcn[MCHgcn$age1==8,]
MCHgcn8$r_mass <- MCHgcn8$mass_change - mean89
MCHgcn9 <- MCHgcn[MCHgcn$age1==9,]
MCHgcn9$r_mass <- MCHgcn9$mass_change - mean910

MCHgcn <- rbind(MCHgcn4, MCHgcn5, MCHgcn6, MCHgcn7, MCHgcn8, MCHgcn9) 
unique(factor(MCHgcn$id)) #33 obs / 21 ind

# As ID as a random effect creates a singularity, we need to keep only 1 obs / ind
# We keep the observation for which the age is the closest to the mean age for ind with only 1 obs
# mean age:
MCHunique <- MCHgcn %>% group_by(id) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
MCHrepeti <- MCHgcn %>% group_by(id) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(MCHrepeti$id))
mean(MCHunique$age1) # 5.31 yo on average for individuals with 1 obs
MCHrepeti$diff <- abs(MCHrepeti$age1 - 5.31) # calculate difference between age at observation and the mean age of individuals with 1 obs
MCHrepet <-  MCHrepeti %>% 
  group_by(id) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest
MCHrepet <- MCHrepet[,-23] # remove "diff" column so MCHunique and MCHrepet have the same number of column
MCHgcn <- rbind(MCHunique, MCHrepet)
unique(factor(MCHgcn$id)) # 21 individuals


# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
MCHgcn1 <- MCHgcn[!is.na(MCHgcn$FGM1),] # 15
unique(factor(MCHgcn1$id))
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
MCHgcn2 <- MCHgcn[!is.na(MCHgcn$FGM2),] # 18
unique(factor(MCHgcn2$id))
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
MCHgcmean <- MCHgcn[!is.na(MCHgcn$FGM1),] 
MCHgcmean <- MCHgcmean[!is.na(MCHgcmean$FGM2),] # 13
MCHgcmean$mean_fgm <- (MCHgcmean$FGM1 + MCHgcmean$FGM2)/2 # mean FGM during two consecutive years
unique(factor(MCHgcmean$id))

# b) FCHgc: Females, Chiz? ----
FCH <- datagc[datagc$pop == "C" & datagc$sexe == "F",]
FCHgc <- FCH[FCH$ageannee >= 4 & FCH$ageannee <= 10,] # only individuals between 4 and 10 yo
unique(factor(FCHgc$numind)) # 254 obs / 119 ind

# Ordering dataset by individuals and measurement date
FCHgc <- FCHgc[with(FCHgc, order(numind, date)),]
rownames(FCHgc) <- 1:nrow(FCHgc)

# Remove last measurement for each individuals in a first dataset
FCHgc1 <- FCHgc %>% group_by(numind) %>% slice(-n())
FCHgc1 <- FCHgc1[,-c(13,14,15,16,17,19)]

# Remove first measurement for each individuals in a second dataset
FCHgc2 <- FCHgc %>% group_by(numind) %>% slice(-1)
FCHgc2 <- FCHgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
FCHgcn <- cbind(FCHgc1, FCHgc2)
FCHgcn <- FCHgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(FCHgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

# Remove measures with interval >1 year
FCHgcn <- FCHgcn[FCHgcn$age2-FCHgcn$age1==1,]

# Mass change between year 1 and year 2
FCHgcn$mass_change <- FCHgcn$mass2 - FCHgcn$mass1

# Mean mass for all individuals of year1 - Mean mass for all individuals of year2
FCH4 <- FCH[FCH$ageannee==4,]
mean4 <- mean(FCH4$masse)
FCH5 <- FCH[FCH$ageannee==5,]
mean5 <- mean(FCH5$masse)
mean45 <- mean5-mean4
FCH6 <- FCH[FCH$ageannee==6,]
mean6 <- mean(FCH6$masse)
mean56 <- mean6-mean5
FCH7 <- FCH[FCH$ageannee==7,]
mean7 <- mean(FCH7$masse)
mean67 <- mean7-mean6
FCH8 <- FCH[FCH$ageannee==8,]
mean8 <- mean(FCH8$masse)
mean78 <- mean8-mean7
FCH9 <- FCH[FCH$ageannee==9,]
mean9 <- mean(FCH9$masse)
mean89 <- mean9-mean8
FCH10 <- FCH[FCH$ageannee==10,]
mean10 <- mean(FCH10$masse)
mean910 <- mean10-mean9

FCHgcn4 <- FCHgcn[FCHgcn$age1==4,]
FCHgcn4$r_mass <- FCHgcn4$mass_change - mean45
FCHgcn5 <- FCHgcn[FCHgcn$age1==5,]
FCHgcn5$r_mass <- FCHgcn5$mass_change - mean56
FCHgcn6 <- FCHgcn[FCHgcn$age1==6,]
FCHgcn6$r_mass <- FCHgcn6$mass_change - mean67
FCHgcn7 <- FCHgcn[FCHgcn$age1==7,]
FCHgcn7$r_mass <- FCHgcn7$mass_change - mean78
FCHgcn8 <- FCHgcn[FCHgcn$age1==8,]
FCHgcn8$r_mass <- FCHgcn8$mass_change - mean89
FCHgcn9 <- FCHgcn[FCHgcn$age1==9,]
FCHgcn9$r_mass <- FCHgcn9$mass_change - mean910

FCHgcn <- rbind(FCHgcn4, FCHgcn5, FCHgcn6, FCHgcn7, FCHgcn8, FCHgcn9) # 77 obs / 43 ind
unique(factor(FCHgcn$id))

# As ID as a random effect creates a singularity, we need to keep only 1 obs / ind
# We keep the observation for which the age is the closest to the mean age for ind with only 1 obs
# mean age:
FCHunique <- FCHgcn %>% group_by(id) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
FCHrepeti <- FCHgcn %>% group_by(id) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(FCHrepeti$id))
mean(FCHunique$age1) # 5.56 yo on average for individuals with 1 obs
FCHrepeti$diff <- abs(FCHrepeti$age1 - 5.56) # calculate difference between age at observation and the mean age of individuals with 1 obs
FCHrepet <-  FCHrepeti %>% 
  group_by(id) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest
FCHrepet <- FCHrepet[,-23] # remove "diff" column so FCHunique and FCHrepet have the same number of column
FCHgcn <- rbind(FCHunique, FCHrepet)
unique(factor(FCHgcn$id)) # 43 individuals

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
FCHgcn1 <- FCHgcn[!is.na(FCHgcn$FGM1),] # 25
unique(factor(FCHgcn1$id))
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
FCHgcn2 <- FCHgcn[!is.na(FCHgcn$FGM2),] # 31
unique(factor(FCHgcn2$id))
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
FCHgcmean <- FCHgcn[!is.na(FCHgcn$FGM1),]
FCHgcmean <- FCHgcmean[!is.na(FCHgcmean$FGM2),] # 19
FCHgcmean$mean_fgm <- (FCHgcmean$FGM1 + FCHgcmean$FGM2)/2 # mean FGM during two consecutive years
unique(factor(FCHgcmean$id))

# c) MTFgc: Males, Trois-Fontaine ----
MTF <- datagc[datagc$pop == "3F" & datagc$sexe == "M",]
MTFgc <- MTF[MTF$ageannee >= 4 & MTF$ageannee <= 10,] # only individuals between 4 and 10 yo
unique(factor(MTFgc$numind)) # 166 obs / 77 ind

# Ordering dataset by individuals and measurement date
MTFgc <- MTFgc[with(MTFgc, order(numind, date)),]
rownames(MTFgc) <- 1:nrow(MTFgc)

# Remove last measurement for each individuals in a first dataset
MTFgc1 <- MTFgc %>% group_by(numind) %>% slice(-n())
MTFgc1 <- MTFgc1[,-c(13,14,15,16,17,19)]
unique(factor(MTFgc1$numind))

# Remove first measurement for each individuals in a second dataset
MTFgc2 <- MTFgc %>% group_by(numind) %>% slice(-1)
unique(factor(MTFgc2$numind))
MTFgc2 <- MTFgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
MTFgcn <- cbind(MTFgc1, MTFgc2)
MTFgcn <- MTFgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(MTFgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

# Remove measures with interval >1 year
MTFgcn <- MTFgcn[MTFgcn$age2-MTFgcn$age1==1,]

# Mass change between year 1 and year 2
MTFgcn$mass_change <- MTFgcn$mass2 - MTFgcn$mass1

# Mean mass for all individuals of year1 - Mean mass for all individuals of year2
MTF4 <- MTF[MTF$ageannee==4,]
mean4 <- mean(MTF4$masse)
MTF5 <- MTF[MTF$ageannee==5,]
mean5 <- mean(MTF5$masse)
mean45 <- mean5-mean4
MTF6 <- MTF[MTF$ageannee==6,]
mean6 <- mean(MTF6$masse)
mean56 <- mean6-mean5
MTF7 <- MTF[MTF$ageannee==7,]
mean7 <- mean(MTF7$masse)
mean67 <- mean7-mean6
MTF8 <- MTF[MTF$ageannee==8,]
mean8 <- mean(MTF8$masse)
mean78 <- mean8-mean7
MTF9 <- MTF[MTF$ageannee==9,]
mean9 <- mean(MTF9$masse)
mean89 <- mean9-mean8
MTF10 <- MTF[MTF$ageannee==10,]
mean10 <- mean(MTF10$masse)
mean910 <- mean10-mean9

MTFgcn4 <- MTFgcn[MTFgcn$age1==4,]
MTFgcn4$r_mass <- MTFgcn4$mass_change - mean45
MTFgcn5 <- MTFgcn[MTFgcn$age1==5,]
MTFgcn5$r_mass <- MTFgcn5$mass_change - mean56
MTFgcn6 <- MTFgcn[MTFgcn$age1==6,]
MTFgcn6$r_mass <- MTFgcn6$mass_change - mean67
MTFgcn7 <- MTFgcn[MTFgcn$age1==7,]
MTFgcn7$r_mass <- MTFgcn7$mass_change - mean78
MTFgcn8 <- MTFgcn[MTFgcn$age1==8,]
MTFgcn8$r_mass <- MTFgcn8$mass_change - mean89
MTFgcn9 <- MTFgcn[MTFgcn$age1==9,]
MTFgcn9$r_mass <- MTFgcn9$mass_change - mean910

MTFgcn <- rbind(MTFgcn4, MTFgcn5, MTFgcn6, MTFgcn7, MTFgcn8, MTFgcn9) 
unique(factor(MTFgcn$id)) # 67 obs / 37 ind

# As ID as a random effect creates a singularity, we need to keep only 1 obs / ind
# We keep the observation for which the age is the closest to the mean age for ind with only 1 obs
# mean age:
MTFunique <- MTFgcn %>% group_by(id) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
MTFrepeti <- MTFgcn %>% group_by(id) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(MTFrepeti$id))
mean(MTFunique$age1) # 5.67 yo on average for individuals with 1 obs
MTFrepeti$diff <- abs(MTFrepeti$age1 - 5.67) # calculate difference between age at observation and the mean age of individuals with 1 obs
MTFrepet <-  MTFrepeti %>% 
  group_by(id) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest
MTFrepet <- MTFrepet[,-23] # remove "diff" column so MCHunique and MCHrepet have the same number of column
MTFgcn <- rbind(MTFunique, MTFrepet)
unique(factor(MTFgcn$id)) # 37 individuals

# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
MTFgcn1 <- MTFgcn[!is.na(MTFgcn$FGM1),] # 17
unique(factor(MTFgcn1$id))
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
MTFgcn2 <- MTFgcn[!is.na(MTFgcn$FGM2),] # 18
unique(factor(MTFgcn2$id))
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
MTFgcmean <- MTFgcn[!is.na(MTFgcn$FGM1),]
MTFgcmean <- MTFgcmean[!is.na(MTFgcmean$FGM2),] # 9
MTFgcmean$mean_fgm <- (MTFgcmean$FGM1 + MTFgcmean$FGM2)/2 # mean FGM during two consecutive years
unique(factor(MTFgcmean$id))

# d) FTFgc: Females, Trois-Fontaine ----
FTF <- datagc[datagc$pop == "3F" & datagc$sexe == "F",]
FTFgc <- FTF[FTF$ageannee >= 4 & FTF$ageannee <= 10,] # only individuals between 4 and 10 yo
unique(factor(FTFgc$numind)) # 182 obs / 94 ind

# Ordering dataset by individuals and measurement date
FTFgc <- FTFgc[with(FTFgc, order(numind, date)),]
rownames(FTFgc) <- 1:nrow(FTFgc)

# Remove last measurement for each individuals in a first dataset
FTFgc1 <- FTFgc %>% group_by(numind) %>% slice(-n())
FTFgc1 <- FTFgc1[,-c(13,14,15,16,17,19)]
unique(factor(FTFgc1$numind))

# Remove first measurement for each individuals in a second dataset
FTFgc2 <- FTFgc %>% group_by(numind) %>% slice(-1)
unique(factor(FTFgc2$numind))
FTFgc2 <- FTFgc2[,c(1,7,8,10,11,12,18)]

# Bind the two datasets to have one measure and the next one on the same row
FTFgcn <- cbind(FTFgc1, FTFgc2)
FTFgcn <- FTFgcn[,c(2,1,3,6,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(FTFgcn)[] <- c("id", "idkit1", "pop", "sex", "cohorte", "quali_coh", "annee1", "date1", "julian1", "qualite1", "age1", "mass1", "FGM1", "idkit2", "annee2", "date2", "qualite2", "age2", "mass2","FGM2")

# Remove measures with interval >1 year
FTFgcn <- FTFgcn[FTFgcn$age2-FTFgcn$age1==1,]

# Mass change between year 1 and year 2
FTFgcn$mass_change <- FTFgcn$mass2 - FTFgcn$mass1

# Mean mass for all individuals of year1 - Mean mass for all individuals of year2
FTF4 <- FTF[FTF$ageannee==4,]
mean4 <- mean(FTF4$masse)
FTF5 <- FTF[FTF$ageannee==5,]
mean5 <- mean(FTF5$masse)
mean45 <- mean5-mean4
FTF6 <- FTF[FTF$ageannee==6,]
mean6 <- mean(FTF6$masse)
mean56 <- mean6-mean5
FTF7 <- FTF[FTF$ageannee==7,]
mean7 <- mean(FTF7$masse)
mean67 <- mean7-mean6
FTF8 <- FTF[FTF$ageannee==8,]
mean8 <- mean(FTF8$masse)
mean78 <- mean8-mean7
FTF9 <- FTF[FTF$ageannee==9,]
mean9 <- mean(FTF9$masse)
mean89 <- mean9-mean8
FTF10 <- FTF[FTF$ageannee==10,]
mean10 <- mean(FTF10$masse)
mean910 <- mean10-mean9


FTFgcn4 <- FTFgcn[FTFgcn$age1==4,]
FTFgcn4$r_mass <- FTFgcn4$mass_change - mean45
FTFgcn5 <- FTFgcn[FTFgcn$age1==5,]
FTFgcn5$r_mass <- FTFgcn5$mass_change - mean56
FTFgcn6 <- FTFgcn[FTFgcn$age1==6,]
FTFgcn6$r_mass <- FTFgcn6$mass_change - mean67
FTFgcn7 <- FTFgcn[FTFgcn$age1==7,]
FTFgcn7$r_mass <- FTFgcn7$mass_change - mean78
FTFgcn8 <- FTFgcn[FTFgcn$age1==8,]
FTFgcn8$r_mass <- FTFgcn8$mass_change - mean89
FTFgcn9 <- FTFgcn[FTFgcn$age1==9,]
FTFgcn9$r_mass <- FTFgcn9$mass_change - mean910

FTFgcn <- rbind(FTFgcn4, FTFgcn5, FTFgcn6, FTFgcn7, FTFgcn8, FTFgcn9) 
unique(factor(FTFgcn$id)) # 61 obs / 31 ind

# As ID as a random effect creates a singularity, we need to keep only 1 obs / ind
# We keep the observation for which the age is the closest to the mean age for ind with only 1 obs
# mean age:
FTFunique <- FTFgcn %>% group_by(id) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
FTFrepeti <- FTFgcn %>% group_by(id) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(FTFrepeti$id))
mean(FTFunique$age1) # 5.71 yo on average for individuals with 1 obs
FTFrepeti$diff <- abs(FTFrepeti$age1 - 5.71) # calculate difference between age at observation and the mean age of individuals with 1 obs
FTFrepet <-  FTFrepeti %>% 
  group_by(id) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest
FTFrepet <- FTFrepet[,-23] # remove "diff" column so FTFunique and FTFrepet have the same number of column
FTFgcn <- rbind(FTFunique, FTFrepet)
unique(factor(FTFgcn$id)) # 31 individuals


# To test for the effect of the FGM in year 1, we are only limited by the FGM on year 1
FTFgcn1 <- FTFgcn[!is.na(FTFgcn$FGM1),] # 15
unique(factor(FTFgcn1$id))
# To test for the effect of the FGM in year 2, we are only limited by the FGM on year 2
FTFgcn2 <- FTFgcn[!is.na(FTFgcn$FGM2),] # 17
unique(factor(FTFgcn2$id))
# To test for the effect of the mean FGM between the 2 years, we are limited by the FGM on both year
FTFgcmean <- FTFgcn[!is.na(FTFgcn$FGM1),]
FTFgcmean <- FTFgcmean[!is.na(FTFgcmean$FGM2),] # 11
FTFgcmean$mean_fgm <- (FTFgcmean$FGM1 + FTFgcmean$FGM2)/2 # mean FGM during two consecutive years
unique(factor(FTFgcmean$id))

# MODELS ----
# Pool the 4 subsets:
datagcn1 <- rbind(FCHgcn1, MCHgcn1, FTFgcn1, MTFgcn1)
unique(factor(datagcn1$id)) # 131 obs / 90 ind
datagcn2 <- rbind(FCHgcn2, MCHgcn2, FTFgcn2, MTFgcn2)
unique(factor(datagcn2$id)) # 147 obs / 101 ind
datagcmean <- rbind(FCHgcmean, MCHgcmean, FTFgcmean, MTFgcmean)
unique(factor(datagcmean$id)) # 91 obs / 71 ind

# FIGURES ----
# Fig2B - Prime age adults
datagcn1$logFGM1 <- log(datagcn1$FGM1) 
modb1 <- lmer(r_mass ~ scale(logFGM1, scale=F) + (1|cohorte), data = datagcn1)
summary(modb1)

newdatab1=expand.grid(logFGM1=seq(min(log(datagcn1$FGM1)),
                                  max(log(datagcn1$FGM1)), 0.005))

predb1=predictSE(modb1,newdata=newdatab1,se.fit=TRUE)
newdatab1=cbind(newdatab1,predb1)
newdatab1$low=newdatab1$fit-1.44*newdatab1$se.fit
newdatab1$upp=newdatab1$fit+1.44*newdatab1$se.fit

datagcmean$logmeanfgm <- log(datagcmean$mean_fgm) 
modb2 <- lmer(r_mass ~ scale(logmeanfgm, scale=F) + (1|cohorte), data = datagcmean)
summary(modb2)

newdatab2=expand.grid(logmeanfgm=seq(min(log(datagcmean$mean_fgm)),
                                     max(log(datagcmean$mean_fgm)), 0.005))

predb2=predictSE(modb2,newdata=newdatab2,se.fit=TRUE)
newdatab2=cbind(newdatab2,predb2)
newdatab2$low=newdatab2$fit-1.44*newdatab2$se.fit
newdatab2$upp=newdatab2$fit+1.44*newdatab2$se.fit

grobb <- grobTree(textGrob("Prime age adults", x=0.05,  y=0.97, hjust=0,
                           gp=gpar(col="black", fontsize=14, face="bold")))
fig2b <- ggplot(data = datagcn1, aes(y=r_mass, x=logFGM1)) + 
  geom_point(data = datagcn1, aes(y=r_mass, x=logFGM1, color="#5FB55F")) +
  geom_point(data = datagcmean, aes(y=r_mass, x=logmeanfgm, color="#D35FB7")) +
  scale_y_continuous(name = "Change in relative mass (kg)", labels = label_number(accuracy = 0.1)) +
  xlab("FGM levels (log-transformed)") +
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
        legend.background = element_rect(fill='transparent'))+
  scale_color_manual(values=c("#D35FB7", "#5FB55F"), 
                     name="FGM measurement",
                     labels=c("1st year", "Mean 1st and 2nd year"))
fig2b

f2 <- arrangeGrob(fig2a, fig2b, ncol=2, nrow=1)
grid.newpage()
grid.draw(f2)

# Save this figure
ggsave(filename = "results/figures/fig2.png", 
       f2, width = 30, height = 15, dpi = 300, units = "cm", device='png')
