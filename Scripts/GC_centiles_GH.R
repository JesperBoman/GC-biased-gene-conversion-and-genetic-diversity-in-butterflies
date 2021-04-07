# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# Main script for linear models and making plots for the across the genome part of 
# empirical inferences of gBGC in Leptidea butterflies. 
# ==========================================================================================================
# Jesper Boman                      5 mar 2021                            Boman et al. 2021
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### NOTE ####
#Here we used windows from a gradient of GC content across the Leptidea genome 
#but you could adapt the code presented here to other cases were DAFs are created on the basis of e.g. 
#genomic windows or chromosomes.

##### SETUP ENVIRONMENT ####

install.packages("wesanderson")
library(wesanderson)
col <- wes_palette("Darjeeling1")
col2 <-wes_palette("Darjeeling2")

install.packages("ggplot2")
library(ggplot2)

library(car)

#### INPUT DATA ####

#Here we load two files per population. One from the main output of the script making our GC centiles (i.e. contains the DAF spectra), e.g. swe_sin_GC_centiles.
# The second files is obtained from the output of the inference of gBGC using the "gBGC_Estimation_GC_centiles_GH.nb" in Mathematica, e.g. swe_sin_GC_centiles.results.
swe_sin_GC_centiles <- read.table(file = file.choose(),header = T)
swe_sin_GC_centiles$population <- "Swe-sin"
swe_sin_GC_centiles.results <- read.table(file = file.choose(),header = T)
swe_sin_GC_centiles.results$population <- "Swe-sin"

spa_sin_GC_centiles <- read.table(file = file.choose(),header = T)
spa_sin_GC_centiles$population <- "Spa-sin"
spa_sin_GC_centiles.results <- read.table(file = file.choose(),header = T)
spa_sin_GC_centiles.results$population <- "Spa-sin"

kaz_sin_GC_centiles <- read.table(file = file.choose(),header = T)
kaz_sin_GC_centiles$population <- "Kaz-sin"
kaz_sin_GC_centiles.results <- read.table(file = file.choose(),header = T)
kaz_sin_GC_centiles.results$population <- "Kaz-sin"

kaz_juv_GC_centiles <- read.table(file = file.choose(),header = T)
kaz_juv_GC_centiles$population <- "Kaz-juv"
kaz_juv_GC_centiles.results <- read.table(file = file.choose(),header = T)
kaz_juv_GC_centiles.results$population <- "Kaz-juv"

ire_juv_GC_centiles <- read.table(file = file.choose(),header = T)
ire_juv_GC_centiles$population <- "Ire-juv"
ire_juv_GC_centiles.results <- read.table(file = file.choose(),header = T)
ire_juv_GC_centiles.results$population <- "Ire-juv"

spa_rea_GC_centiles <- read.table(file = file.choose(),header = T)
spa_rea_GC_centiles$population <- "Spa-rea"
spa_rea_GC_centiles.results <- read.table(file = file.choose(),header = T)
spa_rea_GC_centiles.results$population <- "Spa-rea"

#!#
GC_centiles <- rbind(swe_sin_GC_centiles, spa_sin_GC_centiles, kaz_sin_GC_centiles, kaz_juv_GC_centiles, ire_juv_GC_centiles, spa_rea_GC_centiles)
GC_centiles.results <- rbind(swe_sin_GC_centiles.results, spa_sin_GC_centiles.results, kaz_sin_GC_centiles.results, kaz_juv_GC_centiles.results, ire_juv_GC_centiles.results, spa_rea_GC_centiles.results)
#!#









##### Centiles ####


## Likelihood-ratio test between model M1 and M0 ##
pchisq(-2*(mean(GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]$lnL0)-mean(GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]$lnL1)), df=1, lower.tail=F)

## Likelihood-ratio test between model M1 and M1corr (aka M1*) ##
pchisq(-2*(mean(GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]$lnL1)-mean(GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]$lnL1cor)), df=3, lower.tail=F)
##
 



ggplot(GC_centiles.results, aes(x=localGC, y=total_pi, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  #stat_smooth(method ="lm", formula=y~poly(x,3))+
  ylab(expression(pi)) +
  ylim(0,0.008)+
  xlab("Observed GC content") +
  labs(title="Intergenic sites") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


#### Mean DAF, a measure that can be used as a proxy for the strength of gBGC. Not better than inferring using the Glémin et al. 2015 method.
GC_centiles$segregating_sites_anchor_per_class <- GC_centiles$X1+
  GC_centiles$X2+GC_centiles$X3+GC_centiles$X4+GC_centiles$X5+GC_centiles$X6+GC_centiles$X7+GC_centiles$X8+GC_centiles$X9+GC_centiles$X10+
  GC_centiles$X11+GC_centiles$X12+GC_centiles$X13+GC_centiles$X14+GC_centiles$X15+GC_centiles$X16+GC_centiles$X17+GC_centiles$X18+GC_centiles$X19

GC_centiles$mean_daf_per_class <- NA

for (i in 1:length(GC_centiles$class)){
  GC_centiles$mean_daf_per_class[i] <- sum(GC_centiles$X1[i]*0.05,
                                         GC_centiles$X2[i]*0.1,GC_centiles$X3[i]*0.15,GC_centiles$X4[i]*0.2,GC_centiles$X5[i]*0.25,GC_centiles$X6[i]*0.3,GC_centiles$X7[i]*0.35,GC_centiles$X8[i]*0.4,GC_centiles$X9[i]*0.45,GC_centiles$X10[i]*0.5,
                                         GC_centiles$X11[i]*0.55,GC_centiles$X12[i]*0.6,GC_centiles$X13[i]*0.65,GC_centiles$X14[i]*0.7,GC_centiles$X15[i]*0.75,GC_centiles$X16[i]*0.8,GC_centiles$X17[i]*0.85,GC_centiles$X18[i]*0.9,GC_centiles$X19[i]*0.95)/GC_centiles$segregating_sites_anchor_per_class[i]
}

GC_centiles.results$meanDAF_N <- GC_centiles[GC_centiles$class == "Neutre",]$mean_daf_per_class
GC_centiles.results$meanDAF_SW <- GC_centiles[GC_centiles$class == "SW",]$mean_daf_per_class
GC_centiles.results$meanDAF_WS <- GC_centiles[GC_centiles$class == "WS",]$mean_daf_per_class
GC_centiles.results$meanDAF_WS_m_SW <- GC_centiles[GC_centiles$class == "WS",]$mean_daf_per_class - GC_centiles[GC_centiles$class == "SW",]$mean_daf_per_class


#Linear models

par(mfrow=c(2,2))
summary(lm(meanDAF_WS_m_SW~localGC+population, data=GC_centiles.results))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))
plot(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(B_M1~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(B_M1~total_pi, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
plot(lm(B_M1~total_pi, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))


par(mfrow=c(6,4))
plot(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
plot(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
plot(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
plot(lm(log(lambda_M1)~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
plot(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
plot(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
plot(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))

##FIGURE 3 plots  ####

"Mean DAF difference between W-->S and S-->W. 6.04 x 5"
ggplot(GC_centiles.results, aes(x=localGC, y=meanDAF_WS_m_SW, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  geom_hline(yintercept = 0)+
  geom_abline(slope=coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[1], col = col[1]) +
  geom_abline(slope=coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[1], col = col[3]) +
  geom_abline(slope=coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[1], col = col[4]) +
  geom_abline(slope=coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[2], intercept = coef(lm(meanDAF_WS_m_SW~localGC, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[1], col = col[5]) +
  
  #stat_smooth(method ="lm", formula=y~poly(x,1))+
  ylab("Mean DAF WS - Mean DAF SW")+
  xlab("Observed GC content") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


"B v GC L. sinapis and L. reali"
ggplot(GC_centiles.results[c(grep("sin", GC_centiles.results$population),grep("rea", GC_centiles.results$population)),], aes(x=B_M1, y=localGC, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[4],col2[2], col[1], col[3])) + 
  theme_classic() +
  geom_abline(slope=coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[1], col = col[1]) +
  geom_abline(slope=coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[1], col = col[3]) +
  geom_abline(slope=coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[1], col = col[4]) +

  #stat_smooth(method ="lm", formula=y~x, alpha=0.3)+
  xlab("B_M1") +
  ylab("Observed GC content") +
  #labs(title="Intergenic sites") +
  theme(legend.position=c(0.101, 0.841), legend.background =element_rect(fill=NA, size=0.5, linetype="solid", colour ="black"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


"B v GC L. juvernica"
ggplot(GC_centiles.results[grep("juv", GC_centiles.results$population),], aes(y=localGC, x=B_M1, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5])) + 
  theme_classic() +
  geom_abline(slope=coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[2], intercept = coef(lm(localGC~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[1], col = col[5]) +
  
  #stat_smooth(method ="lm", formula=y~x, alpha=0.3)+
  xlab("B_M1") +
  ylab("Observed GC content") +
  #labs(title="Intergenic sites") +
  theme(legend.position=c(0.101, 0.831), legend.background =element_rect(fill=NA, size=0.5, linetype="solid", colour ="black"), panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


"B vs pi, 6.04x 5"
ggplot(GC_centiles.results, aes(x=total_pi, y=B_M1, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  #stat_smooth(method ="lm", formula=y~poly(x,1))+
  ylab("B_M1")+
  geom_abline(slope=coef(lm(B_M1~total_pi, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(B_M1~total_pi, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[1], col = col[3]) +
  geom_abline(slope=coef(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(B_M1~total_pi, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[1], col = col[4]) +
  
  xlab(expression(pi)) +
  theme(legend.position = c(0.9,0.8), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))



"B vs lambda, 6.04 x 5"
ggplot(GC_centiles.results, aes(y=lambda_M1, x=B_M1, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  geom_abline(slope=coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[1], col = col[4]) +
  geom_abline(slope=coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))[1], col = col[5]) +
  geom_abline(slope=coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))[1], col = col[2]) +
  geom_abline(slope=coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~B_M1, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))[1], col = col2[2]) +
  
  #stat_smooth(method ="lm", formula=y~poly(x,1))+
  xlab("B_M1")+
  ylab(expression(lambda)) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


"GC content vs lambda, L. sinapis and L. reali 6.04 x 5"
ggplot(GC_centiles.results[c(grep("sin", GC_centiles.results$population),grep("rea", GC_centiles.results$population)),], aes(y=lambda_M1, x=localGC, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[4],  col2[2], col[1], col[3])) + 
  geom_abline(slope=coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))[1], col = col[3]) +
  geom_abline(slope=coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))[1], col = col[1]) +
  geom_abline(slope=coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))[1], col = col[4]) +
  geom_abline(slope=coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))[2], intercept = coef(lm(lambda_M1~localGC, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))[1], col = col2[2]) +
  theme_classic() + 
  #stat_smooth(method ="lm", formula=y~poly(x,1))+
  xlab("Observed GC content")+
  ylab(expression(lambda)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

"GC content vs lambda, L. juvernica 6.04 x 5"
ggplot(GC_centiles.results[grep("juv", GC_centiles.results$population),], aes(y=lambda_M1, x=localGC, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5])) + 
  theme_classic() + 
  #stat_smooth(method ="lm", formula=y~poly(x,1))+
  xlab("Observed GC content")+
  ylab(expression(lambda)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))






#####Table 1####
aggregate(eWS~population, GC_centiles.results, mean)
aggregate(eSW~population, GC_centiles.results, mean)

aggregate(localGC~population, GC_centiles.results, mean)
aggregate(lambda_M1~population, GC_centiles.results, mean)
aggregate(B_M1~population, GC_centiles.results, mean)

#Standard errors of B and lambda
cbind(aggregate(B_M1~population, GC_centiles.results, sd), aggregate(B_M1~population, GC_centiles.results, sd)[,2]/10)
cbind(aggregate(lambda_M1~population, GC_centiles.results, sd), aggregate(lambda_M1~population, GC_centiles.results, sd)[,2]/10)

aggregate(1/(1+lambda_M1)~population, GC_centiles.results, mean)
aggregate(1/(1+lambda_M1*exp(-B_M1))~population, GC_centiles.results, mean)

#GC π-max
pi_temp <- aggregate(total_pi~population, GC_centiles.results, max)
for (i in unique(GC_centiles.results$population)){
  gc_temp <- GC_centiles.results[GC_centiles.results$population == i & GC_centiles.results$total_pi == pi_temp[pi_temp$population == i,]$total_pi,]$localGC
  cat(i, gc_temp, "\n")
}
aggregate(total_pi~population, GC_centiles.results, max)

#GC CDS-min
CDS_temp <- aggregate(cds_total_density~population, GC_centiles.results, min)
for (i in unique(GC_centiles.results$population)){
  gc_temp <- GC_centiles.results[GC_centiles.results$population == i & GC_centiles.results$cds_total_density == CDS_temp[CDS_temp$population == i,]$cds_total_density,]$localGC
  cat(i, gc_temp, "\n")
}


#### Figure 4 prep ####
GC_centiles.results$cds_total_density <- GC_centiles.results$cds_total/(GC_centiles.results$invar_L+GC_centiles.results$anchor_L+GC_centiles.results$unanchor_L)
GC_centiles.results$LiBulmerGCeq <- 1/(1+GC_centiles.results$lambda_M1*exp(-GC_centiles.results$B_M1))
GC_centiles.results$gcDist_from_eq <- GC_centiles.results$localGC-GC_centiles.results$LiBulmerGCeq

GC_centiles.results$predicted_pi <- ((1/(1+GC_centiles.results$lambda_M1*exp(-GC_centiles.results$B_M1)))*GC_centiles.results$lambda_M1*2*(1/(1-exp(GC_centiles.results$B_M1))+1/(GC_centiles.results$B_M1)) +
   (1-((1/(1+GC_centiles.results$lambda_M1*exp(-GC_centiles.results$B_M1)))))*2*(1/(1-exp(-GC_centiles.results$B_M1))-1/(GC_centiles.results$B_M1)) + GC_centiles.results$Theta/GC_centiles.results$ThetaWS) / 
  (((2*GC_centiles.results$lambda_M1)/(GC_centiles.results$lambda_M1+1) + GC_centiles.results$Theta/GC_centiles.results$ThetaWS))

GC_centiles.results$predicted_pi
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]))
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]))
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]))
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]))
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]))
summary(lm(predicted_pi~total_pi, data=GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]))


summary(GC_centiles.results[GC_centiles.results$population == "Swe-sin",]$predicted_pi)
summary(GC_centiles.results[GC_centiles.results$population == "Spa-sin",]$predicted_pi)
summary(GC_centiles.results[GC_centiles.results$population == "Kaz-sin",]$predicted_pi)
summary(GC_centiles.results[GC_centiles.results$population == "Kaz-juv",]$predicted_pi)
summary(GC_centiles.results[GC_centiles.results$population == "Ire-juv",]$predicted_pi)
summary(GC_centiles.results[GC_centiles.results$population == "Spa-rea",]$predicted_pi)






#### FIGURE 4 PLOTS ####
GC_centiles.results <- GC_centiles.results[order(GC_centiles.results$population),]


# GC equilibrium plots

#Mutational equilibrium
ggplot(GC_centiles.results, aes(y=1/(1+lambda_M1), x=localGC, col = population )) + geom_point() + 
  scale_color_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  ylab("Equilibrium GC content") +
  xlab("Observed GC content") +
  xlim(0.15,0.55)+
  ylim(0.15,0.55)+
  geom_abline(intercept=0, slope=1, lty=3)+
  theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

#gBGC-mutation equilibrium
ggplot(GC_centiles.results, aes(y=1/(1+lambda_M1*exp(-B_M1)), x=localGC, col = population )) + geom_point() + 
  scale_color_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  ylab("Equilibrium GC content") +
  xlab("Observed GC content") +
  xlim(0.15,0.55)+
  ylim(0.15,0.55)+
  geom_abline(intercept=0, slope=1, lty=3)+
  theme(legend.position="none" ,panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


#Skewness
ggplot(GC_centiles.results, aes(x=localGC, y=skewness, col = population )) + geom_point() + 
  scale_color_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3]), labels = c("Ire-juv:  0.22", "Kaz-juv: 0.19", "Kaz-sin: 0.25", "Spa-rea: 0.25", "Spa-sin:0.25 & 0.30", "Swe-sin: 0.25 & 0.29")) + 
  theme_classic() + 
  ylab("Skewness") +
  xlab("Observed GC content") +
  #stat_smooth(method ="lm", formula=y~poly(x,2))+
  geom_hline(yintercept = 0, lty=1)+
  geom_vline(col = col[3], xintercept = 1/(1+mean(GC_centiles.results[grep("Swe-sin", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  geom_vline(col = col[1], xintercept = 1/(1+mean(GC_centiles.results[grep("Spa-sin", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  geom_vline(col = col[4], xintercept = 1/(1+mean(GC_centiles.results[grep("Kaz-sin", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  geom_vline(col = col[5], xintercept = 1/(1+mean(GC_centiles.results[grep("Kaz-juv", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  geom_vline(col = col[2], xintercept = 1/(1+mean(GC_centiles.results[grep("Ire-juv", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  geom_vline(col = col2[2], xintercept = 1/(1+mean(GC_centiles.results[grep("Spa-rea", GC_centiles.results$population),]$lambda_M1)), lty=3) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position="none",legend.background = element_rect(fill = "lightgray"), legend.title = element_text(size=14), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12))


#GC vs π
ggplot(GC_centiles.results, aes(x=localGC, y=total_pi, col = population)) + geom_point() + 
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  xlim(0.2,0.5)+
  ylab(expression(pi))+
  xlab("Observed GC content") +
  theme( legend.box.background = element_rect(colour = "black", size = 0.5), legend.position = c(0.895,0.786),panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))











#### ALLMUT - DAF spectra for separate mutation categories ####

#NOTE
#To create the input file here I separated mutations of the type T-->A and A-->T in the data. 
#However, we do not actually know just from simple resequencing data what strand the mutation occured on so T-->A and A-->T etc. need to be grouped.
#

allMut <- read.table(file = file.choose(), header=T)
allMut <- allMut[order(allMut$population),]
allMut$segregating_sites_anchor_per_class <- allMut$X1+
  allMut$X2+allMut$X3+allMut$X4+allMut$X5+allMut$X6+allMut$X7+allMut$X8+allMut$X9+allMut$X10+
  allMut$X11+allMut$X12+allMut$X13+allMut$X14+allMut$X15+allMut$X16+allMut$X17+allMut$X18+allMut$X19

#Labeling of different mutation categories
for (i in 1:length(allMut$class)){
  if(allMut$class[i] == "A_T" | allMut$class[i] == "T_A"){
    allMut$class4cat[i] <- "WW"
    allMut$classTT[i] <- "Transversion"
  }
  if(allMut$class[i] == "C_G" | allMut$class[i] == "G_C"){
    allMut$class4cat[i] <- "SS"
    allMut$classTT[i] <- "Transversion"
  }
  if(allMut$class[i] == "C_A" | allMut$class[i] == "C_T" | allMut$class[i] == "G_A" | allMut$class[i] == "G_T"){
    allMut$class4cat[i] <- "SW"
    if(allMut$class[i] == "C_A" | allMut$class[i] == "G_T"){
      allMut$classTT[i] <- "Transversion"
    }
    if(allMut$class[i] == "C_T" | allMut$class[i] == "G_A"){
      allMut$classTT[i] <- "Transition"
    }
    
  }
  if(allMut$class[i] == "A_C" | allMut$class[i] == "A_G" | allMut$class[i] == "T_C" | allMut$class[i] == "T_G"){
    allMut$class4cat[i] <- "WS"
    if(allMut$class[i] == "A_C" | allMut$class[i] == "T_G"){
      allMut$classTT[i] <- "Transversion"
    }
    if(allMut$class[i] == "A_G" | allMut$class[i] == "T_C"){
      allMut$classTT[i] <- "Transition"
    }
  }
}

#Calculating the summed pairwise differences per mutation category
allMut$DAF_pi <- (allMut$X1*19)/190 + 
  (allMut$X2*18*2)/190 + 
  (allMut$X3*17*3)/190 + 
  (allMut$X4*16*4)/190 + 
  (allMut$X5*15*5)/190 + 
  (allMut$X6*14*6)/190 + 
  (allMut$X7*13*7)/190 +
  (allMut$X8*12*8)/190 +
  (allMut$X9*11*9)/190 +
  (allMut$X10*10*10)/190 +
  (allMut$X11*9*11)/190 +
  (allMut$X12*8*12)/190 +
  (allMut$X13*7*13)/190 +
  (allMut$X14*6*14)/190 +
  (allMut$X15*5*15)/190 +
  (allMut$X16*4*16)/190 +
  (allMut$X17*3*17)/190 +
  (allMut$X18*2*18)/190 +
  (allMut$X19*19)/190

#Aggregating the summed pairwise differences per class4cat (i.e. S-->S, W-->W, W-->S and S-->W)
daf_pi_class4 <- aggregate(DAF_pi~centile+class4cat+population, allMut, sum)


#Here we calculate the demography-corrected π-values (not used in Boman et al. 2021)
allMut$X2r <-NA
allMut$X3r <-NA
allMut$X4r <-NA
allMut$X5r <-NA
allMut$X6r <-NA
allMut$X7r <-NA
allMut$X8r <-NA
allMut$X9r <-NA
allMut$X10r <-NA
allMut$X11r <-NA
allMut$X12r <-NA
allMut$X13r <-NA
allMut$X14r <-NA
allMut$X15r <-NA
allMut$X16r <-NA
allMut$X17r <-NA
allMut$X18r <-NA
allMut$X19r <-NA

for (i in 1:length(GC_centiles.results$r2)){
  matches <- (GC_centiles.results$region[i] == allMut$centile & GC_centiles.results$population[i] == allMut$population)
  m_list <- (1:length(matches))[matches]
  allMut$X2r[m_list] <-allMut$X2[m_list] * GC_centiles.results$r2[i]
  allMut$X3r[m_list] <-allMut$X3[m_list] * GC_centiles.results$r3[i]
  allMut$X4r[m_list] <-allMut$X4[m_list] * GC_centiles.results$r4[i]
  allMut$X5r[m_list] <-allMut$X5[m_list] * GC_centiles.results$r5[i]
  allMut$X6r[m_list] <-allMut$X6[m_list] * GC_centiles.results$r6[i]
  allMut$X7r[m_list] <-allMut$X7[m_list] * GC_centiles.results$r7[i]
  allMut$X8r[m_list] <-allMut$X8[m_list] * GC_centiles.results$r8[i]
  allMut$X9r[m_list] <-allMut$X9[m_list] * GC_centiles.results$r9[i]
  allMut$X10r[m_list] <-allMut$X10[m_list] * GC_centiles.results$r10[i]
  allMut$X11r[m_list] <-allMut$X11[m_list] * GC_centiles.results$r11[i]
  allMut$X12r[m_list] <-allMut$X12[m_list] * GC_centiles.results$r12[i]
  allMut$X13r[m_list] <-allMut$X13[m_list] * GC_centiles.results$r13[i]
  allMut$X14r[m_list] <-allMut$X14[m_list] * GC_centiles.results$r14[i]
  allMut$X15r[m_list] <-allMut$X15[m_list] * GC_centiles.results$r15[i]
  allMut$X16r[m_list] <-allMut$X16[m_list] * GC_centiles.results$r16[i]
  allMut$X17r[m_list] <-allMut$X17[m_list] * GC_centiles.results$r17[i]
  allMut$X18r[m_list] <-allMut$X18[m_list] * GC_centiles.results$r18[i]
  allMut$X19r[m_list] <-allMut$X19[m_list] * GC_centiles.results$r19[i]
}

allMut$DAF_pi_corr <- (allMut$X1*19)/190 + 
  (allMut$X2r*18*2)/190 + 
  (allMut$X3r*17*3)/190 + 
  (allMut$X4r*16*4)/190 + 
  (allMut$X5r*15*5)/190 + 
  (allMut$X6r*14*6)/190 + 
  (allMut$X7r*13*7)/190 +
  (allMut$X8r*12*8)/190 +
  (allMut$X9r*11*9)/190 +
  (allMut$X10r*10*10)/190 +
  (allMut$X11r*9*11)/190 +
  (allMut$X12r*8*12)/190 +
  (allMut$X13r*7*13)/190 +
  (allMut$X14r*6*14)/190 +
  (allMut$X15r*5*15)/190 +
  (allMut$X16r*4*16)/190 +
  (allMut$X17r*3*17)/190 +
  (allMut$X18r*2*18)/190 +
  (allMut$X19r*19)/190

daf_pi_corr_class4 <- aggregate(DAF_pi_corr~centile+class4cat+population, allMut, sum)


#Different data frames of parameters of interest
all_freq_class4<- aggregate((allMut$X1+allMut$X2+allMut$X3+allMut$X4+allMut$X5+allMut$X6+allMut$X7+allMut$X8+allMut$X9+allMut$X10+allMut$X11+allMut$X12+allMut$X13+allMut$X14+allMut$X15+allMut$X16+allMut$X17+allMut$X18+allMut$X19)~centile+class4cat+population, allMut, sum)


colnames(all_freq_class4) <- c("centile", "class4cat", "population", "S")
watterson_class4 <- all_freq_class4$S/sum(1/seq(1:19))
all_freq_class4$watterson<- watterson_class4

singleton_freq_class4<- aggregate((allMut$X1)~centile+class4cat+population, allMut, sum)
colnames(singleton_freq_class4) <- c("centile", "class4cat", "population", "X1")


gc_class4 <- aggregate(localGC~centile+class4cat+population, allMut, mean)

gc_class4$gcDist_from_eq <- NA
gc_class4$LiBulmerGCeq <- NA

for (i in 1:length(GC_centiles.results$gcDist_from_eq)){
  matches <- (GC_centiles.results$region[i] == gc_class4$centile & GC_centiles.results$population[i] == gc_class4$population)
  m_list <- (1:length(matches))[matches]
  gc_class4$gcDist_from_eq[m_list] <- GC_centiles.results$gcDist_from_eq[i]
  gc_class4$LiBulmerGCeq[m_list] <- GC_centiles.results$LiBulmerGCeq[i]
  }


cds_class4 <- aggregate(cds_total/(invar_L+unanchor_L+anchor_L)~centile+class4cat+population, allMut, mean)
cds_anchor_class4 <- aggregate(cds_anchor/anchor_L~centile+class4cat+population, allMut, sum)

L_class4_frW <- aggregate(((invar_L+unanchor_L+anchor_L)-(invar_GC+unanchor_GC+anchor_GC))~centile+class4cat+population, allMut[c(grep("A_", allMut$class), grep("T_", allMut$class)), ], mean)
L_class4_frS <- aggregate((invar_GC+unanchor_GC+anchor_GC)~centile+class4cat+population, allMut[c(grep("C_", allMut$class), grep("G_", allMut$class)), ], mean)

colnames(L_class4_frW) <- c("centile", "class4cat", "population", "L")
colnames(L_class4_frS) <- c("centile", "class4cat", "population", "L")
L_class4 <- rbind(L_class4_frS, L_class4_frW)
L_class4 <- L_class4[order(L_class4[,3], L_class4[,2]),]

#Tajima's D calculations
#The sample size
n=20

#Notations here follow Tajima (1989)
a1 <- sum(1/seq(1:(n-1)))
a2 <- sum(1/(seq(1:(n-1))^2))
b1 <- (n+1)/(3*(n-1))
b2 <- (2*((n^2)+n+3))/(9*n*(n-1))
c1 <- b1 - 1/a1
c2 <- b2 - (n+2)/(a1*n) +(a2/(a1^2))
e1 <- c1/a1
e2 <- c2/((a1^2)+a2)

tajimas_D_class4 <- (daf_pi_class4$DAF_pi - all_freq_class4$watterson)/sqrt(e1*all_freq_class4$S+e2*all_freq_class4$S*(all_freq_class4$S+1))
all_freq_class4$tajimasD <- tajimas_D_class4


col3 <- wes_palette("Cavalcanti1")


#GC vs π - class4cat
par(mfrow=c(1,1))
pop="Swe-sin" #can be any here
ggplot(daf_pi_class4[daf_pi_class4$population == pop,], aes(y=daf_pi_class4[daf_pi_class4$population == pop,4]/L_class4[L_class4$population == pop,4], x=gc_class4[daf_pi_class4$population == pop,4], col=class4cat, shape=class4cat))+geom_point(size=2)+
  
  scale_shape_manual(name="", values = c(SS=4, SW=16, WS=17,  WW=3)) +
  scale_color_manual(name="", values = c(SS=col3[1], SW=col3[2], WS=col3[5],  WW=col3[4])) + 
  ylab(expression(pi))+
  xlab("Observed GC content")+
  theme_classic()+
  theme(legend.position = c(0.9,0.83), panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

#GC vs π - class4cat - demography corrected
ggplot(daf_pi_corr_class4[daf_pi_corr_class4$population == pop,], aes(y=daf_pi_corr_class4[daf_pi_corr_class4$population == pop,4]/L_class4[L_class4$population == pop,4], x=gc_class4[daf_pi_class4$population == pop,4], col=class4cat, shape=class4cat))+geom_point(size=2)+
  
  scale_shape_manual(name="", values = c(SS=4, SW=16, WS=17,  WW=3)) +
  scale_color_manual(name="", values = c(SS=col3[1], SW=col3[2], WS=col3[5],  WW=col3[4])) + 
  ylab(expression(pi))+
  xlab("Observed GC content")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))






#### Coplot prep ####

treated_allMut <- cbind(daf_pi_class4, L=L_class4[4], cds_density=cds_class4[,4], localGC=gc_class4[,4], DAF_pi_corr=daf_pi_corr_class4[,4], watterson=all_freq_class4[,5], tajD= all_freq_class4[,6])
treated_allMut$watterson_per_bp <- treated_allMut$watterson/treated_allMut$L
treated_allMut$pi <- treated_allMut$DAF_pi/treated_allMut$L
treated_allMut$pi_corr <- treated_allMut$DAF_pi_corr/treated_allMut$L







#### ------- Full coplot fig ---------####
### GC vs CDS ###
pop="Swe-sin" #Can be any here
category="SS" #Can be any here
treated_allMut_cat_pop<- treated_allMut[treated_allMut$population == pop,]
treated_allMut_GC_split <- split(treated_allMut[treated_allMut$class4cat== category & treated_allMut$population == pop,], cut(treated_allMut[treated_allMut$class4cat == category & treated_allMut$population == pop, 7], breaks = 5))

regInfo = function(x, y, ...) {
  if(lmp(lm(y ~ x))<0.05){
    abline(lm(y ~ x))
    r2 = round(summary(lm(y ~ x))$adj.r.squared, digits = 2)
    my.s = round(summary(lm(y ~ x))$coefficients[2,1], digits = 2) #$coefficients[2,1]*1000
    mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    rp = vector('expression',2)
    rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                       list(MYVALUE = format(r2,dig=2)))[2]
    rp[2] = substitute(expression(italic(k) == MYOTHERVALUE), 
                       list(MYOTHERVALUE = format(my.s, digits = 2)))[2]
    legend('topright', legend = rp, bty = 'n')
  }
}


par(mfrow=c(5,4), mai = c(0.1, 0.3, 0.3, 0.1))
for(i in 1:4){
  if (i == 1){par(mai = c(0.1, 0.3, 0.3, 0.1))
    yAX="s"
  }
  main_text=paste("GC: ",names(treated_allMut_GC_split)[i], sep="")
  plot(treated_allMut_GC_split[[i]]$cds_density, treated_allMut_GC_split[[i]]$localGC, xlim = c(min(treated_allMut_cat_pop$cds_density), max(treated_allMut_cat_pop$cds_density)), ylim = c(min(treated_allMut_cat_pop$localGC), max(treated_allMut_cat_pop$localGC)), pch=1, col="black", yaxt=yAX, ylab = "Observed GC content", xaxt="s", main=main_text)
  grid(5, 5, lwd = 1)
  regInfo(x=treated_allMut_GC_split[[i]]$cds_density, y=treated_allMut_GC_split[[i]]$localGC)
  yAX="n"
  par(mai = c(0.1, 0.1, 0.3, 0.1))
}


### CDS vs pi ###
regInfo = function(x, y, ...) {
  if(lmp(lm(y ~ x))<0.05){
    abline(lm(y ~ x))
    r2 = round(summary(lm(y ~ x))$adj.r.squared, digits = 2)
    my.s = round(summary(lm(y ~ x))$coefficients[2,1]*1000, digits = 2)
    mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    rp = vector('expression',2)
    rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                       list(MYVALUE = format(r2,dig=2)))[2]
    rp[2] = substitute(expression(italic(k) == MYOTHERVALUE), 
                       list(MYOTHERVALUE = format(my.s, digits = 2)))[2]
    legend('topright', legend = rp, bty = 'n')
  }
}


#par(mfrow=c(5,4), mai = c(0.1, 0.3, 0.1, 0.1))
#pop="Swe-sin"
for (category in c("SS","WW","SW","WS")){
  if(category == "SS"){catcol=col3[1]
  par(mai = c(0.1, 0.1, 0.3, 0.1))
  }
  else if(category == "WW"){catcol=col3[4]}
  else if(category == "SW"){catcol=col3[2]}
  else{catcol=col3[5]}
  treated_allMut_cat_pop<- treated_allMut[treated_allMut$population == pop,]
  treated_allMut_GC_split <- split(treated_allMut[treated_allMut$class4cat== category & treated_allMut$population == pop,], cut(treated_allMut[treated_allMut$class4cat == category & treated_allMut$population == pop, 7], breaks = 5))
  for(i in 1:4){
    if (category != "WS" & category != "SS" ){
      par(mai = c(0.1, 0.1, 0.1, 0.1))
      if (i == 1){par(mai = c(0.1, 0.3, 0.1, 0.1))
        yAX="s"
      }
        plot(treated_allMut_GC_split[[i]]$cds_density, treated_allMut_GC_split[[i]]$pi, xlim = c(min(treated_allMut_cat_pop$cds_density), max(treated_allMut_cat_pop$cds_density)), ylim = c(min(treated_allMut_cat_pop$pi), max(treated_allMut_cat_pop$pi)), pch=16, col=catcol, yaxt=yAX, xaxt="n", ann=F)
        regInfo(x=treated_allMut_GC_split[[i]]$cds_density, y=treated_allMut_GC_split[[i]]$pi)
        grid(5, 5, lwd = 1)
        yAX="n"
    }
    else if(category == "SS"){
      par(mai = c(0.1, 0.1, 0.3, 0.1))
      if (i == 1){par(mai = c(0.1, 0.3, 0.3, 0.1))
        yAX="s"
      }
      #main_text=paste("GC: ",names(treated_allMut_GC_split)[i], sep="")
      #plot(treated_allMut_GC_split[[i]]$cds_density, treated_allMut_GC_split[[i]]$pi, xlim = c(min(treated_allMut_cat_pop$cds_density), max(treated_allMut_cat_pop$cds_density)), ylim = c(min(treated_allMut_cat_pop$pi), max(treated_allMut_cat_pop$pi)), pch=16, col=catcol, yaxt=yAX, xaxt="n", main=main_text)
      plot(treated_allMut_GC_split[[i]]$cds_density, treated_allMut_GC_split[[i]]$pi, xlim = c(min(treated_allMut_cat_pop$cds_density), max(treated_allMut_cat_pop$cds_density)), ylim = c(min(treated_allMut_cat_pop$pi), max(treated_allMut_cat_pop$pi)), pch=16, col=catcol, yaxt=yAX, xaxt="n")
      regInfo(x=treated_allMut_GC_split[[i]]$cds_density, y=treated_allMut_GC_split[[i]]$pi)
      grid(5, 5, lwd = 1)
      yAX="n"
      }
    else{
      par(mai = c(0.3, 0.1, 0.1, 0.1))
      if (i == 1){par(mai = c(0.3, 0.3, 0.1, 0.1))
        yAX="s"
        }
      plot(treated_allMut_GC_split[[i]]$cds_density, treated_allMut_GC_split[[i]]$pi, xlim = c(min(treated_allMut_cat_pop$cds_density), max(treated_allMut_cat_pop$cds_density)), ylim = c(min(treated_allMut_cat_pop$pi), max(treated_allMut_cat_pop$pi)), pch=16, col=catcol, yaxt=yAX, xaxt="s", ann=F)
      regInfo(x=treated_allMut_GC_split[[i]]$cds_density, y=treated_allMut_GC_split[[i]]$pi)
      grid(5, 5, lwd = 1)
      yAX="n"
    }

  }
  
}
dev.off()




#### ALLMUT - Multiple Linear Regression model ####
library("car")


treated_allMut_cat_pop<- treated_allMut[treated_allMut$population == pop,]
treated_allMut_GC_split <- split(treated_allMut[treated_allMut$class4cat== category & treated_allMut$population == pop,], cut(treated_allMut[treated_allMut$class4cat == category & treated_allMut$population == pop, 7], breaks = 5))

pop="Spa-sin"
categ="WS"
summary(lm(pi~poly(cds_density,2) + poly(localGC,2) + localGC:cds_density, data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))
summary(lm(pi~poly(cds_density,1) + poly(localGC,1), data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))
summary(lm(pi~ poly(localGC,1), data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))

summary(lm(pi~ poly(localGC,1), data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))


model_full <- (lm(pi~poly(cds_density,2) + poly(localGC,2) + localGC:cds_density , data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))
model_sparse <- (lm(pi~poly(cds_density,1) + poly(localGC,2), data=treated_allMut[treated_allMut$population == pop & treated_allMut$class4cat == categ,]))

vif(model_sparse)
vif(model_full)



#### Centile read coverage #### 
# Created 14/5 -2020 
centile_cov<- read.csv(file = file.choose() , header = T, sep = ";", dec =".", stringsAsFactors = F)

centile_cov$Centile_num <- as.integer(gsub("C_([0-9]+)", replacement = "\\1", centile_cov$Centile))

centile_cov$localGC <- NA
centile_cov$pi <- NA

for (i in 1:length(GC_centiles$region)){
  matches <- (GC_centiles$region[i] == centile_cov$Centile & GC_centiles$population[i] == centile_cov$Population)
  m_list <- (1:length(matches))[matches]
  centile_cov$localGC[m_list] <- GC_centiles$localGC[i]
  centile_cov$pi[m_list] <- GC_centiles$total_pi[i]
}

#GC content vs coverage
ggplot(centile_cov, aes(x=localGC, y=Coverage, col = Population)) +
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) +
  geom_point(alpha=0.2)+
  ylim(0,40)+
  geom_smooth(method="loess")+
  xlab("GC content")+
  ylab("Coverage (average reads per bp)")+
  theme_classic() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

#π vs coverage
ggplot(centile_cov, aes(x=pi, y=Coverage, col = Population)) +
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) +
  geom_point(alpha=0.2)+
  ylim(0,40)+
  geom_smooth(method = "lm", formula = y~poly(x,2))+
  xlab(expression(pi))+
  ylab("Coverage (average reads per bp)")+
  theme_classic() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))




