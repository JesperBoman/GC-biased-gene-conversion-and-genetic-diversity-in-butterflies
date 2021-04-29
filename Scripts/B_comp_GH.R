# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# Comparative analyses of GC-biased gene conversion (B)
# ==========================================================================================================
# Jesper Boman                      5 mar 2021
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #




#### gBGC point estimates LRT tests ####
gBGC_genomewide_stats <- read.csv(file = file.choose() , header = T, sep = ";", dec =".", stringsAsFactors = F)

gBGC_genomewide_stats

library(wesanderson)
library(ggplot2)
col <- wes_palette("Darjeeling1")
col2 <-wes_palette("Darjeeling2")


#M1 vs M0
pchisq(-2*(gBGC_genomewide_stats[grep("Swe-sin", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Swe-sin", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Spa-sin", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Spa-sin", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Kaz-sin", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Kaz-sin", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Kaz-juv", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Kaz-juv", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Ire-juv", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Ire-juv", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Spa-rea", gBGC_genomewide_stats$Population),]$lnL0-gBGC_genomewide_stats[grep("Spa-rea", gBGC_genomewide_stats$Population),]$lnL1), df=1, lower.tail=F)
#

#M1* vs M1
pchisq(-2*(gBGC_genomewide_stats[grep("Swe-sin", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Swe-sin", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Spa-sin", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Spa-sin", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Kaz-sin", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Kaz-sin", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Kaz-juv", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Kaz-juv", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Ire-juv", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Ire-juv", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
pchisq(-2*(gBGC_genomewide_stats[grep("Spa-rea", gBGC_genomewide_stats$Population),]$lnL1-gBGC_genomewide_stats[grep("Spa-rea", gBGC_genomewide_stats$Population),]$lnL1cor), df=3, lower.tail=F)
#

#Examples of data exploration (Comparing different filter sets)
t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$eWS, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$eWS, paired = T)
t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$eSW, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$eSW, paired = T)

t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$B_M1, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$B_M1, paired = T)
t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$Bcor_M1, paired = T)

t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_noCpG",]$B_M1, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$B_M1, paired = T)
t.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_noCpG",]$Bcor_M1, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$Bcor_M1, paired = T)


1/(1+gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$lambdacor_M1*exp(-gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_lib",]$Bcor_M1))

gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$lambdacor_M1-gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_noCpG",]$lambdacor_M1
gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$lambdacor_M1-gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_onlyCpG",]$lambdacor_M1
gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_noCpG",]$lambdacor_M1-gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons_onlyCpG",]$lambdacor_M1


#### Genomewide - bootstrapped sites, no exons ####
swe_sin_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
swe_sin_gw_bs$Population <- "Swe_sin"

spa_sin_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_sin_gw_bs$Population <- "Spa_sin"

kaz_sin_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_sin_gw_bs$Population <- "Kaz_sin"

kaz_juv_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_juv_gw_bs$Population <- "Kaz_juv"

ire_juv_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
ire_juv_gw_bs$Population <- "Ire_juv"

spa_rea_gw_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_rea_gw_bs$Population <- "Spa_rea"

pchisq(-2*(mean(swe_sin_gw_bs$lnL1)-mean(swe_sin_gw_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_sin_gw_bs$lnL1)-mean(spa_sin_gw_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_sin_gw_bs$lnL1)-mean(kaz_sin_gw_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_juv_gw_bs$lnL1)-mean(kaz_juv_gw_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(ire_juv_gw_bs$lnL1)-mean(ire_juv_gw_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_rea_gw_bs$lnL1)-mean(spa_rea_gw_bs$lnL1cor)), df=3, lower.tail=F)

comb_gw_bs <- rbind(swe_sin_gw_bs, spa_sin_gw_bs, kaz_sin_gw_bs, kaz_juv_gw_bs, ire_juv_gw_bs, spa_rea_gw_bs)


leptidea_Bcor_M1_estimates_NE <- data.frame(Population = unique(factor(comb_gw_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin"))), Bcor_M1 =gBGC_genomewide_stats[gBGC_genomewide_stats$Filter=="noExons",]$Bcor_M1 )

"6.04 x 5"
ggplot(comb_gw_bs, aes(x=factor(comb_gw_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=Bcor_M1, fill = Population)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + theme_classic() + 
  xlab("Population") +
  ylab("B") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text=element_text(size=11, colour="black"), axis.title=element_text(size=14)) +
  labs(title="All non-exonic sites") +
  ylim(0,1) +
  geom_violin(data=comb_gw_bs, aes(x=factor(comb_gw_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=B_M1, fill=Population), alpha=0.5, draw_quantiles = 0.5) +
  geom_violin(draw_quantiles = 0.5)+
  geom_point(data=leptidea_Bcor_M1_estimates_NE, shape=1, aes(x=Population, y=Bcor_M1))+
  scale_x_discrete(labels=c("Ire-juv","Kaz-juv", "Spa-rea", "Kaz-sin", "Swe-sin", "Spa-sin"))+
  geom_text(x=6.17, y=0.9, label="B M1*")+
  annotate(geom="point", x=5.8, y=0.9, shape=15, size = 5)+
  geom_text(x=6.15, y=0.85, label="B M1")+
  annotate(geom="point", x=5.8, y=0.85, shape=15, alpha=0.5, size = 5)

#### Genomewide - bootstrapped sites, no exons no CpG ####
swe_sin_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
swe_sin_gw_filt_bs$Population <- "Swe_sin"

spa_sin_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_sin_gw_filt_bs$Population <- "Spa_sin"

kaz_sin_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_sin_gw_filt_bs$Population <- "Kaz_sin"

kaz_juv_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_juv_gw_filt_bs$Population <- "Kaz_juv"

ire_juv_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
ire_juv_gw_filt_bs$Population <- "Ire_juv"

spa_rea_gw_filt_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_rea_gw_filt_bs$Population <- "Spa_rea"


pchisq(-2*(mean(swe_sin_gw_filt_bs$lnL1)-mean(swe_sin_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_sin_gw_filt_bs$lnL1)-mean(spa_sin_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_sin_gw_filt_bs$lnL1)-mean(kaz_sin_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_juv_gw_filt_bs$lnL1)-mean(kaz_juv_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(ire_juv_gw_filt_bs$lnL1)-mean(ire_juv_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_rea_gw_filt_bs$lnL1)-mean(spa_rea_gw_filt_bs$lnL1cor)), df=3, lower.tail=F)

comb_gw_filt_bs <- rbind(swe_sin_gw_filt_bs, spa_sin_gw_filt_bs, kaz_sin_gw_filt_bs, kaz_juv_gw_filt_bs, ire_juv_gw_filt_bs, spa_rea_gw_filt_bs)

leptidea_Bcor_M1_estimates_noCpG <- data.frame(Population = unique(factor(comb_gw_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin"))), Bcor_M1 =gBGC_genomewide_stats[gBGC_genomewide_stats$Filter=="noExons_noCpG",]$Bcor_M1 )


ggplot(comb_gw_filt_bs, aes(x=factor(comb_gw_filt_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=Bcor_M1, fill = Population)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + theme_classic() + 
  xlab("Population") +
  ylab("B") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text=element_text(size=11, colour="black"), axis.title=element_text(size=14)) +
  labs(title="No ancestral CpG-prone sites") +
  ylim(0,1.5) +
  geom_violin(data=comb_gw_filt_bs, aes(x=factor(comb_gw_filt_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=B_M1, fill=Population, alpha=Population), draw_quantiles = 0.5) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(data=leptidea_Bcor_M1_estimates_noCpG, shape=1, aes(x=Population, y=Bcor_M1))+
  scale_x_discrete(labels=c("Ire-juv","Kaz-juv", "Spa-rea", "Kaz-sin", "Swe-sin", "Spa-sin"))


#### Genomewide - bootstrapped sites, no exons CpGonly ####
swe_sin_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
swe_sin_gw_CpG_bs$Population <- "Swe_sin"

spa_sin_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_sin_gw_CpG_bs$Population <- "Spa_sin"

kaz_sin_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_sin_gw_CpG_bs$Population <- "Kaz_sin"

kaz_juv_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
kaz_juv_gw_CpG_bs$Population <- "Kaz_juv"

ire_juv_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
ire_juv_gw_CpG_bs$Population <- "Ire_juv"

spa_rea_gw_CpG_bs <- read.table(file = file.choose() , header = T, dec =".", stringsAsFactors = F)
spa_rea_gw_CpG_bs$Population <- "Spa_rea"

pchisq(-2*(mean(swe_sin_gw_CpG_bs$lnL1)-mean(swe_sin_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_sin_gw_CpG_bs$lnL1)-mean(spa_sin_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_sin_gw_CpG_bs$lnL1)-mean(kaz_sin_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(kaz_juv_gw_CpG_bs$lnL1)-mean(kaz_juv_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(ire_juv_gw_CpG_bs$lnL1)-mean(ire_juv_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)
pchisq(-2*(mean(spa_rea_gw_CpG_bs$lnL1)-mean(spa_rea_gw_CpG_bs$lnL1cor)), df=3, lower.tail=F)


comb_gw_CpG_bs <- rbind(swe_sin_gw_CpG_bs, spa_sin_gw_CpG_bs, kaz_sin_gw_CpG_bs, kaz_juv_gw_CpG_bs, ire_juv_gw_CpG_bs, spa_rea_gw_CpG_bs)

#gBGC_genomewide_stats <- read.csv(file = file.choose() , header = T, sep = ";", dec =".", stringsAsFactors = F)
leptidea_Bcor_M1_estimates_CpG <- data.frame(Population = unique(factor(comb_gw_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin"))), Bcor_M1 =gBGC_genomewide_stats[gBGC_genomewide_stats$Filter=="noExons_onlyCpG",]$Bcor_M1 )

"6.04 x 5"
ggplot(comb_gw_CpG_bs, aes(x=factor(comb_gw_CpG_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=Bcor_M1, fill = Population)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + theme_classic() + 
  xlab("Population") +
  ylab("B") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text=element_text(size=11, colour="black"), axis.title=element_text(size=14)) +
  labs(title="Ancestral CpG-prone sites") +
  ylim(0,1.5) +
  geom_violin(data=comb_gw_CpG_bs, aes(x=factor(comb_gw_CpG_bs$Population, levels=c("Ire_juv","Kaz_juv", "Spa_rea", "Kaz_sin", "Swe_sin", "Spa_sin")), y=B_M1, fill=Population, alpha=Population), draw_quantiles = 0.5) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(data=leptidea_Bcor_M1_estimates_CpG, shape=1, aes(x=Population, y=Bcor_M1))+
  scale_x_discrete(labels=c("Ire-juv","Kaz-juv", "Spa-rea", "Kaz-sin", "Swe-sin", "Spa-sin"))


#### DAF spectrum plot ####

DAF_swe_sin <- read.table(file.choose(), stringsAsFactors = F, header=T, dec=",")


ggplot(data=DAF_swe_sin, aes(x=Freq_cat_num, y=Freq, fill=Class)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=c("Black", "Grey", "White"))+
  theme_classic()+
  ylab("Density") +
  xlab("Derived allele count") +
  labs(title="All non-exonic sites: Swe-sin") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))



#### Leptidea phylogeny ####


#Read in phylogeny
library(ape)
leptidea_phylogeny <- read.tree(file = file.choose())
plot(leptidea_phylogeny)
tips <-  leptidea_phylogeny$tip.label
tips <- tips[-(grep("leptidea", tips))]

leptidea_phylogeny_pruned <- drop.tip(leptidea_phylogeny, leptidea_phylogeny$tip.label[match(as.vector(tips), leptidea_phylogeny$tip.label)])
tips <-  leptidea_phylogeny_pruned$tip.label
leptidea_phylogeny_pruned$tip.label <- c("Leptidea juvernica - Ireland", "Leptidea juvernica - Kazakhstan", "Leptidea reali - Spain", "Leptidea sinapis - Spain", "Leptidea sinapis - Kazakhstan", "Leptidea sinapis - Sweden" )

#Plot phylogeny
plot(leptidea_phylogeny_pruned, show.node.label =T, direction="rightwards", tip.color=c(col[2], col[5], col2[2],  col[1], col[4], col[3]))
plot(leptidea_phylogeny_pruned)
leptidea_phylogeny_pruned$tip.label <- c("Ire-juv", "Kaz-juv", "Spa-rea", "Spa-sin", "Kaz-sin", "Swe-sin" )





#### Leptidea phylostats ####

library("ape")

### Population comparative (dependent on pi) ####
#Needs GC_centiles.results from the script GC_centiles.R, or just add it to the gBGC_genomewide_stats input file
genomewide_pi <- cbind(aggregate((pi_sum_unanchor+pi_sum_anchor)~population, GC_centiles.results, sum), aggregate((invar_L+unanchor_L+anchor_L)~population, GC_centiles.results, sum)[,2])
colnames(genomewide_pi) <- c("Population", "Pi_sum", "L")
genomewide_pi$pi <- genomewide_pi$Pi_sum/genomewide_pi$L

gBGC_genomewide_stats <- read.csv(file = file.choose() , header = T, sep = ";", dec =".", stringsAsFactors = F)
gBGC_genomewide_stats$Chr_num_mean <- (gBGC_genomewide_stats$Chr_num_high+gBGC_genomewide_stats$Chr_num_low)/2

cor.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$pi, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1)
cor.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Chr_num_low, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1)
cor.test(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Chr_num_high, gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1)


#Plot π vs B - N.B needs plot p2 from below for inset
"6.04 x 5"
ggplot(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",], aes(x=pi, y=Bcor_M1, col=Population)) + geom_point(size=4) + 
  ylim(0,0.85)+
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  ylab("B") +
  xlab(expression(pi)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14)) +
  annotation_custom(ggplotGrob(p2), xmin = 0.005, xmax = 0.00703, 
                    ymin = 0.45, ymax = 0.908)

#Plot chromosome number vs B- N.B needs plot p3 from below for inset
"6.04 x 5"
ggplot(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",], aes(x=Chr_num_low, y=Bcor_M1, col=Population)) + geom_point(size=4) + geom_point(size=4, aes(x=Chr_num_high))+  
  ylim(0,0.85)+
  scale_color_manual(name="Population", values = c(col[2], col[5], col[4],  col2[2], col[1], col[3])) + 
  theme_classic() + 
  xlim(50, 125)+
  ylab("B") +
  xlab("Diploid chromosome number") +
  geom_segment(aes(x=Chr_num_low, y=Bcor_M1, xend=Chr_num_high, yend=Bcor_M1), lty=3) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position ="none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))+
  annotation_custom(ggplotGrob(p3), xmin = 94, xmax = 129.9, 
                  ymin = 0.45, ymax = 0.908)


#Produce phylogenetic independent contrasts using the method described by Felsenstein (1985)
pic(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1,leptidea_phylogeny_pruned)
pic(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$pi,leptidea_phylogeny_pruned)

pipic <- as.data.frame(pic(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$pi,leptidea_phylogeny_pruned, var.contrasts = T, scaled=F))
chrpic <- as.data.frame(pic(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Chr_num_mean,leptidea_phylogeny_pruned, var.contrasts = T, scaled=F))
Bpic <- as.data.frame(pic(gBGC_genomewide_stats[gBGC_genomewide_stats$Filter == "noExons",]$Bcor_M1,leptidea_phylogeny_pruned, var.contrasts = T, scaled=F))


summary(lm(scale(Bpic$contrasts)~scale(pipic$contrasts)))
summary(lm(scale(Bpic$contrasts)~scale(chrpic$contrasts)))
summary(lm(scale(Bpic$contrasts)~scale(chrpic$contrasts)))
summary(lm(Bpic$contrasts~chrpic$contrasts+pipic$contrasts))

p2 <- ggplot(Bpic, aes(x=pipic$contrasts, y=Bpic$contrasts))+geom_point(size=3)+
  theme_classic()+
  xlim(-0.0022, 0.0025)+
  xlab("") +
  ylab("") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position ="none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=8, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

p3 <- ggplot(Bpic, aes(x=chrpic$contrasts, y=Bpic$contrasts))+geom_point(size=3)+
  theme_classic()+
  xlim(-50,40)+
  xlab("") +
  ylab("") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position ="none", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=8, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))



#Same procedure as above but excluding the Spanish L. sinapis population
gBGC_genomewide_stats_noSpaSin <- gBGC_genomewide_stats[gBGC_genomewide_stats$Population != "Spa-sin",]

cor.test(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Chr_num_low, gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Bcor_M1)
cor.test(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Chr_num_high, gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Bcor_M1)
cor.test(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Chr_num_mean, gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Bcor_M1)

leptidea_phylogeny_pruned_noSpasin <- drop.tip(leptidea_phylogeny_pruned,"Spa-sin")

plot(leptidea_phylogeny_pruned_noSpasin)

pipic_noSpasin <- as.data.frame(pic(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$pi,leptidea_phylogeny_pruned_noSpasin, var.contrasts = T, scaled=F))
chrpic_noSpasin <- as.data.frame(pic(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Chr_num_mean,leptidea_phylogeny_pruned_noSpasin, var.contrasts = T, scaled=F))
Bpic_noSpasin <- as.data.frame(pic(gBGC_genomewide_stats_noSpaSin[gBGC_genomewide_stats_noSpaSin$Filter == "noExons",]$Bcor_M1,leptidea_phylogeny_pruned_noSpasin, var.contrasts = T, scaled=F))

cor.test(chrpic_noSpasin$contrasts, Bpic_noSpasin$contrasts)
summary(lm(scale(Bpic_noSpasin$contrasts)~scale(chrpic_noSpasin$contrasts)))
summary(lm(scale(Bpic_noSpasin$contrasts)~scale(chrpic_noSpasin$contrasts)+scale(pipic_noSpasin$contrasts)))
plot(chrpic_noSpasin$contrasts, Bpic_noSpasin$contrasts)
