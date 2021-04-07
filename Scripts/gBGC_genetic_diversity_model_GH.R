# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A model for the effect of gBGC and mutation bias on genetic diversity
# ==========================================================================================================
# Jesper Boman                      5 mar 2021
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Colour for plots
library(wesanderson)
col <- wes_palette("Darjeeling1")
col2 <-wes_palette("Darjeeling2")
col3 <- wes_palette("Cavalcanti1")


#set parameter values for mutation bias and gBGC
lambda<-list(1,1.2,1.8,2,3)
B_seq<-seq(0,8,0.01)

# compute equilibrium GC content based on a mutation-selection balanace
x_gc<-lapply(lambda, function (x) 1/(1+x*exp(-B_seq)))

pi<-x_gc
piSW<-x_gc
piWS<-x_gc
piWW <- x_gc
piSS <- x_gc

#This parameter (a) represents the ratio of N-->N to W-->S diversity
a=1

#Standardization function
pi_0<-lapply(lambda, function (l,a=1) (2*l+a*(1+l))/(1+l) )


for(ii in 1:length(lambda)){
  
  pi[[ii]]<-((x_gc[[ii]])*lambda[[ii]]*2*(1/(1-exp(B_seq))+1/(B_seq)) +
               (1-(x_gc[[ii]]))*2*(1/(1-exp(-B_seq))-1/(B_seq)) + a) 

  piSW[[ii]]<-lambda[[ii]]*2*(1/(1-exp(B_seq))+1/(B_seq))
  piWS[[ii]]<-2*(1/(1-exp(-B_seq))-1/(B_seq))
  piWW[[ii]]<-rep(1, times=801)
  piSS[[ii]]<-rep(1, times=801)
  
  pi[[ii]]<-pi[[ii]]/pi_0[[ii]]
  piSW[[ii]]<-piSW[[ii]]/pi_0[[ii]]
  piWS[[ii]]<-piWS[[ii]]/pi_0[[ii]]
  piWW[[ii]]<-piWW[[ii]]/pi_0[[ii]]
  piSS[[ii]]<-piSS[[ii]]/pi_0[[ii]]
}


#Main model plot: GC vs π-rel
#We are basically varying the equilibrium GC content by varying the strength of B, since the equilibrium GC content is described by the balance of B and lambda according to the Li-Bulmer equation. 
#See Boman et al. 2021 for more details.

plot(x_gc[[5]],pi[[5]],type="l",xlim=c(0,1), col=col[1], ylim =c(0,1.5), lwd=3,  cex.lab=1.3,  xlab = "Equilibrium GC content", ylab = expression(pi~ ""["rel"]), main = "gBGC-mutation-drift model") 
lines(x_gc[[1]],pi[[1]],type="l",lty=1, lwd=3)
lines(x_gc[[4]],pi[[4]],type="l",lty=1, col=col[4], lwd=3)
abline(h=1,col="grey",lty=2)
legend("topright", legend =c(expression(lambda~ " = 3"), expression(lambda~ " = 2"), expression(lambda~ " = 1")), col=c(col[1], col[4], "black"), lty=1, lwd=3)


#Plot of B vs π-rel
plot(B_seq,pi[[5]],type="l",ylim=c(0,1.5))
abline(h=1,col="grey",lty=2)



#Plot for separate mutation categories
plot(x_gc[[5]], pi[[5]], type="l",xlim=c(0,1), col=col[1], ylim =c(0,1.5), lwd=3, lty=2,  cex.lab=1.3,  xlab = "Equilibrium GC content", ylab = expression(pi~ ""["rel"]), main = "gBGC-mutation-drift model") 
lines(x_gc[[5]], piSW[[5]], col=col3[2], lwd=3)
lines(x_gc[[5]], piWS[[5]], col=col3[5], lwd=3)
lines(x_gc[[5]], piWW[[5]], col=col3[4], lwd=3, lty=1)
lines(x_gc[[5]], piSS[[5]], col=col3[1], lwd=3, lty=1)
abline(h=1,col="grey",lty=2)
legend("topright", legend =c("All","SS", "SW", "WS", "WW"), col=c(col[1], col3[1], col3[2], col3[5], col3[4]), lwd=3, pch=c(NA,16,17,3,4), lty=c(2,rep(1,4)))
points(x_gc[[5]][seq(2, length(x_gc[[5]]), 60)], piSW[[5]][seq(2, length(x_gc[[5]]), 60)], col=col3[2], pch=16) #SW
points(x_gc[[5]][seq(2, length(x_gc[[5]]), 60)], piWS[[5]][seq(2, length(x_gc[[5]]), 60)], col=col3[5], pch=17) #WS
points(x_gc[[5]][seq(2, length(x_gc[[5]]), 60)], piWW[[5]][seq(2, length(x_gc[[5]]), 60)], col=col3[4], pch=3)  #WW
points(x_gc[[5]][seq(2, length(x_gc[[5]]), 83)], piSS[[5]][seq(3, length(x_gc[[5]]), 83)], col=col3[1], pch=4)  #SS









##### Extra: Investigating the effect of random noise for B and lambda using empirical values ####

#set parameter values for mutation bias and gBGC
lambda_mean<-list(1,1.2,1.8,2,3)
B_seq_mean<-seq(0,8,0.01)



#Here we use the standard deviation from our empirical data (see GC_centiles script)
#You could e.g. use variation in B and lambda among chromosomes or genomic windows
lambda_stdev=sd(GC_centiles.results[GC_centiles.results$population == "Ire-juv", ]$lambda_M1)
B_stdev=sd(GC_centiles.results[GC_centiles.results$population == "Ire-juv", ]$B_M1)



x_gc<-lapply(lambda, function (x) 1/(1+x*exp(-B_seq)))


#Sampling of B and lambda from normal distributions
lambda_rand <- rnorm(length(B_seq), mean=3, sd = lambda_stdev)
B_seq_rand <- rnorm(length(B_seq), mean=B_seq, sd = B_stdev)
x_gc_rand <- 1/(1 + lambda_rand * exp(- B_seq_rand))





pi<-x_gc_rand
piSW<-x_gc_rand
piWS<-x_gc_rand
piWW <- x_gc_rand
piSS <- x_gc_rand

#This parameter represents the ratio of N-->N to W-->S diversity
a=1

#Standardization function (same as above, just a factored out)
pi_0_rand <- (2*lambda_rand)/(1+lambda_rand) + a



for(i in 1:length(x_gc_rand)){
  pi[i]<-((x_gc_rand[i])*lambda_rand[i]*2*(1/(1-exp(B_seq_rand[i]))+1/(B_seq_rand[i])) +
            (1-(x_gc_rand[i]))*2*(1/(1-exp(-B_seq_rand[i]))-1/(B_seq_rand[i])) + 1)
  
  piSW[i]<-lambda_rand[i]*2*(1/(1-exp(B_seq_rand[i]))+1/(B_seq_rand[i]))
  piWS[i]<-2*(1/(1-exp(-B_seq_rand[i]))-1/(B_seq_rand[i]))
  piWW[i]<-rep(0.5, times=801)
  piSS[i]<-rep(0.5, times=801)
  
  pi[i] <- pi[i]/pi_0_rand[i]
  piSW[i]<-piSW[i]/pi_0_rand[i]
  piWS[i]<-piWS[i]/pi_0_rand[i]
  piWW[i]<-piWW[i]/pi_0_rand[i]
  piSS[i]<-piSS[i]/pi_0_rand[i]
}

#Plot for equilibrium GC content and π-rel when B and lambda have been sampled from a normal distribution
plot(x_gc_rand, pi,xlim=c(0,1), col=col[1], pch=16, ylim =c(0,1.5), lwd=3,  cex.lab=1.3,  xlab = "Equilibrium GC content", ylab = expression(pi~ ""["rel"]), main = "gBGC-mutation-drift model") 
points(x_gc_rand, piSW, col=col3[2], lwd=3, pch=16)
points(x_gc_rand, piWS, col=col3[5], lwd=3, pch=17)
points(x_gc_rand, piWW, col=col3[4], lwd=3, lty=1, pch=3)
points(x_gc_rand, piSS, col=col3[1], lwd=3, lty=1, pch=4)
abline(h=1,col="grey",lty=2)
legend("topright", legend =c("All","SS", "SW", "WS", "WW"), col=c(col[1], col3[1], col3[2], col3[5], col3[4]), pch=c(NA,4,16,17,3))

