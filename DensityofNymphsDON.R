#Density of nymphs###

library(lme4)
library(MuMIn)
library(aod)
library(dplyr)
library(vcd)
library(glmmTMB)
library(ggeffects)
library(ggplot2)
library(DHARMa)

all<-read.csv("all_final.csv",header=TRUE)
View(all2)

#data cleaning####
all2<-all %>%
  mutate(Site=as.factor(Site)) %>%
  mutate(Rodent_SR=as.numeric(Rodent_SR))%>%
  mutate(Camera_SR=as.numeric(Camera_SR))%>%
  mutate(SR=as.numeric(SR))%>%
  mutate(lagged_DON=as.numeric(lagged_DON)) %>%
  mutate(Lagged_Pos=as.numeric(Lagged_Pos))%>%
  mutate(Lizard=as.factor(Lizard))
str(all2)

hist(all2$Lagged_DON)

all2<-na.exclude(all2)

#Traditional metrics####
DON_SR<-glmmTMB(lagged_DON~SR+(1|Site), 
                data=all2, family = nbinom1(),
                ziformula=~1)
summary(DON_SR)

mod1<-glmmTMB(lagged_DON~Camera_SR+Camera_SD+Rodent_SR+Rodent_SD+
                (1|Year)+(1|Site), 
              data=all2, family = nbinom1(),
              ziformula=~1)
summary(mod1)
DREDGEdotap <- dredge(mod1, rank=AICc, trace=2)
#only keeps rodent sr
r.squaredGLMM(mod1)
#poor

mod11<-glmmTMB(lagged_DON~Rodent_SR+
                (1|Site), 
              data=all2, family = nbinom1(),
              ziformula=~1)
summary(mod11)
#278.9, significant
#backtransform estimate
exp(0.40637)
r.squaredGLMM(mod11) #0.92

#Phylogenetic metrics####
mod2<-glmmTMB(lagged_DON~mpd_taxa+(1|Site),family = nbinom1, 
                    ziformula=~1, data=all2)
summary(mod2)
#282, significant
#backtransform estimate
exp(-0.5571)
r.squaredGLMM(mod2) #0.97

mod3<-glmmTMB(lagged_DON~mntd_taxa+(1|Site),family = nbinom1,na.action=na.exclude, 
              ziformula=~1, data=all2)
summary(mod3)
#285.7
r.squaredGLMM(mod3) #0.89

mod4<-glmmTMB(lagged_DON~PD+(1|Site),family = nbinom1, 
              ziformula=~1, data=all2)
summary(mod4)
#287.1
r.squaredGLMM(mod4) #0.96

#Best combination####

#mpd
#best model
mod55<-glmmTMB(lagged_DON~Rodent_SR+mpd_taxa+
                (1|Site),data=all2, family = nbinom1(),
              ziformula=~1)
summary(mod55)
#317.8, rodent very significant
DREDGEdotap <- dredge(mod55, rank=AICc, trace=2)
r.squaredGLMM(mod55) #92.8

#mntd
mod66<-glmmTMB(lagged_DON~Rodent_SR+mntd_taxa+ 
                 (1|Site),data=all2, family = nbinom1(),
               ziformula=~1)
summary(mod66)
DREDGEdotap <- dredge(mod6, rank=AICc, trace=2)
#only rodent sr
r.squaredGLMM(mod66)

#PD can't be combined with SR so tried with SD
mod7<-glmmTMB(lagged_DON~PD+Rodent_SD+Camera_SD+
                    (1|Year)+(1|Site),data=all2, family = nbinom1())
summary(mod7)

DREDGEdotap <- dredge(mod7, rank=AICc, trace=2)
#none
r.squaredGLMM(mod7)

##Lizards####
mod5<-glmmTMB(lagged_DON~Lizard+
                (1|Site),data=all2,na.action = na.exclude, family = nbinom1(),
              ziformula=~1)
summary(mod5) #316
r.squaredGLMM(mod5) #93
#worth putting in?

lizm<-glmmTMB(lagged_DON~Rodent_SR+Lizard+mpd_taxa+
                    (1|Site),data=all2,na.action = na.exclude, family = nbinom1(),
                  ziformula=~1)
summary(lizm) #310
#backtransform estimates
exp(0.19407)
exp(0.63774)
exp(-0.22859)
r.squaredGLMM(lizm) #95

#lizard graph
plot(all2$lagged_DON, all2$Lizard)
Lizard<-as.factor(all2$Lizard)
cdplot(all2$Lizard, all2$lagged_DON)


#best models, corrected p####
mod44<-glmmTMB(lagged_DON~Rodent_SR+Lizard+mpd_taxa+
                 (1|Site),data=all2,na.action = na.exclude, family = nbinom1(),
               ziformula=~1)
summary(mod44)
#310
r.squaredGLMM(mod44) #95.7

mod11<-glmmTMB(lagged_DON~Rodent_SR+
                 (1|Site), 
               data=all2, family = nbinom1(),
               ziformula=~1)
summary(mod11)

mod2<-glmmTMB(lagged_DON~mpd_taxa+(1|Site),family = nbinom1, 
              ziformula=~1, data=all2)
summary(mod2)

p<-c(0.031602,0.142010,4.06e-07,0.00592)
p.adjust(p, method = "holm", n = length(p))
1.6240e-06

#Graph#####
model<-glmmTMB(lagged_DON~Rodent_SR+
                 (1|Site), 
               data=all2, family = nbinom1(),
               ziformula=~1)
summary(model)

rplot <- ggpredict(model, c("Rodent_SR"))

rodent<-plot(rplot)

Rplot<-rodent + theme_test() +geom_point(data = all2, aes(x = Rodent_SR, y = lagged_DON), 
                                         alpha = 0.5, colour = "black", shape = "circle") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "lightgreen", alpha = 0.5)+
  geom_line(color="darkgreen")+
  xlab("Rodent Species Richness") +
  ylab("Density of Nymphs") +
  ggtitle("Relationship between DON and Rodent SR")
Rplot

mod<-glmmTMB(lagged_DON~mpd_taxa+ 
               (1|Site)+(1|Year),data=all2, family = nbinom1(),
             ziformula=~1)
mplot <- ggpredict(mod, c("mpd_taxa"))

mpd<-plot(mplot)

MPDplot<-mpd + theme_test() +geom_point(data = all2, aes(x = mpd_taxa, y = lagged_DON), 
                                        alpha = 0.5, colour = "black", shape = "circle") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "violet", alpha = 0.5)+
  geom_line(color="purple")+
  xlab("Mean Pairwise distance") +
  ylab("Density of Nymphs") +
  ggtitle("Relationship between DON and Rodent SR")
MPDplot

