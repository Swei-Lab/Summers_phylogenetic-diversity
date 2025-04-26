#Density of Infected Nymphs###

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
View(all)

#data cleaning####
all2<-all %>%
  mutate(Site=as.factor(Site)) %>%
  mutate(Rodent_SR=as.numeric(Rodent_SR))%>%
  mutate(Camera_SR=as.numeric(Camera_SR))%>%
  mutate(SR=as.numeric(SR))%>%
  mutate(lagged_DON=as.numeric(lagged_DON)) %>%
  mutate(Lagged_Pos=as.numeric(Lagged_Pos))
str(all2)

hist(all2$Lagged_Pos)

all2<-na.exclude(all2)

#Best traditonal metric####
mod1<-glmmTMB(Lagged_Pos~SR+(1|Site),family = nbinom2, 
                ziformula=~1, data=all2)
summary(mod1)
#168

mod2<-glmmTMB(Lagged_Pos~Rodent_SD+Rodent_SR+Camera_SD+Camera_SR+
                (1|Site),family = nbinom2, na.action = na.exclude,
            ziformula=~1, data=all2)
summary(mod2)
DREDGEdotap <- dredge(mod2, rank=AICc, trace=2)
#kept only rodent sr
r.squaredGLMM(mod2)
#171
#p=0.065

sim<-simulateResiduals(fittedModel=mod22, family=nbinom2)
plot(sim)

mod22<-glmmTMB(Lagged_Pos~Rodent_SR+
                (1|Site),family = nbinom2,na.action=na.exclude, 
              ziformula=~1, data=all2)
summary(mod22)
#167.4
r.squaredGLMM(mod22)

#Phylogenetic metrics

mod3<-glmmTMB(Lagged_Pos~PD+(1|Site),family = nbinom2, 
                ziformula=~1, data=all2)
summary(mod3)
#166
r.squaredGLMM(mod3) #0.7

mod4<-glmmTMB(Lagged_Pos~mpd_taxa+(1|Site),family = nbinom2,
         na.action=na.exclude,
                ziformula=~1, data=all2)
summary(mod4)
#163, significant 0.0386
r.squaredGLMM(mod4)
#0.77

mod5<-glmmTMB(Lagged_Pos~mntd_taxa+(1|Site),family = nbinom2, 
              ziformula=~1, data=all2)
summary(mod5)
#165.8
r.squaredGLMM(mod5)
#0.71

#Best combination####
#shannon and pd
mod6<-glmmTMB(Lagged_Pos~PD + Rodent_SD+Camera_SD+
                 (1|Site),ziformula=~1,data=all2, family = nbinom2())
summary(mod6)
DREDGEdotap <- dredge(mod6, rank=AICc, trace=2)
#none
r.squaredGLMM(mod6)

#mpd
mod7<-glmmTMB(Lagged_Pos~Rodent_SR+Rodent_SD+Camera_SR+Camera_SD+
                   mpd_taxa+ (1|Site),data=all2, 
                 family = nbinom2(),
                 ziformula=~1)
summary(mod7)
DREDGEdotap <- dredge(mod7, rank=AICc, trace=2)
#kept mpd only

#mntd
mod8<-glmmTMB(Lagged_Pos~Rodent_SR+Rodent_SD+Camera_SR+Camera_SD+
                mntd_taxa+ (1|Site),data=all2, 
              family = nbinom2(),
              ziformula=~1)
DREDGEdotap <- dredge(mod8, rank=AICc, trace=2)
#kept rodent sr only but not sign

#best model#
model<-glmmTMB(Lagged_Pos~mpd_taxa+(1|Site),family = nbinom2, 
               ziformula=~1, data=all2)
summary(model)
exp(-0.89051) #back transform
r.squaredGLMM(model)
#163, significant 0.0386
#0.77 fit

#graph####
model<-glmmTMB(Lagged_Pos~mpd_taxa+(1|Site),family = nbinom2, 
                     ziformula=~1, data=all2)

mpdplot <- ggpredict(model, c("mpd_taxa"))

mpd<-plot(mpdplot)

Mplot<-mpd + theme_test() +geom_point(data = all2, aes(x = mpd_taxa, y = Lagged_Pos), 
                                     alpha = 0.5, colour = "black", shape = "circle") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "lightyellow", alpha = 0.5)+
  geom_line(color="orange")+
  xlab("Mean pairwise distance") +
  ylab("Density of Infected Nymphs (Bbss)")
Mplot

#lizards####

mod8<-glmmTMB(Lagged_Pos~mpd_taxa+Lizard+(1|Site),data=all2, 
              family = nbinom2(),na.action = na.exclude,
              ziformula=~1)
summary(mod8)
r.squaredGLMM(mod8)

mod8<-glmmTMB(Lagged_Pos~Lizard+(1|Site)+(1|mpd_taxa),data=all2, 
              family = nbinom2(),na.action = na.exclude,
              ziformula=~1)
summary(mod8)
