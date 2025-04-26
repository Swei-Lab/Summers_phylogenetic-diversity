#Nymphal Infection Prevalence models####

library(MuMIn)
install.packages("MuMIn")
library(aod)
library(dplyr)
library(vcd)
library(glmmTMB)
library(ggeffects)
library(ggplot2)
library(DHARMa)
library(car)


setwd("C:/Users/shann/OneDrive/Documents/Swei lab/R")

all<-read.csv("all_final.csv",header=TRUE)

#data cleaning####
all2<-all %>%
  mutate(Site=as.factor(Site)) %>%
  mutate(Rodent_SR=as.numeric(Rodent_SR))%>%
  mutate(Camera_SR=as.numeric(Camera_SR))%>%
  mutate(SR=as.numeric(SR))%>%
  mutate(lagged_DON=as.numeric(lagged_DON)) %>%
  mutate(Lagged_Pos=as.numeric(Lagged_Pos)) %>%
  mutate(Lagged_test=as.numeric(Lagged_test))
str(all2)

all2<-na.exclude(all2)
View(all2)

#Traditional metrics#####
#test of total species richness and shannon diversity
model1<-glmmTMB(Lagged_Pos/Lagged_test~SR+Rodent_SD+Camera_SD+(1|Site),
                data=all2, family = binomial(link="logit"),
                weights = Lagged_test)
summary(model1)
#aic 142
#not good r squared

options(na.action="na.fail")

DREDGEdotap <- dredge(model1, rank=AICc, trace=2)
#Dropped SR
r.squaredGLMM(model1)

##Shannon and richness
model2<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SR+Rodent_SD+Camera_SR+Camera_SD
                +(1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model2)

sim<-simulateResiduals(fittedModel=model2.2, family=binomial)
plot(sim)

DREDGEdotap <- dredge(model2, rank=AICc, trace=2)
#only keeps rodent SD
r.squaredGLMM(model2)

model2.2<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+
                    (1|Site),data=all2, 
                  family = binomial(link="logit"), weights = Lagged_test)
summary(model2.2)
#141 aic
r.squaredGLMM(model2.2)

#dropping randoms
model2<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SR+Rodent_SD+Camera_SR+Camera_SD+
                  (1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model2)
DREDGEdotap <- dredge(model2, rank=AICc, trace=2)
#keeps rodent SD and camera sr w/ site
#keeps camera sd and rodent sd w/ year
r.squaredGLMM(model2)


model22<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SR+
                   (1|Site),data=all2,na.action = na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model22)
#aic 143,0.046, 
r.squaredGLMM(model22)
#.20
exp(0.6480) #backtransfrom
exp(-0.1759)

sim<-simulateResiduals(fittedModel=model22, family=binomial)
plot(sim)

plot(all2$Rodent_SD,all2$Lagged_Bbss)
plot(all2$Camera_SR,all2$Lagged_Bbss)

model25<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+
                   (1|Site),data=all2, na.action = na.exclude,
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model25)
r.squaredGLMM(model25)

model23<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SD+
                   (1|Site),data=all2, na.action = na.exclude,
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model23)
#rodent 0.029, 146.3
r.squaredGLMM(model23)
#aic 149, sign but i think year is driving that

#Phylogenetic Metrics####

model3<-glmmTMB(Lagged_Pos/Lagged_test~PD+
                  (1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model3)
r.squaredGLMM(model3)
#aic lower without year
#142.5 

model<-lm(Lagged_Bbss~PD+mntd_taxa, data=all2)
vif(model)
#do all seperately

model4<-glmmTMB(Lagged_Pos/Lagged_test~mpd_taxa+
                  (1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model4)
r.squaredGLMM(model4)
#144.2

model5<-glmmTMB(Lagged_Pos/Lagged_test~mntd_taxa+
                  (1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model5)
r.squaredGLMM(model5)
#144.8

#Best combination####
#PD
model<-lm(Lagged_Bbss~PD+Rodent_SR+Camera_SR, data=all2)
vif(model) #no sr w/ PD

model6<-glmmTMB(Lagged_Pos/Lagged_test~PD+Rodent_SD+Camera_SD+(1|Site)
                ,data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model6)
r.squaredGLMM(model6)
DREDGEdotap <- dredge(model6, rank=AICc, trace=2)
#year drops PD, better r2 but higher AIC

model66<-glmmTMB(Lagged_Pos/Lagged_test~PD+Rodent_SD+(1|Site)
                 ,data=all2, na.action=na.exclude,
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model66)
#aic 143.9, 0.0248
r.squaredGLMM(model66) #0.19
DREDGEdotap <- dredge(model66, rank=AICc, trace=2)
#backtransform estimates
exp(-0.002544)
exp(0.735176)

sim<-simulateResiduals(fittedModel=model66, family=binomial)
plot(sim)


#mpd
model<-lm(Lagged_Bbss~mpd_taxa+Rodent_SD+Camera_SR, data=all2)
vif(model)

model7<-glmmTMB(Lagged_Pos/Lagged_test~mpd_taxa+Rodent_SD+Rodent_SR+Camera_SR
                +Camera_SD+(1|Year),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model7)
r.squaredGLMM(model7)
DREDGEdotap <- dredge(model7, rank=AICc, trace=2)
#with site only keep rodent sd and camera sr
#with year keep cam, rodent SD and mpd

model72<-glmmTMB(Lagged_Pos/Lagged_test~mpd_taxa+Rodent_SD
                 +Camera_SR+(1|Site),data=all2,na.action=na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model72)
#148.3
r.squaredGLMM(model72) #0.26

model73<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD
                 +mpd_taxa+(1|Site),data=all2, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model73)
#139.4
r.squaredGLMM(model73) #0.22

#mntd

model8<-glmmTMB(Lagged_Pos/Lagged_test~mntd_taxa+Rodent_SD+Rodent_SR+Camera_SR
                +Camera_SD+(1|Year),data=all2, na.action=na.exclude,
                family = binomial(link="logit"), weights = Lagged_test)
summary(model8)
r.squaredGLMM(model8)
DREDGEdotap <- dredge(model8, rank=AICc, trace=2)
#drops mntd for both year and site as random

model8<-glmmTMB(Lagged_Pos/Lagged_test~mntd_taxa+Rodent_SD
                +Camera_SR+(1|Site),data=all2, na.action=na.exclude,
                family = binomial(link="logit"), weights = Lagged_test)
summary(model8)

###best model##
model66<-glmmTMB(Lagged_Pos/Lagged_test~PD+Rodent_SD+(1|Site)
                 ,data=all2, na.action = na.exclude,
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model66)
View(all2)

#graph####
model<-glmmTMB(Lagged_Bbss~Rodent_SD+(1|Site)
               ,data=all2, 
               family = binomial(link="logit"), weights = Lagged_test)
summary(model)
Rodentplot <- ggpredict(model, c("Rodent_SD [all]"))

rodent<-plot(Rodentplot)

Rplot<-rodent + theme_test() +geom_point(data = all2, aes(x = Rodent_SD, y = Lagged_Bbss), 
                                         alpha = 0.5,color="black", shape = "circle") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "skyblue", alpha = 0.5)+
  geom_line(color="darkblue")+
  xlab("Rodent Shannon diversity") +
  ylab("Nymphal infection prevalence(Bbss)")
Rplot


model3<-glmmTMB(Lagged_Bbss~PD+
                  (1|Site),data=all2, 
                family = binomial(link="logit"), weights = Lagged_test)

plot <- ggpredict(model3, c("PD"))

PD<-plot(plot)

PDplot<-PD + theme_test() +geom_point(data = all2, aes(x = PD, y = Lagged_Bbss), 
                                      alpha = 0.5, colour = "black", shape = "circle") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "orange", alpha = 0.5)+
  geom_line(color="darkred")+
  xlab("Faith's PD") +
  ylab("Nymphal infection prevalence(Bbss)")
PDplot

#interaction graph####
library(sjPlot)
library(sjmisc)
library(ggplot2)

mod<-glmmTMB(Lagged_Bbss~PD+Rodent_SD
             ,data=all2, family = binomial(link="logit"), weights = Lagged_test)
summary(mod)

plot_model(mod, type = "pred", terms = c("PD", "Rodent_SD"))

#Random effect testing####
model22<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SR+
                   (1|Site)+(1|mpd_taxa),data=all2,na.action = na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model22)
r.squaredGLMM(model22)

sim<-simulateResiduals(fittedModel=model22, family=binomial)
plot(sim)

model23<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SR+
                   (1|Site)+(1|mntd_taxa),data=all2,na.action = na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model23)
r.squaredGLMM(model23)

##Lizards####
model22<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SR+Lizard+
                   (1|Site),data=all2,na.action = na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model22)
r.squaredGLMM(model22)

model2<-glmmTMB(Lagged_Pos/Lagged_test~Lizard+PD+
                  (1|Site),data=all2,na.action = na.exclude, 
                family = binomial(link="logit"), weights = Lagged_test)
summary(model2)


##best models, adjust pvalue####
model22<-glmmTMB(Lagged_Pos/Lagged_test~Rodent_SD+Camera_SR+
                   (1|Site),data=all2,na.action = na.exclude, 
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model22)
#aic 143,0.046, 
r.squaredGLMM(model22)

model66<-glmmTMB(Lagged_Pos/Lagged_test~PD+Rodent_SD+(1|Site)
                 ,data=all2, na.action=na.exclude,
                 family = binomial(link="logit"), weights = Lagged_test)
summary(model66)

p<-c(0.0248, 0.046510)
p.adjust(p, method = "holm", n = length(p))
p.adjust(p, method = "BH", n = length(p))



