
rm(list = ls())
library(openxlsx)
library(plyr)
library(tidyverse)
library(lme4)
library(ggpubr) 

FigDir = '../Figs/'
pd <- position_dodge(0.6) # set dodge width
# some common ggplot settings
common = theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(colour = "black"), # axis color
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5), # center title
        text=element_text(size=16)) # text size & font

use.col.conds = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")
use.col.sess = c("1" = "pink","2" = "red","3" = "darkred")
use.col.conds.day1 = c("sham" = "gray","cTBS" = "green")
use.col.conds.day2 = c("sham" = "gray","cTBS" = "blue")

###############################
#### analyze conditioning data ####
###############################

conditioning_dat = read.xlsx('../ProcessedData/Conditioning.xlsx') %>%
  mutate(SelectedCue = if_else(SelectedCue == 0, NA_integer_, SelectedCue)) %>% # replace zero with NA
  mutate(Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('Sweet','Savory')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                        to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                             'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('C-S-S','C-S-S','S-C-S','S-S-C','S-C-S','S-S-C')),
         StimOrder_day2 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('S-C-S','S-S-C','C-S-S','C-S-S','S-S-C','S-C-S'))) %>%
  mutate(SubID=factor(SubID),
         StimOrder=factor(StimOrder),
         StimOrder_day1=factor(StimOrder_day1),
         StimOrder_day2=factor(StimOrder_day2),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Set=factor(Set),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         ChoiceType=factor(ChoiceType),
         StimLoc=factor(StimLoc))

save(conditioning_dat,file = '../ProcessedData/Conditioning.RData')

# does conditioning improve over runs?
# collapsing across 3 conditions/sessions
p1=conditioning_dat %>%
  group_by(SubID,Run) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=SubID)) +
  geom_line(show.legend = F,alpha=0.5) +
  geom_line(stat = "summary", fun = "mean", linewidth = 2, color = "black") +
  labs(x = "Runs", y = "P(Choosing Odors)",
       title = "Conditioning (all subjects)") +
  common

model_conditioning_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_0 <- glmer(OdorChosen ~ (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_1,model_conditioning_0)
summary(model_conditioning_1)
# yes conditioning improves over runs



# does different TMS expt condition affect learning?
p3=conditioning_dat %>%
  group_by(SubID,Run,Cond) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common

# same as p3 but using error bar
p3a=conditioning_dat %>%
  group_by(SubID,Run,Cond) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "errorbar", width = 0.2) +  
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common

model_conditioning_2 <- glmer(OdorChosen ~ Run + Cond + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_3 <- glmer(OdorChosen ~ Run * Cond + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_2,model_conditioning_1)
anova(model_conditioning_2,model_conditioning_3)
summary(model_conditioning_3)
# interaction effect between cTBS-sham/sham-cTBS conditions and runs


# is this modulated by different TMS locations?
p4=conditioning_dat %>%
  group_by(SubID,Run,Cond,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) + 
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common +
  theme(legend.position = c(0.85,0.3))

model_conditioning_4 <- glmer(OdorChosen ~ Run * Cond + StimLoc + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_4,model_conditioning_3) # need StimLoc? p=0.056
model_conditioning_5 <- glmer(OdorChosen ~ Run * Cond * StimLoc + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_4,model_conditioning_5)
anova(model_conditioning_3,model_conditioning_5)
summary(model_conditioning_5)

# any difference from savory vs sweet odors?
conditioning_dat %>%
  group_by(SubID,Run,Cond,ChoiceType) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  facet_wrap(~ChoiceType) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +
  labs(x = "Runs", y = "P(Choosing Odors)", 
       color = "", fill = "") +
  common

model_conditioning_6 <- glmer(OdorChosen ~ Run * Cond + ChoiceType + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_6,model_conditioning_3)
# lack of ChoiceType main effect
# no need to add ChoiceType

# any difference from two sets?
conditioning_dat %>%
  group_by(SubID,Run,Cond,Set) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond)) +
  facet_wrap(~Set) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       color = "", fill = "") +
  common

model_conditioning_7 <- glmer(OdorChosen ~ Run * Cond + Set + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_7,model_conditioning_3)
# lack of Set main effect
# no need to add Set

# put Set and StimLoc factorially
conditioning_dat %>%
  group_by(SubID,Run,Cond,Set,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond)) +
  facet_wrap(~Set + StimLoc) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  stat_summary(fun = mean, geom = "line", size = 1, na.rm = T) +  # Average line
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "ribbon", alpha = 0.3) +
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common

####################################
# focus on cond_day1: cTBS vs sham
####################################

model_conditioning1_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_2 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_3 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_2,model_conditioning1_3)
anova(model_conditioning1_1,model_conditioning1_2)
summary(model_conditioning1_3) # interactions b/t runs and Day 1 conditions

# conditioning performance across runs, by diff day1 TMS
pp1=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", linewidth = 0, show.legend = F,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "Choosing Odors", color = "Day1 TMS") +
  common

# conditioning performance across runs, by diff day1 TMS & stim locs
pp2=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  common +
  theme(legend.position = c(0.85,0.3))

# same to pp2, but using error bar
pp2a=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "errorbar", linewidth = 1, width = 0.3) +  
  labs(x = "Runs", y = "Choosing Odors", color = "Day1 TMS") +
  common +
  theme(legend.position = c(0.85,0.3))

model_conditioning1_4 <- glmer(OdorChosen ~ Run * Cond_day1 + StimLoc + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_4,model_conditioning1_3) # need StimLoc? p=0.055
model_conditioning1_5 <- glmer(OdorChosen ~ Run * Cond_day1 + Run * StimLoc + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_4,model_conditioning1_5)
anova(model_conditioning1_3,model_conditioning1_5)
summary(model_conditioning1_5)
# from the models, it seems to suggest an interaction b/t runs & conditions
# and an interaction b.t runs & stim locations

# aOFC only
model_conditioning_aOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Anterior'),
                                   family = 'binomial')
model_conditioning_aOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Anterior'),
                                   family = 'binomial')
model_conditioning_aOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Anterior'),
                                   family = 'binomial')
anova(model_conditioning_aOFC_0,model_conditioning_aOFC_1) # day1 TMS main effect
anova(model_conditioning_aOFC_1,model_conditioning_aOFC_2) # interaction b/t TMS & runs
summary(model_conditioning_aOFC_2)
# describe: faster learning with runs under sham for aOFC stimulation group

# pOFC only
model_conditioning_pOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Posterior'),
                                   family = 'binomial')
model_conditioning_pOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Posterior'),
                                   family = 'binomial')
model_conditioning_pOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = subset(conditioning_dat,StimLoc=='Posterior'),
                                   family = 'binomial')
anova(model_conditioning_pOFC_0,model_conditioning_pOFC_1) # marginal day1 TMS main effect
anova(model_conditioning_pOFC_1,model_conditioning_pOFC_2) # no interaction b/t TMS & runs
summary(model_conditioning_pOFC_1)

# check conditioning performance with Day 2 TMS condition (for fun!)
conditioning_dat %>%
  group_by(SubID,Run,Cond_day2,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day2, fill = Cond_day2)) +
  scale_color_manual(values = use.col.conds.day2) +
  scale_fill_manual(values = use.col.conds.day2) +
  facet_wrap(~StimLoc) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "Choosing Odors", 
       color = "Condition-Day2",fill = "Condition-Day2") +
  common

# look at Sess, Cond, runs together
# note that subjects in each panel are not the same 
conditioning_dat %>%
  group_by(SubID,Run,Cond,Sess,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, 
             color = factor(Sess), linetype = factor(Sess),
             fill = factor(Sess))) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  facet_wrap(~StimLoc + Cond) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       color = "Session", fill = "Session", linetype = "Session") +
  common

# look at Cond_day1 instead of conds
conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,Sess,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, linetype = factor(Sess),
             color = factor(Sess), fill = factor(Sess))) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  facet_wrap(~StimLoc + Cond_day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       title = "Conditioning by stim loc & Day1 TMS",
       color = "Session", fill = "Session", linetype = "Session") +
  common

pp3=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,Sess,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, fill = Cond_day1,
             color = Cond_day1)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  facet_wrap(~StimLoc + factor(Sess)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       title = "Conditioning by stim loc & Session #",
       color = "", fill = "") +
  common

# overall change of conditioning across sessions 1,2,3
conditioning_dat %>%
  group_by(SubID,Run,Sess) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess),fill=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(x = "Runs", y = "P(Choosing Odors)",
       color = 'Session', title = 'Conditioning over runs & sessions') +
  common

pp4_sub_group=conditioning_dat %>%
  group_by(SubID,Run,Sess,StimOrder_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,linetype=factor(Sess),
             color=factor(Sess),fill=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +
  facet_wrap(~StimLoc+StimOrder_day1,nrow = 2) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  labs(x = "Runs", y = "P(Choosing Odors)", fill = 'Session',
       color = 'Session', linetype = 'Session',
       title = 'Conditioning by stim loc & Day 1 stim order') +
  common

pdf(file.path(FigDir,'Conditioning.pdf'),10,8)
print(p1)
print(p3)
print(p4)
print(pp1)
print(pp2)
print(pp3)
print(pp4_sub_group)
dev.off()

pdf(file.path(FigDir,'Fig2a.pdf'),8,4)
print(pp2)
dev.off()

pdf(file.path(FigDir,'Fig2b.pdf'),8,6)
print(pp3)
dev.off()


# check RT btw
prt1 = conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc,Sess) %>%
  dplyr::summarise(meanRT=mean(rt,na.rm = T)) %>%
  ggplot(aes(x = Run, y = meanRT, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "RT (ms)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  common

prt2 = conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc,Sess) %>%
  dplyr::summarise(meanRT=mean(rt,na.rm = T)) %>%
  ggplot(aes(x = Run, y = meanRT, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc + factor(Sess)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "RT (ms)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  common

pdf(file.path(FigDir,'Conditioning_rt.pdf'),10,8)
print(prt1)
print(prt2)
dev.off()

###############################
# analyze pleasantness ratings data
###############################

dat = read.xlsx('../ProcessedData/Pleasant.xlsx') %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre','Post')),
         IfDevalue=mapvalues(IfDevalue,from=c('0','1'),to=c('Non-Dv','Devalued')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre','Post')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         IfDevalue=factor(IfDevalue,levels=c('Non-Dv','Devalued')),
         Cond=factor(Cond)) %>%
  arrange(SubID,Cond,IfDevalue)
  
dat = dat %>% 
  mutate(pair=rep(c(1:(nrow(dat)/2)),each=2))

p1=ggplot(dat,aes(x=IfDevalue,y=Pleasant,fill=PrePost)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(color=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
              size = 1, alpha = 0.8) +
  facet_wrap(~Cond) +
  labs(x = "Odors", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

p2=ggplot(dat,aes(x=IfDevalue,y=Pleasant,fill=PrePost)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(color=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
              size = 1, alpha = 0.8) +
  facet_wrap(~Sess) +
  labs(x = "Odors", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

p3=ggplot(dat,aes(x=IfDevalue,y=Pleasant,fill=PrePost)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(color=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
              size = 1, alpha = 0.8) +
  labs(x = "Odors", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

pdf(file.path(FigDir,'PleasantRatings.pdf'),8,4)
print(p1)
print(p2)
print(p3)
dev.off()


# Fit the mixed-effects model with two factors and their interaction
model1 <- lmer(Pleasant ~ IfDevalue * PrePost + (1|SubID), data = dat)
model2 <- lmer(Pleasant ~ IfDevalue + PrePost + (1|SubID), data = dat)

summary(model1)
summary(model2)
anova(model1,model2)
# collapsing across all 3 sessions, there was an interaction between devalue
# $ test times

model3 <- lmer(Pleasant ~ IfDevalue * PrePost + Sess + (1|SubID), data = dat)
model4 <- lmer(Pleasant ~ IfDevalue * PrePost * Sess + (1|SubID), data = dat)
summary(model3)
anova(model3,model1)
anova(model3, model4)
# this interaction effect does not change across sessions

model3 <- lmer(Pleasant ~ IfDevalue * PrePost + Cond + (1|SubID), data = dat)
model4 <- lmer(Pleasant ~ IfDevalue * PrePost * Cond + (1|SubID), data = dat)
anova(model3, model4)
# this interaction effect does not change across task conditions

# a diff way to plot to show lines pre & post meal
p=ggplot(dat,aes(x=PrePost,y=Pleasant,fill=PrePost)) +
  geom_boxplot(width = 0.6, alpha = 0.4, outlier.alpha = 0) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=PrePost,group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  facet_wrap(~IfDevalue+Cond) +
  labs(x = "", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        legend.position = 'none')

pdf(file.path(FigDir,'PleasantRatings_w_lines.pdf'),8,8)
print(p)
dev.off()

# to plot the change of pleasantness ratings (post meal-pre meal)
# create a new df 'reduced' from 'dat'
pre_plea = subset(dat,PrePost=='Pre')$Pleasant
post_plea = subset(dat,PrePost=='Post')$Pleasant
change_plea = post_plea - pre_plea
reduced = subset(dat,PrePost=='Pre')
reduced$PrePost = 'Change'
reduced$Pleasant = change_plea

reduced = reduced %>% 
  mutate(pair=rep(c(1:(nrow(reduced)/2)),each=2))

q1=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~Cond) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q1a=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc + Cond) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q2=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~Sess) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q3=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.2, color = 'darkgray') +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, position = position_dodge(0.75),
               outlier.alpha = 0) +
  labs(x = "", y = "Pleasantness Ratings (Post-Pre)") + common +
  theme(strip.background = element_blank())

pdf(file.path(FigDir,'PleasantChanges.pdf'),6,4)
print(q1)
print(q2)
print(q1a)
dev.off()

pdf(file.path(FigDir,'PleasantChanges_all_sessions.pdf'),4,5)
print(q3)
dev.off()

# Fit the mixed-effects model
model1 <- lmer(Pleasant ~ IfDevalue + (1|SubID), data = reduced)
model0 <- lmer(Pleasant ~ (1|SubID), data = reduced)

summary(model1)
anova(model1)
summary(model0)
anova(model1,model0)
# collapsing across all 3 sessions, there was a devaluation effect

model2 <- lmer(Pleasant ~ Cond + IfDevalue + (1|SubID), data = reduced)
model3 <- lmer(Pleasant ~ Cond * IfDevalue + (1|SubID), data = reduced)

anova(model2,model1)
summary(model2)
summary(model3)
anova(model2,model3)
anova(model3,model1)
# the devaluation effect was not modulated by conditions


# further create a measure of devaluation effect per session per subject
# change of pleasantness rating (devalued - non-devalued)
# post-pre: more negative value means better devaluation effect

ND = subset(reduced,IfDevalue=='Non-Dv')$Pleasant
D = subset(reduced,IfDevalue=='Devalued')$Pleasant
Didx = D - ND
rereduced = subset(reduced,IfDevalue=='Non-Dv')
rereduced$PrePost = NULL # these cols are meanless now after reducing
rereduced$Pleasant = NULL
rereduced$IfDevalue = NULL
rereduced$OdorType = NULL
rereduced$pair = NULL
rereduced$Didx = Didx

# Compare the mixed-effects model
model0 <- lmer(Didx ~ (1|SubID), data = rereduced)
model1 <- lmer(Didx ~ Cond + (1|SubID), data = rereduced)
model2 <- lmer(Didx ~ Sess + (1|SubID), data = rereduced)
model3 <- lmer(Didx ~ Devalued + (1|SubID), data = rereduced)
model4 <- lmer(Didx ~ StimLoc + (1|SubID), data = rereduced)
summary(model0)
summary(model1)
anova(model0,model1)
anova(model0,model2)
anova(model0,model3)
anova(model0,model4)
# devaluation effect was not modulated by Cond, Sess, DevaluedOdor, StimLoc

p1=ggplot(rereduced,aes(x=Cond,y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Devaluation effect") + common +
  theme(strip.background = element_blank())

p2=ggplot(rereduced,aes(x=factor(Sess),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Devaluation effect") + common +
  theme(strip.background = element_blank())

p3=ggplot(rereduced,aes(x=factor(Devalued),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Devaluation effect") + common +
  theme(strip.background = element_blank())

p4=ggplot(rereduced,aes(x=factor(StimLoc),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Devaluation effect") + common +
  theme(strip.background = element_blank())

pdf(file.path(FigDir,'DevaluationEffect.pdf'),4.5,4)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

###############################
# analyze choice data
###############################

choice_dat = read.xlsx('../ProcessedData/Choices.xlsx') %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre','Post')),
         Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre','Post')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Cond=factor(Cond),
         Set=factor(Set),
         StimLoc=factor(StimLoc))

# only look at choice with one sweet and one savory (denoted as '_ss')
# perhaps where the effect is the strongest
choice_dat_ss = choice_dat %>%
  subset(ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) # Choice==1 means selecting sated odor

summary_choice_ss = choice_dat_ss %>%
  group_by(SubID,Cond,PrePost,Set,StimLoc) %>%
  dplyr::summarise(Choice=mean(Choice,na.rm = T))

# separted by stimulation location
p1=summary_choice_ss %>%
  ggplot(aes(x=Cond,y=Choice,fill=interaction(PrePost,Set))) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(aes(color=interaction(PrePost,Set)),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
              size = 1, alpha = 0.8) +
  facet_wrap(~StimLoc) +
  labs(x = "", y = "Choice of Sated Odor (Sweet vs. Savory)")

# collapsing aOFC and pOFC stimulations
p2=summary_choice_ss %>%
  ggplot(aes(x=Cond,y=Choice,fill=interaction(PrePost,Set))) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(aes(color=interaction(PrePost,Set)),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
              size = 1, alpha = 0.8) +
  labs(x = "", y = "Choice of Sated Odor (Sweet vs. Savory)")

pdf(file.path(FigDir,'ChoiceSatedOdor.pdf'),8,5)
print(p1)
print(p2)
dev.off()

# focus on Set B
# only on pOFC
use.dat = subset(choice_dat_ss, Set=='B' & StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + (1|SubID), data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ (1|SubID), data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)
# lower choice in sham-sham condition

# only on aOFC
use.dat = subset(choice_dat_ss, Set=='B' & StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + (1|SubID), data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ (1|SubID), data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# no sig difference b/t conditions

# if including both aOFC & pOFC
use.dat = subset(choice_dat_ss, Set=='B' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ Cond + StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond * StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3)
# interaction efffect between condition and location
# was significant

# Set A
# if only focused on pOFC
use.dat = subset(choice_dat_ss, Set=='A' &
                   StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond * PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# p = 0.090
# with trend but not sig

# if only focused on aOFC
use.dat = subset(choice_dat_ss, Set=='A' &
                   StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond * PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# no sig

# although not quite legit, want to see if stimloc modulates
# any interaction effect, but no
use.dat = subset(choice_dat_ss, Set=='A' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ Cond * PrePost * StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond * PrePost + StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3)

# improved from above by considering session-wise bias
# using pre(A) as baseline
# then get difference from pre to post
pre_choice = subset(summary_choice_ss,Set=='A' & PrePost=='Pre')$Choice
post_choice_A = subset(summary_choice_ss,Set=='A' & PrePost=='Post')$Choice
post_choice_B = subset(summary_choice_ss,Set=='B' & PrePost=='Post')$Choice
conds = subset(summary_choice_ss,Set=='A' & PrePost=='Pre')$Cond
# check if there is pre-post correlations of choices
df_pre_post_summary = data.frame(
  'conds'=conds,
  'pre_choice'=pre_choice,
  'post_choice_A'=post_choice_A,
  'post_choice_B'=post_choice_B
)

pcorr1=ggplot(df_pre_post_summary,aes(x=pre_choice,y=post_choice_A,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")
pcorr2=ggplot(df_pre_post_summary,aes(x=pre_choice,y=post_choice_B,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")
pcorr3=ggplot(df_pre_post_summary,aes(x=post_choice_A,y=post_choice_B,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")

avg_changeAB = (post_choice_A + post_choice_B)/2 - pre_choice
avg_changeB = post_choice_B - pre_choice
avg_changeA = post_choice_A - pre_choice

summary_choice_corrected = summary_choice_ss %>%
  ungroup() %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  select(SubID,StimLoc,Cond) %>%
  mutate(ChoiceChangeAB=avg_changeAB, # both A and B
         ChoiceChangeA=avg_changeA, # only A
         ChoiceChangeB=avg_changeB) # only B

c1=summary_choice_corrected %>%
  #subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  #geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre',
     y = "Choice of Sated Odor (Sweet vs. Savory)") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

c2=summary_choice_corrected %>%
  #subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  #geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(B only) - Pre',
       y = "Choice of Sated Odor (Sweet vs. Savory)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

c3=summary_choice_corrected %>%
  #subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  #geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A only) - Pre',
       y = "Choice of Sated Odor (Sweet vs. Savory)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

c4=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre',
       y = "Choice of Sated Odor (Sweet vs. Savory)") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

pdf(file.path(FigDir,'ChoiceSatedOdor_Changes.pdf'),6,5)
print(pcorr1)
print(pcorr2)
print(pcorr3)
print(c1)
print(c2)
print(c3)
print(c4)
dev.off()

# stat testing on pOFC
# test using Post(A & B) - Pre
# overall devaluation effect regardless of conditions: yes
# (later realizing this also contains cTBS-sham condition)
value_Y1Z1 = subset(summary_choice_corrected,StimLoc=='Posterior')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.013

value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                     StimLoc=='Posterior')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                   StimLoc=='Posterior')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? No
wilcox.test(value_Z1,alternative = 'l') # lower than zero? Yes, p=0.026
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s

# test using Post(B only) - Pre
value_Y2Z2 = subset(summary_choice_corrected,
                    StimLoc=='Posterior')$ChoiceChangeB
wilcox.test(value_Y2Z2,alternative = 'l') # lower than zero? No

# stat testing on aOFC
# test using Post(A & B) - Pre
# overall devaluation effect regardless of conditions: yes
value_Y1Z1 = subset(summary_choice_corrected,
                    StimLoc=='Anterior')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.010
value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                    StimLoc=='Anterior')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='Anterior')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? p = 0.065
wilcox.test(value_Z1,alternative = 'l') # lower than zero? n.s.
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s

# test using Post(B only) - Pre
value_Y2Z2 = subset(summary_choice_corrected,StimLoc=='Anterior')$ChoiceChangeB
wilcox.test(value_Y2Z2,alternative = 'l') # lower than zero? p=0.080
# overall effect was weaker so not further looking into separate conditions

# go back to choice_dat_ss to do the testing at trial level
# first calculate a session-wise pre-meal bias for each sub
pre_choice_df = summary_choice_ss %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  group_by(SubID,Cond) %>%
  dplyr::summarise(base=Choice)
# insert this base value to the whole choice dat ss df
choice_dat_ss_w_base = choice_dat_ss %>%
  mutate(base=NA)
for(i in 1:nrow(choice_dat_ss_w_base)){
  choice_dat_ss_w_base$base[i] = pre_choice_df$base[
    pre_choice_df$SubID==choice_dat_ss_w_base$SubID[i] &
    pre_choice_df$Cond==choice_dat_ss_w_base$Cond[i]]
}

# models are on the trial level
# with session-wise baseline odor preference as a covariate

# full model with all 3 conditions
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior')
model_choice_full <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
summary(model_choice_full)

# focus on the two conditions of interest
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)
# adding condition sig improved model fit
# and both variables are sig

# repeate this analysis on aOFC
# full model with all 3 conditions
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Anterior')
model_choice_full <- glmer(Choice ~ Cond + base + (1|SubID), 
                           data = use.dat,family = 'binomial')
summary(model_choice_full)

use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)

# put aOFC and pOFC together in the test
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ StimLoc * Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
summary(model_choice_1)

###############################
### correlate the choice updating with devaluation effect idx (plea) ###
###############################

# (having checked rows are aligned)
df_comb = cbind(summary_choice_corrected,
                rereduced[,c('Sess','Devalued','Didx')])
corr1=ggplot(df_comb,aes(x=Didx,y=ChoiceChangeAB)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "",
       x = 'Devaluation on Pleasantness rating',
       y = 'Choice update')
corr2=ggplot(df_comb,aes(x=Didx,y=ChoiceChangeAB,color=Cond)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  facet_wrap(~StimLoc) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "",
       x = 'Devaluation on Pleasantness rating',
       y = 'Choice update')

pdf(file.path(FigDir,'Corr_Devaluation_plea_choice.pdf'),6,5)
print(corr1)
print(corr2)
dev.off()

cor.test(df_comb$ChoiceChangeAB,df_comb$Didx)
# as expected, sig positive corr b/t choice updates and devaluation
# effect on pleasantness ratings
# explored if this correlation was modulated by conditions or stim locs, but no
model_corr0 <- lmer(ChoiceChangeAB ~ Didx + (1|SubID), data = df_comb)
model_corr1 <- lmer(ChoiceChangeAB ~ Didx + Cond + (1|SubID), data = df_comb)
model_corr2 <- lmer(ChoiceChangeAB ~ Didx + StimLoc + (1|SubID), data = df_comb)

summary(model_corr0)
anova(model_corr0)
anova(model_corr0,model_corr1) # n.s.
anova(model_corr0,model_corr2) # n.s.

###############################
## Day1 TMS effect on probe choice
######## look at cTBS-sham & sham-sham comparison ###
# focus on posterior OFC stimulation
###############################

p1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  geom_boxplot(fatten = NULL,width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  labs(x = "", title = 'Post(A) - Pre',
       y = "Choice of Sated Odor (Sweet vs. Savory)") + 
  common +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(fatten = NULL,width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  labs(x = "", title = 'Post(B) - Pre',
       y = "Choice of Sated Odor (Sweet vs. Savory)") + 
  common +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(p1,p2,ncol = 2, nrow = 1)

# on set A
value_XZ_A = subset(summary_choice_corrected,
                    Cond %in% c('sham-sham','cTBS-sham') &
                    StimLoc=='Posterior')$ChoiceChangeA
wilcox.test(value_XZ_A,alternative = 'l') # lower than zero? Yes, p=0.00037
value_X_A = subset(summary_choice_corrected, Cond=='cTBS-sham'& 
                    StimLoc=='Posterior')$ChoiceChangeA
value_Z_A = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='Posterior')$ChoiceChangeA
wilcox.test(value_X_A,alternative = 'l') # lower than zero? Yes, p=0.0088
wilcox.test(value_Z_A,alternative = 'l') # lower than zero? Yes, p=0.011
wilcox.test(value_X_A,value_Z_A,paired = T,alternative = 'g') # n.s

# on set B
value_XZ_B = subset(summary_choice_corrected,
                    Cond %in% c('sham-sham','cTBS-sham') &
                      StimLoc=='Posterior')$ChoiceChangeB
wilcox.test(value_XZ_B,alternative = 'l') # lower than zero? p = 0.064
value_X_B = subset(summary_choice_corrected, Cond=='cTBS-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeB
value_Z_B = subset(summary_choice_corrected, Cond=='sham-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeB
wilcox.test(value_X_B,alternative = 'l') # lower than zero? No, p=0.29
wilcox.test(value_Z_B,alternative = 'l') # lower than zero? No, p=0.092
wilcox.test(value_X_B,value_Z_B,paired = T,alternative = 'g') # n.s

tmp_dat = as.numeric(cbind(value_X_A-value_Z_A,
                value_X_B-value_Z_B))

tmp_df = data.frame('Set'=rep(c('A','B'),23),
                    'ChoiceSated'=tmp_dat)
tmp_df %>%
  ggplot(aes(x=Set,y=ChoiceSated)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4)
  #geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +

# models are on the trial level
# with session-wise baseline odor preference as a covariate
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Set=='A' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition
summary(model_choice_0)
coef_summary <- summary(model_choice_1)
coef_est <- coef_summary$coefficients["CondcTBS-sham", "Estimate"]
coef_se <- coef_summary$coefficients["CondcTBS-sham", "Std. Error"]
z_value <- coef_est / coef_se
p_value <- 1 - pnorm(z_value)
cat("Z-value:", z_value, "\n")
cat("One-sided p-value:", p_value, "\n")

# on Set B
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Set=='B' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p=0.075
summary(model_choice_1)
coef_summary <- summary(model_choice_1)
coef_est <- coef_summary$coefficients["CondcTBS-sham", "Estimate"]
coef_se <- coef_summary$coefficients["CondcTBS-sham", "Std. Error"]
z_value <- coef_est / coef_se
p_value <- 1 - pnorm(z_value)
cat("Z-value:", z_value, "\n")
cat("One-sided p-value:", p_value, "\n")

# both sets
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_2 <- glmer(Choice ~ (1|Set) + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ (1|Cond) + (1|Set) + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_0) # yes
anova(model_choice_1,model_choice_0) # yes
anova(model_choice_1,model_choice_2) # no
summary(model_choice_1)


# want to add conditioning variables into the linear models
# before that, correlate conditioning performance with devaluation effect
# question: does devaluation effect on choice correlate with conditioning performance?
# here I used run5 performance, but need to change to reflect the whole 5 runs

df_cdfinal_bothSets = conditioning_dat %>%
  subset(Run == 5) %>%
  group_by(SubID,Cond) %>%
  dplyr::summarise(cdfinal=mean(OdorChosen,na.rm = T))

df_comb_choice_cd = cbind(summary_choice_corrected,
                          df_cdfinal_bothSets$cdfinal)
colnames(df_comb_choice_cd)[7]='cdfinal'


# even if I use the final run, I should not calculate 
# conditioning for each set, should be for each cue
df_cdfinal = conditioning_dat %>%
  subset(Run == 5) %>%
  group_by(SubID,Cond,Set) %>%
  dplyr::summarise(cdfinal=mean(OdorChosen,na.rm = T))

# insert this cdfinal value to the whole choice dat ss df
choice_dat_ss_w_base_cdfinal = choice_dat_ss_w_base %>%
  mutate(cdfinal=NA)
for(i in 1:nrow(choice_dat_ss_w_base_cdfinal)){
  choice_dat_ss_w_base_cdfinal$cdfinal[i] = 
    df_cdfinal$cdfinal[
      df_cdfinal$SubID==choice_dat_ss_w_base_cdfinal$SubID[i] &
        df_cdfinal$Cond==choice_dat_ss_w_base_cdfinal$Cond[i] &
        df_cdfinal$Set==choice_dat_ss_w_base_cdfinal$Set[i]] 
}

# rerun models with cdfinal variable
# models are on the trial level
# with session-wise baseline odor preference as a covariate
use.dat = subset(choice_dat_ss_w_base_cdfinal, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Set=='A' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + cdfinal + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + cdfinal + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition
summary(model_choice_0)


use.dat = subset(choice_dat_ss_w_base_cdfinal, PrePost=='Post' &
                   StimLoc=='Posterior' &
                   Set=='B' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + logistic(cdfinal) + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + logistic(cdfinal) + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p=0.079
summary(model_choice_0)
summary(model_choice_1)

logistic <- function(x) {
  x[x==1] = 0.9999
  log(x / (1 - x))
}


##############################
# another set of analysis #
# only look at choice with one odor and one air
choice_dat_one = choice_dat %>%
  subset(ChoiceType<3) %>%
  mutate(OdorOption=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(IfDevalue=ifelse(OdorOption==Devalued,'Devalued','Non-Dv')) %>%
  mutate(IfDevalue=factor(IfDevalue,levels=c('Non-Dv','Devalued')))


  # mutate(Chosen = ifelse(ChoiceType == 1 & OdorChosen == 1, 'W', 
  #                        ifelse(ChoiceType==2 & OdorChosen==1,'A','N')))

summary_choice_one = choice_dat_one %>%
  group_by(SubID,Cond,PrePost,Set,StimLoc,IfDevalue) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T))

ggplot(summary_choice_one,
       aes(x=IfDevalue,y=Choice,fill=interaction(Set,PrePost))) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(aes(color=interaction(Set,PrePost)),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.3), 
              size = 1, alpha = 0.8) +
  labs(x = "", y = "Choosing Odor (vs. air)")

ggplot(summary_choice_one,
       aes(x=IfDevalue,y=Choice,fill=interaction(Set,PrePost))) +
  geom_boxplot(alpha=0.5) +
  facet_wrap(~ StimLoc + Cond) +
  labs(x = "", y = "Choosing Odor (vs. air)")




