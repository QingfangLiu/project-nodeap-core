
rm(list = ls())
source('Setup.R')

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

###############################
save(conditioning_dat,file = '../ProcessedData/Conditioning.RData')
###############################

# does conditioning improve over runs?

p01=conditioning_dat %>%
  subset(StimLoc=='Anterior') %>%
  mutate(Sess=factor(Sess)) %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Sess,group=Sess)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.sess) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Anterior Stim)") 
p02=conditioning_dat %>%
  subset(StimLoc=='Posterior') %>%
  mutate(Sess=factor(Sess)) %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Sess,group=Sess)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.sess) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Posterior Stim)") 

p03=conditioning_dat %>%
  subset(StimLoc=='Anterior') %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Cond_day1,group=Cond_day1)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Anterior Stim)") 

p04=conditioning_dat %>%
  subset(StimLoc=='Posterior') %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Cond_day1,group=Cond_day1)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Posterior Stim)") 

pdf(file.path(FigDir,'Conditioning_subs.pdf'),10,8)
print(p01)
print(p02)
print(p03)
print(p04)
dev.off()

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
levels(conditioning_dat$StimLoc)=c('aOFC','pOFC')
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
  labs(x = "Runs", y = "Choosing stimuli predicting odors", 
       color = "Day1 TMS") +
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

pdf(file.path(FigDir,'Conditioning_error_bar.pdf'),7,4)
print(pp2a)
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



