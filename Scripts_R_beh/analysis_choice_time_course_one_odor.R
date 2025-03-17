

rm(list = ls())
source('Setup.R')

#######################################
## update 12/10/24, look at odor-air choice pairs
# reviewed 1/3/2025

# perhaps those trials are more informative?
# pre-meal phase
# odors are labeled as sated & non-sated
# (although not meaningful before the meal)

load(file = '../ProcessedData/choice_dat.RData')
choice_dat_one_pre = choice_dat %>%
  subset(PrePost=='Pre' & ChoiceType<3) %>%
  mutate(ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(OdorType=ifelse(ChoiceType==Devalued,'sated','non-sated')) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
       Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) 

# 12 such trials per session
# split choices into 6 parts and take each mean
# y-axis is whether or not choosing odor, regardless of odor type
choice_dat_one_pre_summary = choice_dat_one_pre %>%
  mutate(TrialCount=rep(1:12,nrow(choice_dat_one_pre)/12)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 6)) %>%
  group_by(SubID, StimLoc, Cond, part) %>%
  reframe(avg_choice = mean(ChosenOdor,na.rm=T)) 

# see if there's learning on selecting odor-predictive cues
k1 = choice_dat_one_pre_summary %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~StimLoc) +
  labs(y='Selecting odor',
       x='Trials',
       title ='Pre-meal task') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "ribbon", alpha = 0.3) 

# did not see learning, but there appears to be some difference
# regarding aOFC/pOFC

pdf(file.path(FigDir,'Pre_meal_choices_one_odor_time_course.pdf'),10,5)
print(k1)
dev.off()

# then also looked at individual results, but no clear pattern identified
pre_aOFC = choice_dat_one_pre_summary %>%
  subset(StimLoc=='aOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting odor',x='Pre-meal task',
       title = 'aOFC stim') +
  scale_color_manual(values = use.col.conds) 

pre_pOFC = choice_dat_one_pre_summary %>%
  subset(StimLoc=='pOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting odor',x='Pre-meal task',
       title = 'pOFC stim') +
  scale_color_manual(values = use.col.conds) 

pdf(file.path(FigDir,'Pre_meal_choices_one_odor_subs.pdf'),10,7)
print(pre_aOFC)
print(pre_pOFC)
dev.off()


# without a clear learning pattern identified
# just look at all trials
choice_dat_one_pre_summary_all = choice_dat_one_pre %>%
  group_by(SubID, StimLoc, Cond, Cond_day1, Cond_day2) %>%
  reframe(avg_choice = mean(ChosenOdor,na.rm=T)) 
pp1 = choice_dat_one_pre_summary_all %>%
  ggplot(aes(x=StimLoc,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~Cond) +
  geom_jitter() +
  geom_boxplot(alpha=0.5) +
  labs(y='P(Selecting odor)',
       x=NULL,
       title ='Pre-meal choices') +
  stat_compare_means(aes(group = StimLoc), 
                     method = "t.test", 
                     label = "p.format", 
                     label.y = 0.2,
                     paired = F) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common

pp2 = choice_dat_one_pre_summary_all %>%
  ggplot(aes(x=Cond,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~StimLoc) +
  geom_jitter(show.legend = F) +
  geom_boxplot(alpha=0.5,show.legend = F) +
  labs(y='P(Selecting odor)',
       x=NULL,
       title ='Pre-meal choices') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common

# collapse sham-sham and sham-cTBS here
pp3 = choice_dat_one_pre_summary_all %>%
  ggplot(aes(x=Cond_day1,y=avg_choice,color=Cond_day1,fill=Cond_day1)) +
  facet_wrap(~StimLoc) +
  geom_jitter() +
  geom_boxplot(alpha=0.5) +
  labs(y='P(Selecting odor)',
       x=NULL,
       title ='Pre-meal choices') +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) + common

model_choice_1 <- glmer(ChosenOdor ~ Cond_day1 * StimLoc + (1|SubID), 
                        data = choice_dat_one_pre,family = 'binomial')
model_choice_0 <- glmer(ChosenOdor ~ Cond_day1 + StimLoc + (1|SubID), 
                        data = choice_dat_one_pre,family = 'binomial')
anova(model_choice_0,model_choice_1)
summary(model_choice_1)
# sig. interaction effect between Cond-day1 and stimloc
# also sig. main effect of pOFC vs. aOFC
# not sig. sham vs. cTBS of day 1 TMS

use.dat = subset(choice_dat_one_pre,StimLoc=='aOFC')
model_choice_aOFC_1 <- glmer(ChosenOdor ~ Cond_day1 + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_aOFC_0 <- glmer(ChosenOdor ~ (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_aOFC_1,model_choice_aOFC_0)
summary(model_choice_aOFC_1)

use.dat = subset(choice_dat_one_pre,StimLoc=='pOFC')
model_choice_pOFC_1 <- glmer(ChosenOdor ~ Cond_day1 + (1|SubID), 
                             data = use.dat,family = 'binomial')
model_choice_pOFC_0 <- glmer(ChosenOdor ~ (1|SubID), 
                             data = use.dat,family = 'binomial')
anova(model_choice_pOFC_1,model_choice_pOFC_0)
summary(model_choice_pOFC_1)

# post-hoc analyses suggest this interaction was driven by pOFC, not
# by aOFC difference, making it even trickier to explain
# I think this is a dead end.

pdf(file.path(FigDir,'Pre_meal_choices_one_odor_all.pdf'),7,4)
print(pp1)
print(pp2)
print(pp3)
dev.off()


# still odor vs. air choices
# post-meal
# after the meal, it makes sense to separate such choices
# into sated and non-sated odors
choice_dat_one_post = choice_dat %>%
  subset(PrePost=='Post' & ChoiceType<3) %>%
  mutate(ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(OdorType=ifelse(ChoiceType==Devalued,'sated','non-sated'))

# 24 such trials per session
# split choices into 6 parts and take each mean
choice_dat_one_post_summary = choice_dat_one_post %>%
  mutate(TrialCount=rep(1:24,nrow(choice_dat_one_post)/24)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 6)) %>%
  group_by(SubID, StimLoc, Cond, part, OdorType) %>%
  reframe(avg_choice = mean(ChosenOdor,na.rm=T)) 

r = subset(choice_dat_one_post_summary,
           Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~StimLoc + OdorType) +
  labs(y='Selecting odor',
       x='Trials',
       title ='Post-meal task') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "ribbon", alpha = 0.3) 

pdf(file.path(FigDir,'Post_meal_choices_one_odor_time_course.pdf'),10,7)
print(r)
dev.off()

# forget about learning
# look at all post-meal trials
choice_dat_one_post_summary_all = choice_dat_one_post %>%
  group_by(SubID, StimLoc, Cond) %>%
  reframe(avg_choice = mean(ChosenOdor,na.rm=T)) 

qq1 = choice_dat_one_post_summary_all %>%
  ggplot(aes(x=StimLoc,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~Cond) +
  geom_jitter() +
  geom_boxplot(alpha=0.5) +
  labs(y='P(Selecting odor)',
       x=NULL,
       title ='Post-meal choices (odor vs air)') +
  stat_compare_means(aes(group = StimLoc), 
                     method = "t.test", 
                     label = "p.format", 
                     label.y = 0.2,
                     paired = F) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common

choice_dat_one_post_summary_all_byOdorType = choice_dat_one_post %>%
  group_by(SubID, StimLoc, Cond, OdorType) %>%
  reframe(avg_choice = mean(ChosenOdor,na.rm=T))

qq2 = choice_dat_one_post_summary_all_byOdorType %>%
  ggplot(aes(x=StimLoc,y=avg_choice,color=Cond,fill=Cond)) +
  facet_wrap(~OdorType + Cond) +
  geom_jitter() +
  geom_boxplot(alpha=0.5) +
  labs(y='P(Selecting odor)',
       x=NULL,
       title ='Post-meal choices (odor vs air)') +
  stat_compare_means(aes(group = StimLoc), 
                     method = "t.test", 
                     label = "p.format", 
                     label.y = 0.2,
                     paired = F) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common

pdf(file.path(FigDir,'Post_meal_choices_one_odor_all.pdf'),7,4)
print(qq1)
print(qq2)
dev.off()


####### include Day 1 learned cue value in the one-odor choice dataset
# here updated from ss choices, because we only have & only need w estimates 
# for the odor cues, while treating clean air cue as w=0

load(file = '../ProcessedData/choice_dat.RData')

choice_dat_one_pre = choice_dat %>%
  subset(SubID != 'NODEAP_17') %>% # this has weird thing when trying to find cue id
  subset(PrePost=='Pre' & ChoiceType<3) %>%
  mutate(ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(OdorType=ifelse(ChoiceType==Devalued,'sated','non-sated'))

load(file = '../ProcessedData/conditioning_cue_mapping.RData')
df_w = read.csv(file = '../BehModeling/df_w_per_sub_sess_cuepair.csv')

find_cue_id = function(odorcue,sub,sess){
  # input: odorcue: cue id in choice task
  # output: cue index during conditioning
  find_dat = subset(conditioning_cue_mapping,CueLeft==odorcue & SubID==sub & Sess==sess)
  return(find_dat$CuePair)
}

find_ws = function(sub,sess,Cue_id_conditioning){
  find_dat = subset(df_w,SubID==sub & Sess==sess & CuePair==Cue_id_conditioning)
  find_dat$w
}

choice_dat_one_pre$ValueOdorCue = NA

for(i in 1:nrow(choice_dat_one_pre)){
  curr_dat = choice_dat_one_pre[i,]
  odorcue = ifelse(curr_dat$OdorLR==1,curr_dat$CueLeft,curr_dat$CueRight)
  Cue_id_conditioning = find_cue_id(odorcue,curr_dat$SubID,curr_dat$Sess)
  w = find_ws(curr_dat$SubID,curr_dat$Sess,Cue_id_conditioning)
  choice_dat_one_pre$ValueOdorCue[i] = w
}

p1 = choice_dat_one_pre %>%
  mutate(ChosenOdor=factor(ChosenOdor)) %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x=ChosenOdor,y=ValueOdorCue)) +
  facet_wrap(~StimLoc + Cond) +
  geom_boxplot(alpha=0.3,outlier.alpha = 0) +
  geom_jitter(alpha=0.5) +
  labs(y='Odor cue values (from Day 1)',
       x='Selecting odor cue',
       title ='Pre-meal choices (odor vs air)') +
  stat_compare_means(aes(group = ChosenOdor), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     label.y.npc = 0.2,
                     color = 'red',
                     paired = F) 

# test it
# aOFC subjects only
use.dat = choice_dat_one_pre %>%
  subset(StimLoc=='aOFC') %>%
  filter(complete.cases(.)) %>%
  mutate(ChosenOdor=factor(ChosenOdor))
model_2 <- lmer(ValueOdorCue ~ Cond * ChosenOdor + (1|SubID), data = use.dat)
model_1 <- lmer(ValueOdorCue ~ Cond + ChosenOdor + (1|SubID), data = use.dat)
model_0 <- lmer(ValueOdorCue ~ ChosenOdor + (1|SubID), data = use.dat)
anova(model_0,model_1)
anova(model_0,model_2)
anova(model_1,model_2)
summary(model_2)
summary(model_1)

# pOFC subjects only
use.dat = choice_dat_one_pre %>%
  subset(StimLoc=='pOFC') %>%
  filter(complete.cases(.)) %>%
  mutate(ChosenOdor=factor(ChosenOdor))
model_2 <- lmer(ValueOdorCue ~ Cond * ChosenOdor + (1|SubID), data = use.dat)
model_1 <- lmer(ValueOdorCue ~ Cond + ChosenOdor + (1|SubID), data = use.dat)
model_0 <- lmer(ValueOdorCue ~ ChosenOdor + (1|SubID), data = use.dat)
anova(model_0,model_1)
anova(model_0,model_2)
anova(model_1,model_2)
summary(model_2)
summary(model_1)

# 12 such trials per session
# split choices into 3 parts and take each mean
choice_dat_one_pre_summary = choice_dat_one_pre %>%
  mutate(TrialCount=rep(1:12,nrow(choice_dat_one_pre)/12)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 3)) 

# see the relationship between Day 1 learned value and odor choices
# change over time
# focus on cTBS-sham condition
use_choice_dat_one_pre_summary = choice_dat_one_pre_summary %>%
  mutate(ChosenOdor=factor(ChosenOdor)) %>%
  filter(complete.cases(ChosenOdor)) %>%
  subset(Cond=='cTBS-sham') 

p2 = use_choice_dat_one_pre_summary %>%
  ggplot(aes(x=ChosenOdor,y=ValueOdorCue)) +
  geom_boxplot(alpha=0.3,outlier.alpha = 0) +
  geom_jitter(alpha=0.5) +
  facet_wrap(~StimLoc+part,nrow = 2) +
  labs(y='Learned value (from Day 1)',
       x='Selecting odor',
       title ='Pre-meal task') +
  stat_compare_means(aes(group = ChosenOdor), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     label.y.npc = 0.1,
                     paired = F)

p3 = choice_dat_one_pre_summary %>%
  filter(complete.cases(ChosenOdor)) %>%
  subset(Cond=='cTBS-sham') %>%
  mutate(part=factor(part)) %>%
  ggplot(aes(x=ValueOdorCue,y=ChosenOdor,
             group = part, color = part, fill = part)) +
  geom_jitter(alpha=0.5,width = 0.1,height = 0) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, alpha = 0.3) + 
  facet_wrap(~StimLoc,nrow = 2)


pdf(file.path(FigDir,'Post_meal_choices_one_odor_w_values.pdf'),7,6)
print(p1)
print(p2)
print(p3)
dev.off()



