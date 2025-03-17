
# this was copied and modifed from the other analysis script
# modifed from the code to Thorsten

library(openxlsx)
library(plyr)
library(tidyverse)
library(lme4)
library(ggpubr) 

FigDir = '../Figs/'
pd <- position_dodge(0.6)
common = theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) 

###############################
#### analyze conditioning data ####
###############################

conditioning_dat = read.xlsx('../ProcessedData/Conditioning.xlsx') %>%
  mutate(Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('Sweet','Savory')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Set=factor(Set),
         ChoiceType=factor(ChoiceType),
         StimLoc=factor(StimLoc))

# if I create new conditions based on Day1 & Day2 TMS separately
conditioning_dat= conditioning_dat %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2))

p1=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  subset(StimLoc %in% 'Posterior') %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  scale_color_manual(values = c("sham" = "gray","cTBS" = "green")) +
  scale_fill_manual(values = c("sham" = "gray","cTBS" = "green")) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  labs(x = "Runs", y = "Selecting odors", color = "", fill = "") +
  common +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7,0.3),
        legend.title =  element_blank(),
        strip.background = element_blank())

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

# to plot the change of pleasantness ratings (post meal-pre meal)
# create a new df 'reduced' from 'dat'
pre_plea = subset(dat,PrePost=='Pre')$Pleasant
post_plea = subset(dat,PrePost=='Post')$Pleasant
change_plea = post_plea - pre_plea
reduced = subset(dat,PrePost=='Pre')
reduced$PrePost = 'Change'
reduced$Pleasant = change_plea

reduced = reduced %>% 
  mutate(pair=rep(c(1:(nrow(reduced)/2)),each=2)) %>%
  subset(StimLoc %in% 'Posterior')

p2=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.2, color = 'darkgray') +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_x_discrete(labels = c("Non-Dv" = "non-sated", "Devalued" = "sated")) +
  labs(x = "Odors", y = "Pleasantness ratings (post-pre)") + 
  common

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

# only look at choice with one sweet and one savory (denoted as ss)
# perhaps where the effect is the strongest
choice_dat_ss = choice_dat %>%
  subset(ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0))

summary_choice = choice_dat_ss %>%
  group_by(SubID,Cond,PrePost,Set,StimLoc) %>%
  dplyr::summarise(Choice=mean(Choice,na.rm = T))

# using pre-A as baseline
# get difference from pre to post
pre_choice = subset(summary_choice,Set=='A' & PrePost=='Pre')$Choice
post_choice_A = subset(summary_choice,Set=='A' & PrePost=='Post')$Choice
post_choice_B = subset(summary_choice,Set=='B' & PrePost=='Post')$Choice
conds = subset(summary_choice,Set=='A' & PrePost=='Pre')$Cond
# check if there is pre-post correlations of choices
df_pre_post_summary = data.frame(
  'conds'=conds,
  'pre_choice'=pre_choice,
  'post_choice_A'=post_choice_A,
  'post_choice_B'=post_choice_B
)

avg_changeAB = (post_choice_A + post_choice_B)/2 - pre_choice
avg_changeB = post_choice_B - pre_choice
avg_changeA = post_choice_A - pre_choice

summary_choice_corrected = summary_choice %>%
  ungroup() %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  select(SubID,StimLoc,Cond) %>%
  mutate(ChoiceChangeAB=avg_changeAB, # both A and B
         ChoiceChangeA=avg_changeA, # only A
         ChoiceChangeB=avg_changeB) # only B

p3=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  subset(StimLoc %in% 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  labs(x = "Day2 TMS\n (Day1 TMS:sham)", title = '',
       y = "Selecting sated odor (post - pre)") + 
  scale_x_discrete(labels = c("sham-sham" = "sham", "sham-cTBS" = "cTBS")) +
  common + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())


p4=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  labs(x = "Day1 TMS\n (Day2 TMS:sham)", title = '',
       y = "Selecting sated odor (post - pre)") + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_x_discrete(labels = c("sham-sham" = "sham", "cTBS-sham" = "cTBS")) +
  common +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  
pdf(file.path(FigDir,'NODEAPSummary_K99.pdf'), 14, 4)
fig = ggpubr::ggarrange(p2,p1,p3,p4,ncol = 4, nrow = 1, widths = c(0.8,1,0.8,0.8))
print(fig)
dev.off()


