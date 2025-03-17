
rm(list = ls())
source('Setup.R')

###############################
#### analyze conditioning data ####
###############################

conditioning_dat = read.xlsx('../ProcessedData/Conditioning_NODEAP_17.xlsx') %>%
  mutate(SelectedCue = if_else(SelectedCue == 0, NA_integer_, SelectedCue)) %>% # replace zero with NA
  mutate(Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('Sweet','Savory')),
         StimOrder=mapvalues(StimOrder,from=c('213'),
                             to=c('SC-CS-SS')),
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
         Sess=factor(Sess),
         Set=factor(Set),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         ChoiceType=factor(ChoiceType),
         StimLoc=factor(StimLoc))

conditioning_dat %>%
  group_by(Run,Cond) %>%
  dplyr::summarise(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,group=Cond,color=Cond)) +
  geom_line() +
  labs(x = "Runs", y = "P(Choosing Odors)",
       title = "Conditioning (NODEAP_17)") + common

### pleasanteness

dat = read.xlsx('../ProcessedData/Pleasant_NODEAP_17.xlsx') %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre-meal','Post-meal')),
         IfDevalue=mapvalues(IfDevalue,from=c('0','1'),to=c('Non-sated','Sated')),
         StimLoc=mapvalues(StimLoc,from=c('Anterior','Posterior'),to=c('aOFC','pOFC')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre-meal','Post-meal')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         IfDevalue=factor(IfDevalue,levels=c('Non-sated','Sated'))) %>%
  arrange(SubID,Cond,IfDevalue)

ggplot(dat,aes(x=IfDevalue,y=Pleasant,fill=PrePost)) +
  geom_jitter(aes(color=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
              size = 2, alpha = 0.8) +
  labs(x = "Odors", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))


###############################
# analyze choice data in the pre-meal
# and post-meal
###############################

choice_dat = read.xlsx('../ProcessedData/Choices_NODEAP_17.xlsx') %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre','Post')),
         Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
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

# summary of the choice_ss data
summary_choice_ss = choice_dat_ss %>%
  group_by(SubID,StimOrder,StimLoc,Cond,Sess,PrePost,Set) %>%
  dplyr::summarise(Choice=mean(Choice,na.rm = T))


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

avg_changeAB = (post_choice_A + post_choice_B)/2 - pre_choice
avg_changeB = post_choice_B - pre_choice
avg_changeA = post_choice_A - pre_choice
avg_BoverA = post_choice_B - post_choice_A

# create a corrected summary dataset
# containing choice update for each subject and each TMS condition
summary_choice_corrected = summary_choice_ss %>%
  ungroup() %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  select(SubID,StimLoc,StimOrder,Cond,Sess) %>%
  mutate(Sess=factor(Sess),
          ChoiceChangeAB=avg_changeAB, # both A and B
         ChoiceChangeA=avg_changeA, # only A
         ChoiceChangeB=avg_changeB, # only B
         ChoiceBoverA=avg_BoverA)

# set A and B 
summary_choice_corrected %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_jitter(aes(color=Sess,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre',
       y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")








