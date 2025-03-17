
rm(list = ls())
source('Setup.R')

#######################################################
# change of choice in the trial course on pre-meal task
# I realized that disruption of identity learning should
# also be shown in the pre-meal choice task (set A)
# (1) if we assume people have a session-wise pre-meal odor
# baseline, then people need to learn which stimuli to choose

# (2) if we believe such a learning process, then we
# can't use the whole trial course to serve as 
# pre-meal odor preference

# updated 12/10/24
# not only identity learning, since we found also value learning
# from Day 1 TMS, pre-meal might provide a remedy for value learning
# as well

#######################################################

# added: I can also look at the choices with odor-ratings
# to see whether people actually made their choices 
# based on the odor preferences
# both before and after the meal

load(file = '../ProcessedData/choice_dat.RData')

# only include SW/SA choices, 18 trials per run (set A only)
pre_meal_choice = choice_dat %>%
  subset(PrePost=='Pre' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) 
  # Choice==1 means selecting sated odor

# split choices into 3 parts and take each mean
pre_meal_choice = pre_meal_choice %>%
  mutate(TrialCount=rep(1:18,nrow(pre_meal_choice)/18)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 3)) %>%
  group_by(SubID, StimLoc, Cond, part) %>%
  reframe(avg_choice = mean(Choice,na.rm=T)) 

pre1=pre_meal_choice %>%
  subset(StimLoc=='aOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Pre-meal task',
       title = 'aOFC stim') +
  scale_color_manual(values = use.col.conds) 

pre2=pre_meal_choice %>%
  subset(StimLoc=='pOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Pre-meal task',
       title = 'pOFC stim') + 
  scale_color_manual(values = use.col.conds) 

pdf(file.path(FigDir,'Pre_meal_choices_subs.pdf'),10,7)
print(pre1)
print(pre2)
dev.off()

#######################################################
# post-meal choices
#######################################################

# only include SW/SA choices, 36 trials per run (set A & B)
post_meal_choice = choice_dat %>%
  subset(PrePost=='Post' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) 
  # Choice==1 means selecting sated odor

# split choices into 3 parts and take each mean
post_meal_choice = post_meal_choice %>%
  mutate(TrialCount=rep(1:36,nrow(post_meal_choice)/36)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 6)) %>%
  group_by(SubID, StimLoc, Cond, part) %>%
  reframe(avg_choice = mean(Choice,na.rm=T)) 

post_meal_choice %>%
  group_by(Cond,part,SubID,StimLoc) %>%
  reframe(Choice=mean(avg_choice,na.rm = T)) %>%
  ggplot(aes(x = part, y = Choice, color = Cond, fill = Cond)) +
  facet_wrap(~StimLoc) +
  stat_summary(fun = mean, geom = "line", size = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "ribbon", alpha = 0.3) + common +
  labs(y='Selecting sated odor',x='Post-meal task') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds)

post1=post_meal_choice %>%
  subset(StimLoc=='aOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Post-meal task',
       title = 'aOFC stim') + 
  scale_color_manual(values = use.col.conds) 

post2=post_meal_choice %>%
  subset(StimLoc=='pOFC') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) + 
  labs(y='Selecting sated odor',x='Post-meal task',
       title = 'pOFC stim') +
  scale_color_manual(values = use.col.conds) 

pdf(file.path(FigDir,'Post_meal_choices_subs.pdf'),10,7)
print(post1)
print(post2)
dev.off()


###### look at corr of choices with odor ratings
# from 48 subjects but one has only two sessions
# look at pre-meal and post-meal separately

load(file = '../ProcessedData/Odor_ratings_dat.RData')

pleasant_diff_pre <- Odor_ratings_dat %>%
  subset(PrePost=='Pre-meal') %>%
  group_by(SubID, Sess) %>%
  reframe(Plea_Diff = Pleasant[IfDevalue == "sated"] - 
            Pleasant[IfDevalue == "non-sated"])

pre_meal_choice = choice_dat %>%
  subset(PrePost=='Pre' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) %>%
  group_by(SubID, StimLoc, Sess, Cond) %>%
  reframe(avg_choice = mean(Choice,na.rm=T)) %>%
  mutate(Plea_Diff = pleasant_diff_pre$Plea_Diff) 
  
corr_pre = pre_meal_choice %>%
  ggplot(aes(x=Plea_Diff,y=avg_choice,color=Cond,fill=Cond)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
  stat_cor(method = 'pearson') +
  labs(y='Selecting sated odor',
       x='Odor pleasant (sated - non-sated)',
       title = 'Pre-meal') + 
  scale_color_manual(values = use.col.conds) + 
  scale_fill_manual(values = use.col.conds) +
  common +
  theme(legend.position = 'none')

corr_pre_by_sess = pre_meal_choice %>%
  mutate(Sess=factor(Sess)) %>%
  ggplot(aes(x=Plea_Diff,y=avg_choice,color=Sess,fill=Sess)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
  stat_cor(method = 'pearson') +
  labs(y='Selecting sated odor',
       x='Odor pleasant (sated - non-sated)',
       title = 'Pre-meal') + 
  scale_color_manual(values = use.col.sess) + 
  scale_fill_manual(values = use.col.sess) +
  common +
  theme(legend.position = 'none')


## post-meal
pleasant_diff_post <- Odor_ratings_dat %>%
  subset(PrePost=='Post-meal') %>%
  group_by(SubID, Sess) %>%
  reframe(Plea_Diff = Pleasant[IfDevalue == "sated"] - 
            Pleasant[IfDevalue == "non-sated"])

post_meal_choice = choice_dat %>%
  subset(PrePost=='Post' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) %>%
  group_by(SubID, StimLoc, Sess, Cond) %>%
  reframe(avg_choice = mean(Choice,na.rm=T)) %>%
  mutate(Plea_Diff = pleasant_diff_post$Plea_Diff)

corr_post = post_meal_choice %>%
  ggplot(aes(x=Plea_Diff,y=avg_choice,color=Cond,fill=Cond)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
  stat_cor(method = 'pearson') +
  labs(y='Selecting sated odor',
       x='Odor pleasant (sated - non-sated)',
       title = 'Post-meal') + 
  scale_color_manual(values = use.col.conds) + 
  scale_fill_manual(values = use.col.conds) +
  common +
  theme(legend.position = 'none')

corr_post_by_sess = post_meal_choice %>%
  mutate(Sess=factor(Sess)) %>%
  ggplot(aes(x=Plea_Diff,y=avg_choice,color=Sess,fill=Sess)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
  stat_cor(method = 'pearson') +
  labs(y='Selecting sated odor',
       x='Odor pleasant (sated - non-sated)',
       title = 'Post-meal') + 
  scale_color_manual(values = use.col.sess) + 
  scale_fill_manual(values = use.col.sess) +
  common +
  theme(legend.position = 'none')

pdf(file.path(FigDir,'Corr_choice_w_odor_ratings_pre_meal.pdf'),6,5)
print(ggMarginal(corr_pre,type = 'density',groupFill = T))
print(ggMarginal(corr_pre_by_sess,type = 'density',groupFill = T),
      newpage = T)
dev.off()

pdf(file.path(FigDir,'Corr_choice_w_odor_ratings_post_meal.pdf'),6,5)
print(ggMarginal(corr_post,type = 'density',groupFill = T))
print(ggMarginal(corr_post_by_sess,type = 'density',groupFill = T),
      newpage = T)
dev.off()


##################################################
# still on choices w/t odor ratings, but separate choices into blocks
load(file = '../ProcessedData/choice_dat.RData')

# pre-meal choices
corr_pre = NULL
fig_name = 'Corr_choice_w_odor_ratings_pre_meal_byparts.pdf'
pdf(file.path(FigDir,fig_name),6,5)

for (p in 1:3){
# only include SW/SA choices, 18 trials per run (set A only)
pre_meal_choice = choice_dat %>%
  subset(PrePost=='Pre' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) 
# Choice==1 means selecting sated odor

# split choices into 3 parts and take each mean
pre_meal_choice = pre_meal_choice %>%
  mutate(TrialCount=rep(1:18,nrow(pre_meal_choice)/18)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 3)) %>%
  group_by(SubID, StimLoc, Sess, Cond, part) %>%
  reframe(avg_choice = mean(Choice,na.rm=T)) %>%
  subset(part == p) %>%
  mutate(Plea_Diff = pleasant_diff_pre$Plea_Diff) 

corr_pre[[p]] = pre_meal_choice %>%
  ggplot(aes(x=Plea_Diff,y=avg_choice,color=Cond,fill=Cond)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
  stat_cor(method = 'pearson') +
  labs(y= 'Selecting sated odor',
       x='Rated pleasant (sated - non-sated)',
       title = paste('Pre-meal part',p)) + 
  scale_color_manual(values = use.col.conds) + 
  scale_fill_manual(values = use.col.conds) +
  common +
  theme(legend.position = 'none')

pp=ggMarginal(corr_pre[[p]],type = 'density',
              groupFill = T)
if(p>1) {
  print(pp,newpage=T)
}else{
  print(pp)
}  
}

dev.off()


# post-meal choices
load(file = '../ProcessedData/choice_dat.RData')
corr_post = NULL
fig_name = 'Corr_choice_w_odor_ratings_post_meal_byparts.pdf'
pdf(file.path(FigDir,fig_name),6,5)

for (p in 1:3){
  # only include SW/SA choices, 36 trials per run (set A and B)
  post_meal_choice = choice_dat %>%
    subset(PrePost=='Post' & ChoiceType==3) %>%
    mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
    mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) 
  # Choice==1 means selecting sated odor
  
  # split choices into 6 parts and take each mean
  post_meal_choice = post_meal_choice %>%
    mutate(TrialCount=rep(1:36,nrow(post_meal_choice)/36)) %>%
    group_by(SubID, Cond) %>%
    mutate(part = ntile(TrialCount, 6)) %>%
    group_by(SubID, StimLoc, Sess, Cond, part) %>%
    reframe(avg_choice = mean(Choice,na.rm=T)) %>%
    subset(part == p) %>%
    mutate(Plea_Diff = pleasant_diff_post$Plea_Diff) 
  
  corr_post[[p]] = post_meal_choice %>%
    ggplot(aes(x=Plea_Diff,y=avg_choice,color=Cond,fill=Cond)) +
    geom_point() +
    geom_smooth(method = rlm, se = TRUE, alpha = 0.2) +
    stat_cor(method = 'pearson') +
    labs(y= 'Selecting sated odor',
         x='Rated pleasant (sated - non-sated)',
         title = paste('Post-meal part',p)) + 
    scale_color_manual(values = use.col.conds) + 
    scale_fill_manual(values = use.col.conds) +
    common +
    theme(legend.position = 'none')
  
  pp=ggMarginal(corr_post[[p]],type = 'density',
                groupFill = T)
  if(p>1) {
    print(pp,newpage=T)
  }else{
    print(pp)
  }  
}

dev.off()





