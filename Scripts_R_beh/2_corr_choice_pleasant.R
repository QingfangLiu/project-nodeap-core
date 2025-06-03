
rm(list = ls())
source('Setup.R')

#######################################################

# check choices with odor-ratings
# to see whether people actually made their choices 
# based on the odor preferences
# both before and after the meal

# from 48 subjects but one has only two sessions
# look at pre-meal and post-meal separately

load(file = '../ProcessedData/choice_dat.RData')
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


pdf(file.path(FigPaperDir,'Corr_choice_w_odor_ratings_pre_meal.pdf'),6,5)
print(ggMarginal(corr_pre,type = 'density',groupFill = T))
dev.off()

pdf(file.path(FigPaperDir,'Corr_choice_w_odor_ratings_post_meal.pdf'),6,5)
print(ggMarginal(corr_post,type = 'density',groupFill = T))
dev.off()


###############################
### correlate the choice updating with selective satiation idx ###
###############################

load(file = file.path(pro_dat_dir,'Summary_Choice_corrected_dat.RData'))
load(file = file.path(pro_dat_dir,'SelectSate_dat.RData'))

# (having checked rows are aligned)
df_comb = cbind(summary_choice_corrected,
                SelectSate_dat[,c('Devalued','Didx')])

corr1=ggplot(df_comb,aes(x=Didx,y=ChoiceChangeAB)) +
  geom_point(aes(color=Cond,fill=Cond)) +
  geom_smooth(method = rlm, se = TRUE, 
              linetype = 2) + 
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_cor(method = 'pearson') +
  common +
  labs(title = NULL,
       x = 'Selective satiation',
       y = 'Choice of sated odor (post - pre)')

# pOFC only
# focus on Day 2 TMS effect
corr2=df_comb %>%
  subset(StimLoc=='pOFC') %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Didx,y=ChoiceChangeAB)) +
  geom_point(aes(color=Cond,fill=Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  #facet_wrap(~StimLoc) +
  geom_smooth(method = rlm, se = T, color = 'black',
              alpha = 0.2, linetype = 3) + common +
  stat_cor(method = 'pearson') +
  labs(title = NULL,
       x = 'Selective satiation',
       y = 'Choice of sated odor \n (post - pre-meal)') +
  theme(legend.position = 'topright')

pdf(file.path(FigDir,'Corr_SelectSate_choice.pdf'),5,5)
print(corr1)
dev.off()

pdf(file.path(FigPaperDir,'Corr_SelectSate_choice_by_Cond.pdf'),5,4)
ggMarginal(corr2,type = 'density',groupFill = T)
dev.off()

use_dat = df_comb %>%
  subset(StimLoc=='pOFC') %>%
  subset(Cond %in% c('sham-sham','sham-cTBS'))
model0 <- lmer(ChoiceChangeAB ~ Didx + (1|SubID), data = df_comb)
model1 <- lmer(ChoiceChangeAB ~ Didx + Cond + (1|SubID), data = df_comb)
anova(model0,model1)






