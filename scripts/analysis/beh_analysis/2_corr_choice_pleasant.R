
rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

#######################################################

# check choices with odor-ratings
# to see whether people actually made their choices 
# based on the odor preferences
# both before and after the meal

# from 48 subjects but one has only two sessions
# look at pre-meal and post-meal separately

load(file.path(processed_dir,'choice_dat.RData'))
load(file.path(processed_dir,'Odor_ratings_dat.RData'))
     
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
  ggplot(aes(x=Plea_Diff,y=avg_choice)) +
  geom_point(aes(colour = Cond, fill=Cond),size=2) +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(y='P(choose sated)',
       x='Odor pleasantness (sated − non-sated), pre-meal',
       title = NULL) + 
  scale_color_manual(values = use.col.conds) + 
  scale_fill_manual(values = use.col.conds) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),  
    limits = c(0, 1)  # ensure the axis runs from 0 to 1
  ) +
  scale_x_continuous(
    breaks = seq(-15, 10, by = 5),  
    limits = c(-15, 15)  
  ) +
  common + theme(legend.position = 'none')


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
  ggplot(aes(x=Plea_Diff,y=avg_choice)) +
  geom_point(aes(colour = Cond, fill=Cond),size=2) +
  geom_smooth(method = rlm, se = TRUE, alpha = 0.2, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(y='P(choose sated)',
       x='Odor pleasantness (sated − non-sated), post-meal',
       title = NULL) + 
  scale_color_manual(values = use.col.conds) + 
  scale_fill_manual(values = use.col.conds) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),  
    limits = c(0, 1)  # ensure the axis runs from 0 to 1
  ) + 
  scale_x_continuous(
    breaks = seq(-15, 10, by = 5),  
    limits = c(-15, 15)  
  ) +
  common +
  theme(legend.position = 'none')

###############################
### correlate the choice updating with selective satiation idx ###
###############################

load(file.path(processed_dir,'Summary_Choice_corrected_dat.RData'))
load(file.path(processed_dir,'SelectSate_dat.RData'))

# (having checked rows are aligned)
df_comb = cbind(summary_choice_corrected,
                SelectSate_dat[,c('Devalued','Didx')])

p_corr=ggplot(df_comb,aes(x=Didx,y=ChoiceChangeAB)) +
  geom_point(aes(color=Cond,fill=Cond),size=2) +
  geom_smooth(method = rlm, se = TRUE, color = 'black') + 
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_cor(method = 'pearson', size = 5) +
  common +
  labs(title = NULL,
       x = 'Selective satiation index',
       y = 'ΔP(choose sated) (post − pre)') +
  theme(legend.position = 'none')

use_dat = df_comb %>%
  subset(StimLoc=='pOFC') %>%
  subset(Cond %in% c('sham-sham','sham-cTBS'))
model0 <- lmer(ChoiceChangeAB ~ Didx + (1|SubID), data = df_comb)
model1 <- lmer(ChoiceChangeAB ~ Didx + Cond + (1|SubID), data = df_comb)
anova(model0,model1)


pdf(file.path(FigPaperDir,'Corr_choice_w_odor_ratings_pre_meal.pdf'),6,5)
print(ggMarginal(corr_pre,type = 'density',groupFill = T))
dev.off()

pdf(file.path(FigPaperDir,'Corr_choice_w_odor_ratings_post_meal.pdf'),6,5)
print(ggMarginal(corr_post,type = 'density',groupFill = T))
dev.off()

pdf(file.path(FigPaperDir,'Corr_SelectSate_choice_by_Cond.pdf'),6,5)
print(ggMarginal(p_corr,type = 'density',groupFill = T))
dev.off()


