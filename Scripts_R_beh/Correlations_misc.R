
rm(list = ls())
source('Setup.R')

pro_dat_dir = '/Users/liuq13/Library/CloudStorage/Box-Box/NODEAP_data_analysis/ProcessedData'

# this code correlates the choices with odor pleasantness ratings
# as well as with TMS ratings

###############################
### correlate the choice updating with selective satiation idx ###
###############################

load(file = file.path(pro_dat_dir,'Summary_Choice_corrected_dat.RData'))
load(file = file.path(pro_dat_dir,'SelectSate_dat.RData'))

# (having checked rows are aligned)
df_comb = cbind(summary_choice_corrected,
                SelectSate_dat[,c('Devalued','Didx')])
corr1=ggplot(df_comb,aes(x=Didx,y=ChoiceChangeAB)) +
  geom_point() +
  geom_smooth(method = rlm, se = TRUE, color = 'black',
              linetype = 2) + common +
  labs(title = NULL,
       x = 'Selective satiation',
       y = 'Choice of sated odor (post - pre)')

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


###############################
### correlate the choice updating with TMS ratings
###############################

###### Day 2
summary_use = summary_choice_corrected %>%
  arrange(SubID,StimLoc,Sess)

survey_dat = read.xlsx('../ProcessedData/Survey_uncomf_strong.xlsx') %>%
  subset(TMSSess %in% c(2,4,6)) %>%
  subset(!(SubID=='NODEAP_17' & sess_name=='S2D2')) %>%
  mutate(rating = (uncom + strong)/2) %>%
  mutate(ChoiceChangeAB = summary_use$ChoiceChangeAB,
         Cond = summary_use$Cond) %>%
  subset(Cond %in% c('sham-sham','sham-cTBS'))

for(roi in c('aOFC','pOFC')){
  if(roi=='aOFC') ylab_text = 'P(SA) \n (post - pre-meal)'
  if(roi=='pOFC') ylab_text = ''
  
  use_df = survey_dat %>% subset(StimLoc==roi)
  corr3 = use_df %>%
  ggplot(aes(x=rating,y=ChoiceChangeAB,color=Cond,fill=Cond)) +
  geom_point() +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_smooth(method = rlm, se = T, 
              alpha = 0.2, linetype = 3) + common +
  stat_cor(method = 'pearson',
           label.x = -1.8,
           label.y = c(0.65,0.8)) +
  labs(x = 'Perceived TMS discomfort & intensity',
       y = ylab_text) +
  coord_cartesian(xlim = c(-1.8,2),ylim = c(-1,1.1)) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, 
    ymin = 1.1, ymax = 0.9,
    fill = use.col.ap.ofc[roi], alpha = 1
  ) +
  annotate("text", x = 0, y = 1, color = 'white',
    label = roi, size = 6, fontface = "bold"
  ) +
  theme(legend.position = 'topright')

figname = paste0('Corr_choice_w_Day2TMS_ratings_',roi,'.pdf')
pdf(file.path(FigDir,figname),5,4)
print(ggMarginal(corr3,type = 'density',groupFill = T))
dev.off()
}

# Subset the data to only include the two conditions and calculate differences
diff_data <- survey_dat %>%
  group_by(SubID) %>%
  reframe(
    StimLoc = StimLoc,
    rating_diff = rating[Cond == 'sham-cTBS'] - rating[Cond == 'sham-sham'],
    ChoiceChangeAB_diff = ChoiceChangeAB[Cond == 'sham-cTBS'] - ChoiceChangeAB[Cond == 'sham-sham']
  )

# Plot the correlation of the differences

corr4 = diff_data %>%
  ggplot(aes(x = rating_diff, y = ChoiceChangeAB_diff, 
             color = StimLoc, fill = StimLoc)) +
  geom_point(color = 'black', show.legend = F) +
  geom_smooth(method = 'lm', se = TRUE, alpha = 0.3,
              linetype = 'dashed', show.legend = F) + 
  stat_cor(method = 'pearson',show.legend = F) +
  scale_color_manual(values = use.col.ap.ofc) +
  scale_fill_manual(values = use.col.ap.ofc) +
  labs(title = NULL,
       x = "Perceived TMS discomfort & intensity \n Day 2 (cTBS - sham)",
       y = "Sated odor update \n (cTBS - sham)") +
  common

figname = 'Corr_choice_w_Day2TMS_ratings_comp.pdf'
pdf(file.path(FigDir,figname),5,4)
print(ggMarginal(corr4,type = 'density'))
dev.off()


###### Day 1
summary_use_Day1TMS = summary_choice_corrected %>%
  subset(!SubID=='NODEAP_17') %>%
  arrange(SubID,StimLoc,Sess)

survey_dat = read.xlsx('../ProcessedData/Survey_uncomf_strong.xlsx') %>%
  subset(TMSSess %in% c(1,3,5)) %>%
  subset(!SubID=='NODEAP_17') %>%
  mutate(rating = (uncom + strong)/2) %>%
  mutate(ChoiceChangeAB = summary_use_Day1TMS$ChoiceChangeAB,
         Cond = summary_use_Day1TMS$Cond) %>%
  subset(Cond %in% c('sham-sham','cTBS-sham'))

for(roi in c('aOFC','pOFC')){
  
  if(roi=='aOFC') ylab_text = 'P(SA) \n (post - pre-meal)'
  if(roi=='pOFC') ylab_text = ''
  
  use_df = survey_dat %>% subset(StimLoc==roi)
  corr5 = use_df %>%
  ggplot(aes(x=rating,y=ChoiceChangeAB,fill = Cond,color=Cond)) +
  geom_point(aes(color=Cond,fill=Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_smooth(method = rlm, se = T, color = 'black',
              alpha = 0.2, linetype = 3) + common +
  stat_cor(method = 'pearson', label.x = -1.8,
           label.y = c(0.65,0.8)) +
  labs(title = NULL,
       x = 'Perceived TMS discomfort & intensity',
       y = ylab_text) +
  coord_cartesian(xlim = c(-1.8,2),ylim = c(-1,1.1)) +
  annotate("rect",
    xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.1,
    fill = use.col.ap.ofc[roi], alpha = 1
  ) +
  annotate(
    "text", x = 0, y = 1, color = 'white',
    label = roi, size = 6, fontface = "bold"
  ) +
  theme(legend.position = 'topright')

figname = paste0('Corr_choice_w_Day1TMS_ratings_',roi,'.pdf')
pdf(file.path(FigDir,figname),5,4)
print(ggMarginal(corr5,type = 'density',groupFill = T))
dev.off()
}

diff_data <- survey_dat %>%
  group_by(SubID) %>%
  reframe(
    StimLoc = StimLoc,
    rating_diff = rating[Cond == 'cTBS-sham'] - rating[Cond == 'sham-sham'],
    ChoiceChangeAB_diff = ChoiceChangeAB[Cond == 'cTBS-sham'] - 
      ChoiceChangeAB[Cond == 'sham-sham']
  )

corr6 = diff_data %>%
  ggplot(aes(x = rating_diff, y = ChoiceChangeAB_diff, 
             color = StimLoc, fill = StimLoc)) +
  geom_point(color = 'black', show.legend = F) +
  geom_smooth(method = 'lm', se = TRUE, alpha = 0.3,
              linetype = 'dashed', show.legend = F) + 
  stat_cor(method = 'pearson',show.legend = F) +
  scale_color_manual(values = use.col.ap.ofc) +
  scale_fill_manual(values = use.col.ap.ofc) +
  labs(title = NULL,
       x = "Perceived TMS discomfort & intensity \n Day 1 (cTBS - sham)",
       y = "") +
  common

figname = 'Corr_choice_w_Day1TMS_ratings_comp.pdf'
pdf(file.path(FigDir,figname),5,4)
print(ggMarginal(corr6,type = 'density'))
dev.off()


