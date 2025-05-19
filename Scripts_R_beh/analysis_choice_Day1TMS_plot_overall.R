

# this code continues with 'analysis_choice.R' code
# with a focus on Day 1 TMS effect on choice

rm(list = ls())
source('Setup.R')
load(file = '../ProcessedData/Summary_Choice_corrected_dat.RData')

#######################################################
# plot using summarized choice data on day 1 TMS effect
# by comparing 'sham-sham' and 'cTBS-sham'
#######################################################

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',face = "bold",size = 16))
summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill=Cond),width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = NULL, y = "Choice of sated odor \n (post-pre meal)") + 
  common +
  coord_cartesian(ylim = c(-1,1.2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme(legend.position = "none")

# note: this one plots the summary stats, but we want the variance
# explained by TMS, after accounting for variance explained by 
# other covariates (value-diff mainly)


