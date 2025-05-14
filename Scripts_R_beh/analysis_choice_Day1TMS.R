

# this code continues with 'analysis_choice.R' code
# with a focus on Day 1 TMS effect on choice

rm(list = ls())
source('Setup.R')
load(file = '../ProcessedData/Summary_Choice_corrected_dat.RData')
load(file = '../ProcessedData/choice_dat.RData')
load(file = '../ProcessedData/choice_dat_ss_w_base.RData')

#######################################################
# stat testing on Day 1 TMS effect
# test using summarized choice data
#######################################################

# need to exclude NODEAP_17 who has missing sessions 
# (this subject has sham-sham sham-cTBS that can contribute one data point
# to Day 2 TMS effect though)

# pOFC subjects
use_summ_dat = subset(summary_choice_corrected,
                      (!SubID=='NODEAP_17') &
                        Cond %in% c('sham-sham','cTBS-sham') & 
                        StimLoc=='pOFC')

#######################################################

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
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

#######################################################
# stat testing Day 1 TMS effect at trial level
#######################################################

survey_dat = read.xlsx('../ProcessedData/Survey_uncomf_strong.xlsx') %>%
  subset(TMSSess %in% c(1,3,5)) %>%
  subset(!SubID=='NODEAP_17') %>%
  mutate(rating = (uncom + strong)/2) %>%
  mutate(Sess=rep(1:3,47))

use_trial_dat = subset(choice_dat_ss_w_base, 
                       (!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         StimLoc=='pOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.)) %>%  # remove all NA trials
  mutate(rating = NA)

# add TMS rating variable to the trial-wise data
for(i in 1:nrow(use_trial_dat)){
  use_trial_dat$rating[i] = survey_dat$rating[survey_dat$SubID==use_trial_dat$SubID[i] &
                                                survey_dat$Sess==use_trial_dat$Sess[i]]
}

###############################
# on aOFC
use_trial_dat = subset(choice_dat_ss_w_base, 
                       PrePost=='Post' &
                         StimLoc=='aOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.))

######## look at cTBS-sham & sham-sham comparison ###
###############################

# p values come from glmer models

p_values_setA <- data.frame(
  StimLoc = c("aOFC", "pOFC"), 
  p_value = c(0.001852, 1) 
)

p_values_setB <- data.frame(
  StimLoc = c("aOFC", "pOFC"), 
  p_value = c(0.0351, p=0.0308) 
)

p_values_setA$label <- ifelse(p_values_setA$p_value < 0.05, 
                              paste("p =", format(p_values_setA$p_value, digits = 2)), 
                              "n.s.")
p_values_setB$label <- ifelse(p_values_setB$p_value < 0.05, 
                              paste("p =", format(p_values_setB$p_value, digits = 2)), 
                              "n.s.")

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))










