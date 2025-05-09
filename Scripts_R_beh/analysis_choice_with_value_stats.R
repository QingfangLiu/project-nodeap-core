
# this code adds the learned cue weights as covariates in similar mixed-effect
# models to analyze choice data
# looked at Day 2 TMS effect first
# then Day 1 TMS effect

# this runs after 'analysis_choice_with_value_visua.R'

rm(list = ls())
source('Setup.R')

# load the choice (ss) with baseline and values
load(file = '../ProcessedData/choice_dat_ss_w_base_values.RData')

# create a var ValueDiff from difference of sated and non-sated
choice_dat_ss_w_base = use_choice_dat_ss %>%
  mutate(ValueDiff=ValueSated-ValueNonSated)

################## Day 2 TMS effect ##################

# focus on pOFC
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' & 
                   StimLoc=='pOFC' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ ValueDiff + Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)

# focus on aOFC
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' & 
                   StimLoc=='aOFC' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ ValueDiff + Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)


# put aOFC and pOFC together 
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ StimLoc * Cond + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ StimLoc + Cond + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)
summary(model_choice_1)
summary(model_choice_2)

# overall, the value diff covariate is sig in the regression models
# and adding this still keeps most of the Day 2 TMS effect intact

# !!!!!!!!!!!!!!!!!!!!
# to update the stat report in the paper, with this value diff 
# as covariate included
# !!!!!!!!!!!!!!!!!!!!


###### then how about Day 1 TMS effect? ##############
# note NODEAP_17 is a posterior cTBS subject
# with only sess 1 (sham-cTBS) and 3 (sham-sham) so remove this
# when analyzing day 1 TMS effect

# split subjects into pOFC and aOFC groups
use_trial_dat_pOFC = choice_dat_ss_w_base %>%
                  subset((!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         StimLoc=='pOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham')) %>%
  mutate(CuePair=factor(CuePair)) %>%
  filter(complete.cases(.))   # remove all NA trials

# check on aOFC-cTBS subjects
use_trial_dat_aOFC = choice_dat_ss_w_base %>%
                            subset(PrePost=='Post' &
                            StimLoc=='aOFC' & 
                            Cond %in% c('cTBS-sham','sham-sham')) %>%
  mutate(CuePair=factor(CuePair)) %>%
  filter(complete.cases(.))   # remove all NA trials





