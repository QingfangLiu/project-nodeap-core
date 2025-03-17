
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


########################################################
# check upon set B first

# pOFC subjects
use_trial_dat_pOFC_B = subset(use_trial_dat_pOFC, Set=='B')
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition, p=0.1746
summary(model_choice_0)
summary(model_choice_1)

model_choice_2 <- glmer(Choice ~ Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
model_choice_3 <- glmer(Choice ~ ValueDiff + base + Cond + Sess + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
anova(model_choice_2,model_choice_3) # no effect of condition, p=0.1828
summary(model_choice_2)
summary(model_choice_3)
anova(model_choice_2,model_choice_0) # sess was not sig

# bummer! adding value difference actually showed no TMS effect on set B
# in the meanwhile the value-diff effect was very postively strong (as it should be)
# but the base was NOT significant

# maybe not using base this time because base comes from set A only
# and now we examine set A and B separately
# also the base was not strong in the above models

model_choice_5 <- glmer(Choice ~ Cond + ValueDiff + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
model_choice_4 <- glmer(Choice ~ ValueDiff + (1|SubID), 
                        data = use_trial_dat_pOFC_B,family = 'binomial')
anova(model_choice_5,model_choice_4) # still no effect of condition
summary(model_choice_4)
summary(model_choice_5)


# still set B
# aOFC subjects
use_trial_dat_aOFC_B = subset(use_trial_dat_aOFC, Set=='B')
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_B,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_B,family = 'binomial')
anova(model_choice_1,model_choice_0) 
summary(model_choice_1)
# marginal effect of condition
# strong value diff effect
# sig & positive base effect
# (this is in contrast with pOFC subjects where base was not significant)
# but the direction is Day1-cTBS-aOFC lowered the sated odor choices
# hard to explain given the hypothesis (damn it!)

model_choice_2 <- glmer(Choice ~ Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_B,family = 'binomial')
model_choice_3 <- glmer(Choice ~ ValueDiff + base + Cond + Sess + (1|SubID), 
                        data = use_trial_dat_aOFC_B,family = 'binomial')
anova(model_choice_2,model_choice_3) 
summary(model_choice_3) 
# still sig & negative effect of condition
# also strong & positive effect of session


################################################
# what if collapse across subjects so not running mixed effect models
# but still separating aOFC and pOFC subjects

glm_choice_pOFC_2 <- glm(Choice ~ Cond * ValueDiff, 
                    data = use_trial_dat_pOFC_B,family = 'binomial')
glm_choice_pOFC_1 <- glm(Choice ~ Cond + ValueDiff, 
                    data = use_trial_dat_pOFC_B,family = 'binomial')
glm_choice_pOFC_0 <- glm(Choice ~ ValueDiff, 
                    data = use_trial_dat_pOFC_B,family = 'binomial')
anova(glm_choice_pOFC_1,glm_choice_pOFC_0) 
anova(glm_choice_pOFC_2,glm_choice_pOFC_0) 
summary(glm_choice_pOFC_0)
summary(glm_choice_pOFC_1)
summary(glm_choice_pOFC_2)


glm_choice_aOFC_2 <- glm(Choice ~ Cond * ValueDiff, 
                         data = use_trial_dat_aOFC_B,family = 'binomial')
glm_choice_aOFC_1 <- glm(Choice ~ Cond + ValueDiff, 
                         data = use_trial_dat_aOFC_B,family = 'binomial')
glm_choice_aOFC_0 <- glm(Choice ~ ValueDiff, 
                         data = use_trial_dat_aOFC_B,family = 'binomial')
anova(glm_choice_aOFC_1,glm_choice_aOFC_0) 
anova(glm_choice_aOFC_2,glm_choice_aOFC_0) 
summary(glm_choice_aOFC_0)
summary(glm_choice_aOFC_1)
summary(glm_choice_aOFC_2)

# without mixed-effects models
# then everything now is not significant
# but the direction of the non-sig effects was the same

# things that are hard to explain from set B results
# (1) for aOFC subjects, day 1 cTBS is associated with less sated choices
# this seems very strong, existing with or without adding other covariates
# (2) aOFC & pOFC subjects show different patterns in terms of with other 
# covariates: session number, baseline preference before the meal (aOFC subjects only)
# both group of subjects' choices show strong correlation with 
# estimated value differences

# what could change these results?
# now the models assume linear additive effects of those variables
# the variable of "base" may not be accurate when considering the dynamics

########################################################
# set A
# pOFC subjects
use_trial_dat_pOFC_A = subset(use_trial_dat_pOFC, Set=='A')
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_A,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_A,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition
summary(model_choice_0) 
# sig effect of value difference & baseline
# but the value diff effect seems much weaker than for set B
# very strong effect of baseline (as expected)

model_choice_2 <- glmer(Choice ~ Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC_A,family = 'binomial')
model_choice_3 <- glmer(Choice ~ ValueDiff + base + Cond + Sess + (1|SubID), 
                        data = use_trial_dat_pOFC_A,family = 'binomial')
anova(model_choice_2,model_choice_3) # still no effect of condition
summary(model_choice_2) # marginal effect of Sess

# aOFC subjects
use_trial_dat_aOFC_A = subset(use_trial_dat_aOFC, Set=='A')
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_A,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_A,family = 'binomial')
anova(model_choice_1,model_choice_0) # strong effect of TMS cond & baseline
summary(model_choice_1) # weak & non-sig effect of value-diff

model_choice_2 <- glmer(Choice ~ Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC_A,family = 'binomial')
model_choice_3 <- glmer(Choice ~ ValueDiff + base + Cond + Sess + (1|SubID), 
                        data = use_trial_dat_aOFC_A,family = 'binomial')
anova(model_choice_2,model_choice_3) 
summary(model_choice_3) # not sig effect of Sess


