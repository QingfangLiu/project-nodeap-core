
# this analyzes choices with only one odor
# w/ considering value estimates

rm(list = ls())
source('Setup.R')

#######################################

# load choice data (one odor vs clean air)
# already included basedline pref from pre-meal
load(file = '../ProcessedData/choice_dat_one_w_base_values.RData')

choice_dat_one_odor = choice_dat_one_w_base %>%
  subset(SubID != 'NODEAP_17') %>% # this has weird thing when trying to find cue id
  subset(PrePost=='Post')

# load the mapping data from cue to ?
load(file = '../ProcessedData/conditioning_cue_mapping.RData')
df_w = read.csv(file = '../BehModeling/df_w_per_sub_sess_cuepair.csv')

find_cue_id = function(odorcue,sub,sess){
  # input: odorcue: cue id in choice task
  # output: cue index during conditioning
  find_dat = subset(conditioning_cue_mapping,CueLeft==odorcue & SubID==sub & Sess==sess)
  return(find_dat$CuePair)
}

find_ws = function(sub,sess,Cue_id_conditioning){
  find_dat = subset(df_w,SubID==sub & Sess==sess & CuePair==Cue_id_conditioning)
  find_dat$w
}

choice_dat_one_odor$ValueOdorCue = NA

for(i in 1:nrow(choice_dat_one_odor)){
  curr_dat = choice_dat_one_odor[i,]
  odorcue = ifelse(curr_dat$OdorLR==1,curr_dat$CueLeft,curr_dat$CueRight)
  Cue_id_conditioning = find_cue_id(odorcue,curr_dat$SubID,curr_dat$Sess)
  w = find_ws(curr_dat$SubID,curr_dat$Sess,Cue_id_conditioning)
  choice_dat_one_odor$ValueOdorCue[i] = w
}

# look for day 2 TMS effect
use.dat = subset(choice_dat_one_odor, 
                 #OdorType == 'sated' &
                   StimLoc == 'aOFC' & 
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_0 <- glmer(ChosenOdor ~ ValueOdorCue + base + OdorType + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(ChosenOdor ~ ValueOdorCue + base + OdorType + Cond + (1|SubID), 
                        data = use.dat,family = 'binomial')
summary(model_choice_0)
anova(model_choice_0,model_choice_1)
summary(model_choice_1)


# look for day 1 TMS effect
use.dat = subset(choice_dat_one_odor, 
                   #OdorType == 'sated' &
                   StimLoc == 'pOFC' & 
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_0 <- glmer(ChosenOdor ~ ValueOdorCue + base + OdorType + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(ChosenOdor ~ ValueOdorCue + base + OdorType + Cond + (1|SubID), 
                        data = use.dat,family = 'binomial')

summary(model_choice_0)
anova(model_choice_0,model_choice_1)
summary(model_choice_1)


