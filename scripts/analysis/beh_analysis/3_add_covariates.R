

# this code adds covariates of value diff and Didx to the choice wise data

rm(list = ls())
# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

#######################

load(file = file.path(processed_dir,'choice_dat_ss_w_base.RData'))

load(file = file.path(processed_dir,'choice_cue_mapping_post.RData'))
load(file = file.path(processed_dir,'conditioning_cue_mapping.RData'))

load(file = file.path(processed_dir,'SelectSate_dat.RData'))

df_w = read.csv(file = 
                  '/Users/liuq13/Library/CloudStorage/Box-Box/NODEAP_data_analysis/BehModeling_output/csv/df_w_per_sub_sess_cuepair.csv')

find_cue_id = function(left,right,sub,sess){
  # input: left, right: cue id in choice task
  # output: cuepair index during conditioning
  find_dat_left = subset(conditioning_cue_mapping,CueLeft==left & SubID==sub & Sess==sess)
  find_dat_right = subset(conditioning_cue_mapping,CueLeft==right & SubID==sub & Sess==sess)
  c(find_dat_left$CuePair,find_dat_right$CuePair)
}

find_ws = function(sub,sess,cuepairs){
  find_dat_left = subset(df_w,SubID==sub & Sess==sess & CuePair==cuepairs[1])
  find_dat_right = subset(df_w,SubID==sub & Sess==sess & CuePair==cuepairs[2])
  c(find_dat_left$w,find_dat_right$w)
}

# to prepare to add values of each option to the choice dataset (ss choices only)
use_choice_dat_ss = choice_dat_ss_w_base 
use_choice_dat_ss$ValueLeft = NA
use_choice_dat_ss$ValueRight = NA

for(i in 1:nrow(use_choice_dat_ss)){
  curr_dat = use_choice_dat_ss[i,]
  CuePairs = find_cue_id(curr_dat$CueLeft,curr_dat$CueRight,curr_dat$SubID,curr_dat$Sess)
  ws = find_ws(curr_dat$SubID,curr_dat$Sess,CuePairs)
  use_choice_dat_ss$ValueLeft[i] = ws[1]
  use_choice_dat_ss$ValueRight[i] = ws[2]
}

# remap Valueleft & ValueRight to ValueSated & ValueNonSated
# OdorLR: 1-Left W, right A; 2-Left A, right W
use_choice_dat_ss <- use_choice_dat_ss %>%
  mutate(
    ValueSated = case_when(
      Devalued == "W" & OdorLR == 1 ~ ValueLeft,
      Devalued == "W" & OdorLR == 2 ~ ValueRight,
      Devalued == "A" & OdorLR == 1 ~ ValueRight,
      Devalued == "A" & OdorLR == 2 ~ ValueLeft,
      TRUE ~ NA_real_  # Default case for unmatched conditions
    )
  ) %>%
  mutate(
    ValueNonSated = case_when(
      Devalued == "W" & OdorLR == 1 ~ ValueRight,
      Devalued == "W" & OdorLR == 2 ~ ValueLeft,
      Devalued == "A" & OdorLR == 1 ~ ValueLeft,
      Devalued == "A" & OdorLR == 2 ~ ValueRight,
      TRUE ~ NA_real_  # Default case for unmatched conditions
    )
  )

# Compute trial-wise value difference (sated - non-sated) and 
# merge session-level Didx from SelectSate_dat based on subject and condition
use_choice_dat_ss <- use_choice_dat_ss %>%
  mutate(ValueDiff = ValueSated - ValueNonSated) %>%
  left_join(
    SelectSate_dat %>% select(SubID, Cond, Didx),
    by = c("SubID", "Cond")
  )

save(use_choice_dat_ss,
     file = file.path(processed_dir,'choice_dat_ss_w_base_values.RData'))


