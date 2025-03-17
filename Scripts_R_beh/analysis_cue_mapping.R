

# to find the odor-predicting cue id associated with each cue pair index
# from this we can map the w values indexed by cue pair index to cue id

# do this for conditioning cues, pre-meal choice cues, post-meal choice cues

rm(list = ls())
source('Setup.R')

load(file = '../ProcessedData/Conditioning.RData')

conditioning_cue_mapping = conditioning_dat %>%
  subset(OdorLR==1 & Run==1) %>%
  select(SubID,Sess,Cond,CuePair,Set,ChoiceType,CueLeft,CueRight) %>%
  arrange(SubID,Sess,CuePair)
# CueLeft should be the column with odor-predicting cue id
# because I only included trials with odor on the left

save(conditioning_cue_mapping,file = '../ProcessedData/conditioning_cue_mapping.RData')
write.xlsx(conditioning_cue_mapping,file = '../ProcessedData/conditioning_cue_mapping.xlsx')


# mapping on pre-meal choice
load(file = '../ProcessedData/choice_dat.RData')
choice_cue_mapping_pre = choice_dat %>%
  subset(PrePost=='Pre' & OdorLR==1) %>%
  select(SubID,Sess,CuePair,Set,ChoiceType,CueLeft,CueRight) %>%
  arrange(SubID,Sess,CuePair)
choice_cue_mapping_post = choice_dat %>%
  subset(PrePost=='Post' & OdorLR==1) %>%
  select(SubID,Sess,CuePair,Set,ChoiceType,CueLeft,CueRight) %>%
  arrange(SubID,Sess,CuePair)

save(choice_cue_mapping_pre,file = '../ProcessedData/choice_cue_mapping_pre.RData')
write.xlsx(choice_cue_mapping_pre,file = '../ProcessedData/choice_cue_mapping_pre.xlsx')

save(choice_cue_mapping_post,file = '../ProcessedData/choice_cue_mapping_post.RData')
write.xlsx(choice_cue_mapping_post,file = '../ProcessedData/choice_cue_mapping_post.xlsx')


