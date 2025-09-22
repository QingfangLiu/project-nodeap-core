

# this code adds covariates of value diff and Didx to the choice wise data
# then it plots choices against them

rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

#######################

load(file = file.path(processed_dir,'choice_dat_ss_w_base.RData'))

load(file = file.path(processed_dir,'choice_cue_mapping_post.RData'))
load(file = file.path(processed_dir,'conditioning_cue_mapping.RData'))

load(file = file.path(processed_dir,'SelectSate_dat.RData'))

df_w = read.csv(file = file.path(beh_model_dir,'results','df_w_per_sub_sess_cuepair.csv'))

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


# to round the learned values
target_values <- seq(0.5, 1, by = 0.05)
target_values
use_choice_dat_ss = use_choice_dat_ss %>%
  subset(SubID != 'NODEAP_17') %>%
  mutate(Rounded_ValueSated = as.numeric(sapply(ValueSated, function(x) 
    target_values[which.min(abs(target_values - x))]))) %>%
  mutate(Rounded_ValueNonSated = as.numeric(lapply(ValueNonSated, function(x) 
    target_values[which.min(abs(target_values - x))]))) %>%
  filter(complete.cases(Choice))

# get a summary across subjects and conditions
summary_use_choice_dat_ss = use_choice_dat_ss %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),
          nsample=n()) 

summary_use_choice_dat_ss_vdiff = use_choice_dat_ss %>%
  mutate(value_diff = Rounded_ValueSated - Rounded_ValueNonSated) %>%
  mutate(value_diff = factor(value_diff)) %>%
  dplyr::group_by(value_diff) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),
          nsample=n()) %>%
  mutate(value_diff = varhandle::unfactor(value_diff))

common_setting = 
  theme(axis.text.x = element_text(colour = "black"), # axis color
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5), # center title
        panel.background = element_rect(fill = "gray90"),  # Lighten the panel background
        panel.grid = element_line(color = "gray85"), # Adjust grid line colors
        #plot.margin = margin(t = 10, r = 5.5, b = 5.5, l = 5.5, unit = 'pt'),
        text = element_text(size=16)) # text size & font

use_xlab = expression(w[SA])
use_ylab = expression(w[NS])

p_choice_1 = ggplot(summary_use_choice_dat_ss,
                    aes(x = Rounded_ValueSated, 
                        y = Rounded_ValueNonSated, 
                        size = log(nsample),  # Circle size indicates log sample size
                        color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = expression(w[SA]), y = expression(w[NS]),
       size = expression(log(n[trial]))) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") + common_setting

p_choice_2 = summary_use_choice_dat_ss_vdiff %>%
  ggplot(aes(x=value_diff,
             y=MeanChoice)) +
  geom_point(aes(size = log(nsample)),
             shape = 16, alpha = 0.6) +
  geom_smooth(method = lm, se = TRUE, color = 'black',
              linetype = 2) + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = expression(w[SA]-w[NS]), 
       size = expression(log(n[trial])),
       y = 'P(SA)') +
  common_setting + 
  theme_classic(base_size = 16)

cor.test(summary_use_choice_dat_ss_vdiff$value_diff, 
         summary_use_choice_dat_ss_vdiff$MeanChoice)

pdf(file.path(FigPaperDir,'Choice_w_values_agg_gradient.pdf'),10,4)
ggarrange(p_choice_1,p_choice_2)
dev.off()

