
rm(list = ls())
source('Setup.R')

#######################

pro_dat_dir = '/Users/liuq13/Library/CloudStorage/Box-Box/NODEAP_data_analysis/ProcessedData'

load(file = file.path(pro_dat_dir,'choice_dat_ss_w_base.RData'))
load(file = file.path(pro_dat_dir,'choice_cue_mapping_post.RData'))
load(file = file.path(pro_dat_dir,'conditioning_cue_mapping.RData'))

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


# save the full choice data set
save(use_choice_dat_ss,
     file = file.path(pro_dat_dir,'choice_dat_ss_w_base_values.RData'))


# showing for each subject
s1 = use_choice_dat_ss %>%
  mutate(Choice=factor(Choice)) %>%
  filter(complete.cases(Choice)) %>%
  ggplot() +
  geom_jitter(aes(x=ValueSated,y=ValueNonSated,
                  shape = Choice,color = Choice),size=2) +
  facet_wrap(~SubID) +
  scale_shape_manual(values = c(21, 15))


pdf(file.path(FigDir,'Choices_w_learned_values_subs.pdf'),10,10)
print(s1)
dev.off()

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

p_size = ggplot(summary_use_choice_dat_ss,
            aes(x = Rounded_ValueSated, 
                y = Rounded_ValueNonSated, 
                color = log(nsample))) +  # Circle color indicates trial number
  geom_point(shape = 16, alpha = 0.8, size = 4) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  labs(x = use_xlab, y = use_ylab) +
  scale_color_gradient2(low = "#578f9f", high = "#003744", 
                        name = expression(log(n[trial]))) +
  theme(legend.position = "right") + common_setting

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

pdf(file.path(FigDir,'Choice_w_values_agg_gradient.pdf'),5,4)
print(p_choice_1)
dev.off()

pdf(file.path(FigDir,'Choice_w_values_agg.pdf'),5,4)
print(p_choice_2)
dev.off()

pdf(file.path(FigDir,'Choice_w_values_size_agg.pdf'),5,4)
print(p_size)
dev.off()




# look at this for each subject
summary_use_choice_dat_ss_subs = 
  use_choice_dat_ss %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,SubID) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n()) 

q1_sub = ggplot(summary_use_choice_dat_ss_subs,
            aes(x = Rounded_ValueSated, 
                y = Rounded_ValueNonSated, 
                color = log(nsample))) +  # Circle color indicates trial number
  geom_point(shape = 16, alpha = 0.8, size = 4) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  facet_wrap(~SubID) +
  labs(x = expression(w[SA]), y = expression(w[NS])) +
  scale_color_gradient2(low = "#578f9f", high = "#003744", 
                        name = expression(log(n[trial]))) +
  theme(legend.position = "right") + common_setting


p1_sub=ggplot(summary_use_choice_dat_ss_subs,
          aes(x = Rounded_ValueSated, 
              y = Rounded_ValueNonSated, 
              size = log(nsample),  # Circle size indicates log sample size
              color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  facet_wrap(~SubID) +
  scale_size_continuous(range = c(2, 6)) +  # Adjust circle size range
  labs(x = expression(w[SA]), y = expression(w[NS]),
       size = expression(log(n[trial]))) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") + common_setting


pdf(file.path(FigDir,'Choices_w_learned_values_subs_replot.pdf'),12,12)
print(q1_sub)
print(p1_sub)
dev.off()


# separated by PrePost and Set
summary_use_choice_dat_ss = use_choice_dat_ss %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,PrePost,Set) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n()) 

p2=ggplot(summary_use_choice_dat_ss,
       aes(x = Rounded_ValueSated, 
           y = Rounded_ValueNonSated, 
           size = log(nsample),  # Circle size indicates log sample size
           color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  facet_wrap(~PrePost+Set) +
  labs(x = "Value of sated option",
       y = "Value non-sated option",
       size = "log (trial number)") +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P (sated odor)") +
  theme(legend.position = "right")

pdf(file.path(FigDir,'Choice_w_values_pre_post_sets.pdf'),12,4)
print(p2)
dev.off()


# focus on post-meal choices (both sets)
summary_use_choice_dat_ss_byCond_Post = use_choice_dat_ss %>%
  subset(PrePost=='Post') %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,Cond,StimLoc) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n())

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     background_y = elem_list_rect(fill = 'transparent'),
                     text_y = elem_list_text(color = use.col.conds,
                                             face = "bold",
                                             size = 16),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

p4=ggplot(summary_use_choice_dat_ss_byCond_Post, 
          aes(x = Rounded_ValueSated, 
              y = Rounded_ValueNonSated, 
              size = log(nsample),  # Circle size indicates log sample size
              color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = expression(w[SA]), y = expression(w[NS]),
       size = expression(log(n[trial]))) +
  facet_grid2(cols = vars(StimLoc), 
              rows = vars(Cond), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") +
  common_setting
  

pdf(file.path(FigDir,'Choice_w_values_TMS_conds.pdf'),9.5,12)
print(p4)
dev.off()

# look at P(SA) and value diff
summary_use_choice_dat_ss_byCond_Post_vdiff = use_choice_dat_ss %>%
  subset(PrePost=='Post') %>%
  mutate(value_diff = Rounded_ValueSated - Rounded_ValueNonSated) %>%
  mutate(value_diff = factor(value_diff)) %>%
  dplyr::group_by(value_diff,Cond,StimLoc) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n()) %>%
  mutate(value_diff = varhandle::unfactor(value_diff))

# look at Day 1 TMS effect
p5 = summary_use_choice_dat_ss_byCond_Post_vdiff %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=value_diff,
             y=MeanChoice, color = Cond, fill = Cond)) +
  geom_point(aes(size = log(nsample)),
             shape = 16, alpha = 0.6) +
  geom_smooth(method = lm, se = TRUE, 
              linetype = 2) + 
  stat_cor(method = 'pearson',show.legend = F) +
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  facet_grid2(cols = vars(StimLoc), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  labs(x = expression(w[SA]-w[NS]), y = 'P(SA)') +
  common_setting + 
  theme_classic(base_size = 16) +
  guides(color = guide_legend(title = NULL), 
         fill = guide_legend(title = NULL),  
         size = guide_legend(title = expression(log(n[trial])))) 

pdf(file.path(FigDir,'Choice_w_values_TMS_conds_Day1_vdiff.pdf'),10,4)
print(p5)
dev.off()

# look at Day 2 TMS effect
p6 = summary_use_choice_dat_ss_byCond_Post_vdiff %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=value_diff,
             y=MeanChoice, color = Cond, fill = Cond)) +
  geom_point(aes(size = log(nsample)),
             shape = 16, alpha = 0.6) +
  geom_smooth(method = lm, se = TRUE, 
              linetype = 2) + 
  stat_cor(method = 'pearson',show.legend = F) +
  facet_grid2(cols = vars(StimLoc), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = expression(w[SA]-w[NS]), y = 'P(SA)') +
  common_setting + 
  theme_classic(base_size = 16) +
  guides(color = guide_legend(title = NULL), 
         fill = guide_legend(title = NULL),  
         size = guide_legend(title = expression(log(n[trial]))))  

pdf(file.path(FigDir,'Choice_w_values_TMS_conds_Day2_vdiff.pdf'),10,4)
print(p6)
dev.off()

