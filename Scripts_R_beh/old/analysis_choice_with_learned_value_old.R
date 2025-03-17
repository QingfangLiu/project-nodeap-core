
rm(list = ls())
source('Setup.R')

#######################

load(file = '../ProcessedData/choice_dat.RData')
load(file = '../ProcessedData/choice_cue_mapping_post.RData')
load(file = '../ProcessedData/conditioning_cue_mapping.RData')

df_w = read.csv(file = '../BehModeling/df_w_per_sub_sess_cuepair.csv')

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

# focus on sweet-savory type of choices for now
use_choice_dat_ss = choice_dat %>% 
  subset(ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) # Choice==1 means selecting sated odor
  
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
     file = '../ProcessedData/choice_dat.RData')


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

p_choice = ggplot(summary_use_choice_dat_ss,
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

pdf(file.path(FigDir,'Choice_w_values_agg.pdf'),10,4)
ggarrange(p_choice,p_size,ncol = 2)
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


summary_use_choice_dat_ss_byCond = use_choice_dat_ss %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,Cond,StimLoc) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n())

p3=ggplot(summary_use_choice_dat_ss_byCond, 
       aes(x = Rounded_ValueSated, 
           y = Rounded_ValueNonSated, 
           size = log(nsample),  # Circle size indicates log sample size
           color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = "Value of sated option",
       y = "Value non-sated option",
       title = 'All choices',
       size = "log (trial number)") +
  facet_wrap(~StimLoc+Cond) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P (sated odor)") +
  theme(legend.position = "right")


summary_use_choice_dat_ss_byCond_Post = use_choice_dat_ss %>%
  subset(PrePost=='Post') %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,Cond,StimLoc) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n())

p4=ggplot(summary_use_choice_dat_ss_byCond_Post, 
          aes(x = Rounded_ValueSated, 
              y = Rounded_ValueNonSated, 
              size = log(nsample),  # Circle size indicates log sample size
              color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = "Value of sated option",
       y = "Value non-sated option",
       title = 'Post-meal choices',
       size = "log (trial number)") +
  facet_wrap(~StimLoc+Cond) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P (sated odor)") +
  theme(legend.position = "right")


summary_use_choice_dat_ss_byCond_PostB = use_choice_dat_ss %>%
  subset(PrePost=='Post' & Set=='B') %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,Cond,StimLoc) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n())

p5=ggplot(summary_use_choice_dat_ss_byCond_PostB, 
          aes(x = Rounded_ValueSated, 
              y = Rounded_ValueNonSated, 
              size = log(nsample),  # Circle size indicates log sample size
              color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = "Value of sated option",
       y = "Value non-sated option",
       title = 'Post-meal choices (B)',
       size = "log (trial number)") +
  facet_wrap(~StimLoc+Cond) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P (sated odor)") +
  theme(legend.position = "right")


pdf(file.path(FigDir,'Choice_w_values_TMS_conds.pdf'),12,8)
print(p3)
print(p4)
print(p5)
dev.off()

# if concerned with Day 2 TMS effect

xx_Day2 = summary_use_choice_dat_ss_byCond_Post %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x = Rounded_ValueSated, 
              y = Rounded_ValueNonSated, 
              size = log(nsample),  # Circle size indicates log sample size
              color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = use_xlab, y = use_ylab,
       title = 'Day 2 TMS effect',
       size = expression(log(n[trial]))) +
  facet_wrap(~StimLoc+Cond) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") + common_setting

summary_use_choice_dat_ss_byCond_Post_bySet = use_choice_dat_ss %>%
  subset(PrePost=='Post') %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated,Cond,StimLoc,Set) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),nsample=n())

xx_Day1 = summary_use_choice_dat_ss_byCond_Post_bySet %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x = Rounded_ValueSated, 
             y = Rounded_ValueNonSated, 
             size = log(nsample),  # Circle size indicates log sample size
             color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = use_xlab, y = use_ylab,
       title = 'Day 1 TMS effect',
       size = expression(log(n[trial]))) +
  facet_wrap(~Set+StimLoc+Cond,nrow = 2) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") + common_setting


pdf(file.path(FigDir,'Choice_w_values_TMS_conds_Day2.pdf'),8,8)
print(xx_Day2)
dev.off()

pdf(file.path(FigDir,'Choice_w_values_TMS_conds_Day1.pdf'),14,8)
print(xx_Day1)
dev.off()

