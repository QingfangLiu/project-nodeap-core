

# this looks at sweet-savory choices
# focus on Day 2 TMS effect

rm(list = ls())
source('Setup.R')

# load the full choice data set

load(file = file.path(pro_dat_dir,'choice_dat.RData'))

# only look at choice with one sweet and one savory (denoted as '_ss')
# perhaps where the effect is the strongest
choice_dat_ss = choice_dat %>%
  subset(ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) 
  # Choice==1 means selecting sated odor

# summary of the choice_ss data
summary_choice_ss = choice_dat_ss %>%
  group_by(SubID,StimLoc,StimOrder,StimOrder_day1,StimOrder_day2,
           Cond,Sess,PrePost,Set) %>%
  reframe(Choice=mean(Choice,na.rm = T),
        rt=mean(rt,na.rm = T)) 

###############################################
# improved from above by considering session-wise odor bias
# using pre(A) as baseline
# then get difference from pre to post
###############################################

pre_choice_A = subset(summary_choice_ss,
                      Set=='A' & PrePost=='Pre')$Choice
post_choice_A = subset(summary_choice_ss,
                       Set=='A' & PrePost=='Post')$Choice
post_choice_B = subset(summary_choice_ss,
                       Set=='B' & PrePost=='Post')$Choice
conds = subset(summary_choice_ss,
               Set=='A' & PrePost=='Pre')$Cond
# check if there is pre-post correlations of choices
df_pre_post_summary = data.frame(
  'conds'=conds,
  'pre_A'=pre_choice_A,
  'post_A'=post_choice_A,
  'post_B'=post_choice_B
)

lab1 = 'pre-meal (set A)'
lab2 = 'post-meal (set A)'
lab3 = 'post-meal (set B)'

pcorr1=ggplot(df_pre_post_summary,
              aes(x=pre_A,y=post_A)) +
  geom_point(aes(color=conds),size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = T, color = 'black') +
  stat_cor(method = 'pearson',size = 5) +
  labs(x=lab1, y=lab2,title = NULL) + common +
  theme(legend.position = 'none')

  theme(legend.text = element_text(size = 14), 
        legend.key.size = unit(1.5, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.spacing.y = unit(0.05, "cm"))

pcorr2=ggplot(df_pre_post_summary,
              aes(x=pre_A,y=post_B)) +
  geom_point(aes(color=conds),size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = T, color = 'black') +
  stat_cor(method = 'pearson',size = 5) +
  labs(x=lab1, y=lab3,title = NULL) + common +
  theme(legend.position = 'none')

pcorr3=ggplot(df_pre_post_summary,
              aes(x=post_A,y=post_B)) +
  geom_point(aes(color=conds),size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = T, color = 'black') +
  stat_cor(method = 'pearson',size = 5) +
  labs(x=lab2, y=lab3,title = NULL) + common +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir,'Choice_corr_paper.pdf'),8,8)
ggarrange(pcorr1,pcorr2,pcorr3,ncol = 2,nrow = 2)
dev.off()

# session-wise measures
avg_changeAB = (post_choice_A + post_choice_B)/2 - pre_choice_A
avg_changeB = post_choice_B - pre_choice_A
avg_changeA = post_choice_A - pre_choice_A
avg_BoverA = post_choice_B - post_choice_A

# create a corrected summary dataset
# containing choice update for each session (from pre to post)
summary_choice_corrected = summary_choice_ss %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  select(SubID,StimLoc,StimOrder,Cond,Sess,
         StimOrder_day1,StimOrder_day2) %>%
  mutate(ChoiceChangeAB=avg_changeAB, # both A and B
         ChoiceChangeA=avg_changeA, # only A
         ChoiceChangeB=avg_changeB, # only B
         ChoiceBoverA=avg_BoverA) # B over A

#######################################################
###  stat testing on Day 2 TMS effect
# test using Post(A & B) - Pre, summarized choice data
# pOFC
# overall devaluation effect regardless of conditions: yes
# (later realizing this also contains cTBS-sham condition)
value_Y1Z1 = subset(summary_choice_corrected,
                    StimLoc=='pOFC')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.0023

value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                    StimLoc=='pOFC')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='pOFC')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? No
wilcox.test(value_Z1,alternative = 'l') # lower than zero? Yes, p=0.026
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s

# stat testing on aOFC
# overall devaluation effect regardless of conditions: yes
value_Y1Z1 = subset(summary_choice_corrected,
                    StimLoc=='aOFC')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.024

value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                    StimLoc=='aOFC')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='aOFC')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? n.s.
wilcox.test(value_Z1,alternative = 'l') # lower than zero? n.s.
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s

# an overall comparison across TMS conditions
# for each of aOFC and pOFC subjects
muh_grob <- rectGrob(
  x=1:2, y=0, gp=gpar(color='black', fill=use.col.ap.ofc, alpha=1))
s1=summary_choice_corrected %>%
  ggplot(aes(x=StimLoc,y=ChoiceChangeAB)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0, alpha = 0.4) +
  scale_color_manual(values = use.col.conds) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  annotate("text", x = 1, y = 0.45, label = "*", size = 4, color = "black") + # aOFC significance
  annotate("text", x = 2, y = 0.75, label = "***", size = 4, color = "black") +
  labs(x = NULL, y = "Choice of sated odor \n (post - pre-meal)",
       title = NULL) + 
  common +
  coord_cartesian(clip='off') +
  theme(axis.text.x = element_text(margin=margin(t=10),color = 'white')) +
  annotation_custom(
     grob=muh_grob, xmin = 0, xmax = 1, ymin = -1.08, ymax= - 0.88
  ) +
  theme(legend.position = 'inside',
    legend.position.inside = c(0.3,0.96),
      legend.key.size = unit(0.28,'cm')) 

pdf(file.path(FigPaperDir,'Day2_overall_ChoiceSatedOdor.pdf'),4,4)
print(s1)
dev.off()



########### plot focus on Day 2 TMS ################

p_values <- data.frame(
  StimLoc = c("aOFC", "pOFC"),  
  p_value = c(0.655, 0.00034),     
  y_start = c(0.74, 0.74),        
  y_end = c(0.81, 0.81),
  y_low = c(0.72, 0.72)
)

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
c1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  facet_wrap2(~StimLoc,scales = 'free',
              strip = strip,
              axes = 'all') +
  geom_segment(
    data = p_values,aes(x = 1, xend = 2, y = y_start, yend = y_start),
    inherit.aes = FALSE, color = "black",linewidth = 0.5) +
  geom_segment(
    data = p_values,aes(x = 1, xend = 1, y = y_low, yend = y_start),
    inherit.aes = FALSE, color = "black",linewidth = 0.5) +
  geom_segment(
    data = p_values,aes(x = 2, xend = 2, y = y_low, yend = y_start),
    inherit.aes = FALSE, color = "black",linewidth = 0.5) +
  geom_text(
    data = p_values,
    aes(
      x = 1.5,         # Position the label between the two groups
      y = y_end,     # Use the pre-calculated y position
      label = paste0("p = ", p_value)
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  annotate("text",x=1.5,y=0.9,label='Mixed-effects logistic model:') +
  scale_fill_manual(values = use.col.conds) +
  scale_color_manual(values = use.col.conds) +
  coord_cartesian(ylim = c(-0.8,0.95)) +
  labs(x = NULL, title = NULL,
       y = "Choice of sated odor\n (post - pre-meal)") + common +
  theme(legend.position = "none")

custom_labeller <- labeller(
  StimOrder_day2 = label_value
)

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
# based on c1
# separated by diff stim orders of Day 2 TMS
c2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.9) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(StimOrder_day2), 
              scales = 'fixed',
              axes = 'all',
              strip = strip,
              labeller = custom_labeller) +
  scale_fill_manual(values = use.col.conds) +
  scale_color_manual(values = use.col.conds) +
  labs(x = NULL, title = NULL, y = "Choice of sated odor\n (post - pre-meal)") + 
  common + theme(legend.position = "none")

pdf(file.path(FigPaperDir,'Day2_TMS_ChoiceSatedOdor_Changes.pdf'),7,4)
print(c1)
dev.off()

pdf(file.path(FigPaperDir,'Day2_TMS_ChoiceSatedOdor_Changes_by_order.pdf'),11,8)
print(c2)
dev.off()


## similarly plot Day 1 TMS effect
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



#######################################################
# testing Day 2 TMS effect at trial level
#######################################################

# go back to choice_dat_ss to do the testing at trial level
# first calculate a session-wise pre-meal bias for each sub
pre_choice_df = summary_choice_ss %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  group_by(SubID,Cond) %>%
  reframe(base=Choice) 
save(pre_choice_df,
     file = '../ProcessedData/pre_choice_df_sated_pref_ss.RData')

###############################

# session-wise summary data, baseline corrected
save(summary_choice_corrected,file = '../ProcessedData/Summary_Choice_corrected_dat.RData')
write.xlsx(summary_choice_corrected,file = '../ProcessedData/Summary_Choice_corrected_dat.xlsx')

# savory-sweet choices only, w/t session-wise baseline
save(choice_dat_ss_w_base,file = '../ProcessedData/choice_dat_ss_w_base.RData')




