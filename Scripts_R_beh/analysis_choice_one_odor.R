
# load choice data (in full)
# containing NODEAP_17 (may be at the end of the dataset)

rm(list = ls())
source('Setup.R')
load(file = '../ProcessedData/choice_dat.RData')

##############################
# only look at choice with one odor and one air

# ChoiceType: 1,2 (1: sweet, 2:savory)
# ChosenOdor: 0,1 (1 for selecting odor)

choice_dat_one = choice_dat %>%
  subset(ChoiceType<3) %>%
  mutate(ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(OdorType=ifelse(ChoiceType==Devalued,'sated','non-sated'))

# focus on Day 2 TMS effect
summary_choice_one = choice_dat_one %>%
  group_by(SubID,StimLoc,Cond,PrePost,OdorType) %>%
  reframe(Choice=mean(ChosenOdor,na.rm = T))

# create a session-wise measure of odor preference called 'pref':
# preference of sated odor over non-sated odors (before eating meal)
pre_choice_df = summary_choice_one %>%
  subset(PrePost=='Pre') %>%
  group_by(SubID,Cond,OdorType) %>%
  reframe(PreChoice=Choice) 
pre_choice_df_sated = subset(pre_choice_df,OdorType=='sated')$PreChoice
pre_choice_df_non_sated = subset(pre_choice_df,OdorType=='non-sated')$PreChoice
pre_choice_df_sated_pref = subset(pre_choice_df,OdorType=='sated')
pre_choice_df_sated_pref$pref = pre_choice_df_sated - pre_choice_df_non_sated
pre_choice_df_sated_pref$OdorType = NULL
pre_choice_df_sated_pref$PreChoice = NULL

save(pre_choice_df_sated_pref,
     file = '../ProcessedData/pre_choice_df_sated_pref_one_odor.RData')


c0=summary_choice_one %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=OdorType,y=Choice)) +
  geom_boxplot(aes(linetype = PrePost),outlier.alpha = 0,
               alpha=0.4) +
  geom_jitter(aes(shape = PrePost), 
              position = position_jitterdodge(
                dodge.width = 0.6, 
                jitter.height = 0.05, # added jitter height to reduce overlapping
                jitter.width = 0.3), 
              size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") + 
  scale_shape_manual(values = c(17, 1)) +
  labs(x = NULL, title = NULL, y = "Choice of odor (vs. air)") + common +
  theme(legend.position = 'none')

pdf(file.path(FigDir,'Day2_ChoiceOneOdor.pdf'),4,4)
print(c0)
dev.off()


c1=summary_choice_one %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=OdorType,y=Choice)) +
  geom_boxplot(aes(linetype = PrePost),outlier.alpha = 0,
    alpha=0.4) +
  geom_jitter(aes(shape = PrePost), 
              position = position_jitterdodge(
                dodge.width = 0.6, 
                jitter.height = 0.05, # added jitter height to reduce overlapping
                jitter.width = 0.3), 
              size = 1.5, alpha = 0.8) +
  scale_shape_manual(values = c(17, 1)) +
  facet_wrap(~Cond,scales = 'free_y') +
  labs(x = NULL, title = NULL, y = "Choice of odor (vs. air)") + common +
  theme(legend.position = 'none')

pdf(file.path(FigDir,'Day2TMS_ChoiceOneOdor.pdf'),8,4)
print(c1)
dev.off()


summary_choice_one_pre = subset(summary_choice_one,PrePost=='Pre')$Choice
summary_choice_one_post = subset(summary_choice_one,PrePost=='Post')$Choice
summary_choice_one_change = subset(summary_choice_one,PrePost=='Pre')
summary_choice_one_change$Change = summary_choice_one_post - 
                                   summary_choice_one_pre
summary_choice_one_change$Choice = NULL
summary_choice_one_change$PrePost = NULL


p_values <- data.frame(
  Cond = c('sham-sham','sham-cTBS'),
  p_value = c(0.018, 0.91),     
  y_start = c(1.05, 1.05),        
  y_end = c(1.15, 1.15),
  y_low = c(1, 1)
)

c2=summary_choice_one_change %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=OdorType,y=Change,fill = Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.2, 
            color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 1, alpha = 0.8) +
  facet_wrap(~ Cond,scales = 'free_y') +
  #facet_wrap(~ StimLoc + Cond, scales = 'free_y',nrow = 1) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
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
    aes(x = 1.5,y = y_end,
        label = paste0("p = ", format(p_value, scientific = FALSE))), 
    inherit.aes = FALSE, size = 3.5) +
  labs(x = NULL, title = NULL,
       y = "Updated choice of odor (vs. air)") + common +
  coord_cartesian(ylim = c(-1, 1.3)) +
  theme(legend.position = "none") 

pdf(file.path(FigDir,'Day2TMS_ChoiceOneOdor_changes.pdf'),7,4)
print(c2)
dev.off()


summary_choice_one_change %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=OdorType,y=Change)) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~ StimLoc + Cond) +
  scale_color_manual(values = use.col.conds) +
  labs(x = "", y = "Choosing Odor (vs. air)") + common +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p")


summary_choice_one_sated = subset(summary_choice_one_change,
                                  OdorType=='sated')$Change
summary_choice_one_nonsated = subset(summary_choice_one_change,
                                     OdorType=='non-sated')$Change
summary_choice_one_comp = subset(summary_choice_one_change,OdorType=='sated')
summary_choice_one_comp$diff = summary_choice_one_sated - 
  summary_choice_one_nonsated
summary_choice_one_comp$OdorType = NULL
summary_choice_one_comp$Change = NULL


muh_grob <- grid::rectGrob(
  x=1:2, y=0, gp=gpar(
    color='black', fill=use.col.ap.ofc, alpha=1))

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

c3=summary_choice_one_comp %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=diff,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.2, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',
              strip = strip,
              axes = 'all') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, title = NULL,
       y = "Preference of sated odor\n (post - pre-meal)") + common +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     method.args = list(alternative='g'),
                     label = "p")

pdf(file.path(FigPaperDir,'Day2TMS_ChoiceOneOdor_Pref_change.pdf'),7,4)
print(c3)
dev.off()

# similar plot to look at Day 1 TMS effect
c4 = summary_choice_one_comp %>%
  subset(SubID != 'NODEAP_17') %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=Cond,y=diff,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.2, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',
              strip = strip,
              axes = 'all') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, title = NULL,
       y = "Preference of sated odor\n (post - pre-meal)") + common +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     method.args = list(alternative='g'),
                     label = "p")

pdf(file.path(FigPaperDir,'Day1TMS_ChoiceOneOdor_Pref_change.pdf'),7,4)
print(c4)
dev.off()


#############################

# use the overall baseline measure to do some analysis
load(file = '../ProcessedData/pre_choice_df_sated_pref_overall.RData')

# change the way of plotting
# odor preference before and after the meal
choice_df_sated = subset(summary_choice_one,OdorType=='sated')$Choice
choice_df_non_sated = subset(summary_choice_one,OdorType=='non-sated')$Choice
choice_df_sated_pref = subset(summary_choice_one,OdorType=='sated')
choice_df_sated_pref$pref = choice_df_sated - choice_df_non_sated
choice_df_sated_pref$OdorType = NULL
choice_df_sated_pref$Choice = NULL


t1=choice_df_sated_pref %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=PrePost,y=pref)) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.2, 
            color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0) +
  facet_wrap(~ Cond,scales = 'free_y') +
  scale_color_manual(values = use.col.conds) +
  labs(x = "", y = "Pref of sated over non-sated (vs. air)") + 
  common +
  coord_cartesian(ylim = c(-1, 1)) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p")

pdf(file.path(FigDir,'Day2TMS_ChoiceOneOdor_Pref.pdf'),6,4)
print(t1)
dev.off()

choice_df_sated_pref_pre = subset(choice_df_sated_pref,PrePost=='Pre')$pref
choice_df_sated_pref_post = subset(choice_df_sated_pref,PrePost=='Post')$pref
choice_df_sated_pref_change = subset(choice_df_sated_pref,PrePost=='Pre')
choice_df_sated_pref_change$change = choice_df_sated_pref_post - 
  choice_df_sated_pref_pre
choice_df_sated_pref_change$change_revised = choice_df_sated_pref_post - 
  (pre_choice_df$base_overall)*2 + 1 # rescale this baseline measure
choice_df_sated_pref_change$PrePost = NULL
choice_df_sated_pref_change$pref = NULL


p_values <- data.frame(
  StimLoc = c("aOFC", "pOFC"),  
  p_value = c(0.43, 0.06),     
  y_start = c(1.05, 1.05),        
  y_end = c(1.15, 1.15),
  y_low = c(1, 1)
)

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
t2=choice_df_sated_pref_change %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=change,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',
              strip = strip,
              axes = 'all') +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  coord_cartesian(ylim = c(-1.2, 1.15)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
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
    aes(x = 1.5,y = y_end,
      label = paste0("p = ", format(p_value, scientific = FALSE))), 
      inherit.aes = FALSE, size = 3.5) +
  labs(x = NULL, y = "Preference of sated odor\n (post - pre-meal)", 
       title = NULL) + common +
  theme(legend.position = "none") 

pdf(file.path(FigDir,'Day2TMS_ChoiceOneOdor_Pref_change1.pdf'),7,4)
print(t2)
dev.off()

#############################################################
# go to trial level
# insert this base value to the whole choice dat one df
#############################################################

choice_dat_one_w_base = choice_dat_one %>%
  mutate(base=NA)
for(i in 1:nrow(choice_dat_one_w_base)){
  choice_dat_one_w_base$base[i] = pre_choice_df$base_one_odor[
    pre_choice_df$SubID==choice_dat_one_w_base$SubID[i] &
      pre_choice_df$Cond==choice_dat_one_w_base$Cond[i]]
}


# put aOFC and pOFC together 
use.dat = subset(choice_dat_one_w_base, PrePost=='Post' &
                   StimLoc == 'pOFC' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_0 <- glmer(ChosenOdor ~ base + OdorType + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(ChosenOdor ~ base + OdorType + Cond + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_2 <- glmer(ChosenOdor ~ base + OdorType * Cond + (1|SubID), 
                        data = use.dat,family = 'binomial')

summary(model_choice_0)
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)
summary(model_choice_1)
summary(model_choice_2)



