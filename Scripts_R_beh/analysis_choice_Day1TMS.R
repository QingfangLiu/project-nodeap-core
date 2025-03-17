

# this code continues with 'analysis_choice.R' code
# with a focus on Day 1 TMS effect on choice

rm(list = ls())
source('Setup.R')
load(file = '../ProcessedData/Summary_Choice_corrected_dat.RData')
load(file = '../ProcessedData/choice_dat.RData')
load(file = '../ProcessedData/choice_dat_ss_w_base.RData')

#######################################################
# stat testing on Day 1 TMS effect
# test using summarized choice data
#######################################################

# need to exclude NODEAP_17 who has missing sessions 
# (this subject has sham-sham sham-cTBS that can contribute one data point
# to Day 2 TMS effect though)

# pOFC subjects
use_summ_dat = subset(summary_choice_corrected,
                      (!SubID=='NODEAP_17') &
                        Cond %in% c('sham-sham','cTBS-sham') & 
                        StimLoc=='pOFC')

# on set A
value_XZ_A = use_summ_dat$ChoiceChangeA
wilcox.test(value_XZ_A,alternative = 'l') # lower than zero? Yes, p-value = 0.0001378
value_X_A = subset(use_summ_dat, Cond=='cTBS-sham')$ChoiceChangeA
value_Z_A = subset(use_summ_dat, Cond=='sham-sham')$ChoiceChangeA
wilcox.test(value_X_A,alternative = 'l') # lower than zero? Yes, p=0.004431
wilcox.test(value_Z_A,alternative = 'l') # lower than zero? Yes, p=0.006773
wilcox.test(value_X_A,value_Z_A,paired = T,alternative = 't') # No, p = 0.8201

# on set B
value_XZ_B = use_summ_dat$ChoiceChangeB
wilcox.test(value_XZ_B,alternative = 'l') # lower than zero? p = 0.03584
value_X_B = subset(use_summ_dat, Cond=='cTBS-sham')$ChoiceChangeB
value_Z_B = subset(use_summ_dat, Cond=='sham-sham')$ChoiceChangeB
wilcox.test(value_X_B,alternative = 'l') # lower than zero? No, p=0.2003
wilcox.test(value_Z_B,alternative = 'l') # lower than zero? No, p=0.06759
wilcox.test(value_X_B,value_Z_B,paired = T,alternative = 't') # No, p=0.5531


#######################################################

# collapsing across set A and B
strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
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

# note: this one plots the summary stats, but we want the variance
# explained by TMS, after accounting for variance explained by 
# other covariates (value-diff mainly)

#######################################################
# stat testing Day 1 TMS effect at trial level
# focus on pOFC stim
#######################################################

survey_dat = read.xlsx('../ProcessedData/Survey_uncomf_strong.xlsx') %>%
  subset(TMSSess %in% c(1,3,5)) %>%
  subset(!SubID=='NODEAP_17') %>%
  mutate(rating = (uncom + strong)/2) %>%
  mutate(Sess=rep(1:3,47))

use_trial_dat = subset(choice_dat_ss_w_base, 
                       (!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         StimLoc=='pOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.)) %>%  # remove all NA trials
  mutate(rating = NA)

# add TMS rating variable to the trial-wise data
for(i in 1:nrow(use_trial_dat)){
  use_trial_dat$rating[i] = survey_dat$rating[survey_dat$SubID==use_trial_dat$SubID[i] &
                                                survey_dat$Sess==use_trial_dat$Sess[i]]
}

# on set A
use.dat = subset(use_trial_dat, Set=='A')
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition

model_choice_2 <- glmer(Choice ~ Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3) # no effect of condition
summary(model_choice_3)

# after including TMS ratings
model_choice_4 <- glmer(Choice ~ rating + Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_5 <- glmer(Choice ~ rating + base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_4,model_choice_5) # no effect of condition
summary(model_choice_5)


# on Set B
use.dat = subset(use_trial_dat, Set=='B')
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p=0.03048, adding Cond helps
summary(model_choice_1)

model_choice_2 <- glmer(Choice ~ Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3) # p=0.0308, adding Cond helps even w/ Sess
summary(model_choice_3)

# if the model contains TMS rating
model_choice_4 <- glmer(Choice ~ rating + Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_5 <- glmer(Choice ~ rating + base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_4,model_choice_5) # p=0.01009
summary(model_choice_5)

ranef(model_choice_3) # diff b/t individual & model intercept
coef(model_choice_3) # eq of each individual
dotplot(ranef(model_choice_3))

plot(model_choice_3)
hist(resid(model_choice_3))
qqnorm(resid(model_choice_3))
qqline(resid(model_choice_3), col = "red", lwd = 3)

res = resid(model_choice_3)

new_data1 <- use.dat %>%
  group_by(SubID,Cond) %>%
  reframe(obs_choice = mean(Choice,na.rm=T))

new_data <- use.dat %>%
  group_by(SubID) %>%
  reframe(Cond = unique(Cond),  # Unique values of Cond for each subject
          base = mean(base),    # Set base to its mean
          Sess = mean(Sess))    # Set Sess to its mean
new_data$obs_choice = new_data1$obs_choice
new_data$StimLoc = 'pOFC'
new_data$predicted_prob <- predict(model_choice_3, newdata = new_data, 
                                   type = "response", re.form = NULL)
ggplot(new_data,aes(x = Cond, y = predicted_prob)) +
  geom_line(aes(group = SubID)) +
  geom_point(data=new_data,aes(x = Cond, y = obs_choice),alpha=0.1) +
  labs(title = "Effect of cTBS",
       x = "",
       y = "Predicted P (Choosing sated odor)") +
  common +
  theme(legend.position = "none")  # Remove legend if there are too many subjects

x_pOFC=ggplot(new_data,aes(x=obs_choice,y=predicted_prob)) +
  geom_point(aes(color=Cond),size=2.5) +
  geom_line(aes(group=SubID)) +
  geom_smooth(method = 'lm',color='black',lty=3) +
  stat_cor(method = 'pearson') +
  scale_color_manual(values = use.col.conds) +
  common +
  labs(x='Experimental set B choices',title = NULL,
       y='Predicted set B choices') +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.7,0.15)) +
  facet_grid(. ~ StimLoc) +
  theme(strip.background = element_rect(fill=use.col.ap.ofc[2]),
        strip.text = element_text(size=15, colour="white", face = 'bold'))



# include both set A and B in trial-wise model
model_choice_0 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Set + Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_2 <- glmer(Choice ~ Set * Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
anova(model_choice_0,model_choice_1) # adding Set helps model fit
summary(model_choice_1) # set B has higher sated choices
anova(model_choice_2,model_choice_1) 
summary(model_choice_2)
# interaction b/t Set and Day1-TMS Cond not detected, p=0.447

model_choice_0 <- glmer(Choice ~ rating + Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ rating + Set + Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_2 <- glmer(Choice ~ rating + Set * Cond + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
anova(model_choice_0,model_choice_1) # adding Set helps model fit
summary(model_choice_1) # set B has higher sated choices
anova(model_choice_2,model_choice_1) 
summary(model_choice_2) # still no interaction detected

###############################
# on aOFC
# it affects both set A and set B
use_trial_dat = subset(choice_dat_ss_w_base, 
                       PrePost=='Post' &
                         StimLoc=='aOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.))

# interaction b/t sets and TMS
model_choice_0 <- glmer(Choice ~ Sess + base + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Set + base + Cond + Sess + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
model_choice_2 <- glmer(Choice ~ Set * Cond + base + Sess + (1|SubID), 
                        data = use_trial_dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0) # p=2.288e-05
summary(model_choice_2)

# on set A
use.dat = subset(use_trial_dat, Set=='A')
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p = 0.00085

model_choice_2 <- glmer(Choice ~ Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3) # p = 0.001852
summary(model_choice_3)

# on Set B
use.dat = subset(use_trial_dat, Set=='B')
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p=0.1198
summary(model_choice_1)

model_choice_2 <- glmer(Choice ~ Sess + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ base + Cond + Sess + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3) # p=0.0351, adding Cond helps even w/ Sess
summary(model_choice_3)
ranef(model_choice_3) # diff b/t individual & model intercept
coef(model_choice_3) # eq of each individual
dotplot(ranef(model_choice_3))

plot(model_choice_3)
hist(resid(model_choice_3))
qqnorm(resid(model_choice_3))
qqline(resid(model_choice_3), col = "red", lwd = 3)

res = resid(model_choice_3)

new_data1 <- use.dat %>%
  group_by(SubID,Cond) %>%
  reframe(obs_choice = mean(Choice,na.rm=T))

new_data <- use.dat %>%
  group_by(SubID) %>%
  reframe(Cond = unique(Cond),  # Unique values of Cond for each subject
          base = mean(base),    # Set base to its mean
          Sess = mean(Sess))    # Set Sess to its mean
new_data$obs_choice = new_data1$obs_choice
new_data$StimLoc = 'aOFC'
new_data$predicted_prob <- predict(model_choice_3, newdata = new_data, 
                                   type = "response", re.form = NULL)
ggplot(new_data,aes(x = Cond, y = predicted_prob)) +
  geom_line(aes(group = SubID)) +
  geom_point(data=new_data,aes(x = Cond, y = obs_choice),alpha=0.1) +
  labs(title = "Effect of cTBS",
       x = "",
       y = "Predicted P (Choosing sated odor)") +
  common +
  theme(legend.position = "none")  # Remove legend if there are too many subjects

x_aOFC=ggplot(new_data,aes(x=obs_choice,y=predicted_prob)) +
  geom_point(aes(color=Cond),size=2.5) +
  geom_line(aes(group=SubID)) +
  geom_smooth(method = 'lm',color='black',lty=3) +
  stat_cor(method = 'pearson') +
  scale_color_manual(values = use.col.conds) +
  common +
  labs(x='Experimental set B choices',title = NULL,
       y='Predicted set B choices') +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.7,0.15)) +
  facet_grid(. ~ StimLoc) +
  theme(strip.background = element_rect(fill=use.col.ap.ofc[1]),
        strip.text = element_text(size=15, colour="white", face = 'bold'))


pdf(file.path(FigDir,'Day1_TMS_setB_GLMM.pdf'),4,8)
ggarrange(x_aOFC,x_pOFC,ncol = 1)
dev.off()


######## set B only, all subjects
use_trial_dat = subset(choice_dat_ss_w_base, 
                       (!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         Set == 'B' &
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.))

model_choice_0_setB <- glmer(Choice ~ StimLoc + Sess + base + (1|SubID), 
                             data = use_trial_dat,family = 'binomial')
model_choice_1_setB <- glmer(Choice ~ StimLoc * Cond + Sess + base + (1|SubID), 
                             data = use_trial_dat,family = 'binomial')
summary(model_choice_0_setB)
anova(model_choice_0_setB,model_choice_1_setB)
summary(model_choice_1_setB)


######## set A only, all subjects
use_trial_dat = subset(choice_dat_ss_w_base, 
                       (!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         Set == 'A' &
                         Cond %in% c('cTBS-sham','sham-sham'))
use_trial_dat = use_trial_dat %>%
  filter(complete.cases(.))

model_choice_0_setA <- glmer(Choice ~ StimLoc + Sess + base + (1|SubID), 
                             data = use_trial_dat,family = 'binomial')
model_choice_1_setA <- glmer(Choice ~ StimLoc * Cond + Sess + base + (1|SubID), 
                             data = use_trial_dat,family = 'binomial')
summary(model_choice_0_setA)
anova(model_choice_0_setA,model_choice_1_setA)
summary(model_choice_1_setA)

########### plot Day 1 TMS ################
######## look at cTBS-sham & sham-sham comparison ###
# look at set A and B separately
###############################

# p values come from glmer models

p_values_setA <- data.frame(
  StimLoc = c("aOFC", "pOFC"), 
  p_value = c(0.001852, 1) 
)

p_values_setB <- data.frame(
  StimLoc = c("aOFC", "pOFC"), 
  p_value = c(0.0351, p=0.0308) 
)

p_values_setA$label <- ifelse(p_values_setA$p_value < 0.05, 
                              paste("p =", format(p_values_setA$p_value, digits = 2)), 
                              "n.s.")
p_values_setB$label <- ifelse(p_values_setB$p_value < 0.05, 
                              paste("p =", format(p_values_setB$p_value, digits = 2)), 
                              "n.s.")

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

n1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill=Cond),width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = 'set A', y = "Choice of sated odor \n (post-pre meal)") + 
  common +
  coord_cartesian(ylim = c(-1,1.2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c('sham-sham', 'cTBS-sham')), 
                     size = 0, label.y = 0.8,
                     method = "t.test") +  
  geom_text(data = p_values_setA %>% filter(StimLoc == "aOFC"), 
            aes(x = 1.5, y = 1, label = label), size = 4) +
  geom_text(data = p_values_setA %>% filter(StimLoc == "pOFC"), 
            aes(x = 1.5, y = 1, label = label), size = 4)

n2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill=Cond),width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = 'set B', y = "Choice of sated odor \n (post-pre meal)") + 
  common +
  coord_cartesian(ylim = c(-1,1.2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c('sham-sham', 'cTBS-sham')), 
                     size = 0, label.y = 0.8,
                     method = "t.test") +  
  geom_text(data = p_values_setB %>% filter(StimLoc == "aOFC"), 
            aes(x = 1.5, y = 1, label = label), size = 4) +
  geom_text(data = p_values_setB %>% filter(StimLoc == "pOFC"), 
            aes(x = 1.5, y = 1, label = label), size = 4)

nn1=ggarrange(n2,n1,ncol = 1, nrow = 2)

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_by_set.pdf'),7,8)
print(nn1)
dev.off()

# focus on pOFC stim
# set A only
p1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & 
           StimLoc == 'pOFC') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = 'set A', y = "selecting sated odor \n (post-pre meal)") + 
  common +
  theme(legend.position = "none")

# set B only
p2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'pOFC') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = '', y = "selecting sated odor \n (post-pre meal)") + 
  common +
  theme(legend.position = "none")

pp1=ggpubr::ggarrange(p1,p2,ncol = 2, nrow = 1)
print(pp1)

# focus on set B
# but separated by stimulation orders
p3=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'pOFC') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  facet_wrap(~StimOrder_day1,scales = 'free_y') +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = NULL, y = "Choice of sated odor \n (post-pre meal)") + 
  common +
  theme(legend.position = "none")

# focus on set A
# but separated by stimulation orders
p4=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'pOFC') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  facet_wrap(~StimOrder_day1) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = 'set A', y = "selecting sated odor \n (post-pre meal)") + 
  common +
  theme(legend.position = "none")

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
q1 = summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  facet_grid2(rows = vars(StimLoc),
              cols = vars(StimOrder_day1),
              scales = 'free',
              strip = strip,
              axes = 'all') +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2.5, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = NULL, title = NULL, y = "Choice of sated odor \n (post-pre meal)") + 
  common +
  theme(legend.position = "none")


pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_Changes.pdf'),8,4)
print(pp1)
print(p3)
print(p4)
dev.off()

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_by_order.pdf'),10,8)
print(q1)
dev.off()

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_Changes_setB.pdf'),12,4)
ggarrange(p2,p3,widths = c(0.9,2))
dev.off()





