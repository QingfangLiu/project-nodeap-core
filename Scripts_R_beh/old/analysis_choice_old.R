
rm(list = ls())
source('Setup.R')

###############################
# analyze choice data in the pre-meal
# and post-meal
###############################

choice_dat = read.xlsx('../ProcessedData/Choices.xlsx') %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre','Post')),
         Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre','Post')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Cond=factor(Cond),
         Set=factor(Set),
         StimLoc=factor(StimLoc))

# only look at choice with one sweet and one savory (denoted as '_ss')
# perhaps where the effect is the strongest
choice_dat_ss = choice_dat %>%
  subset(ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) # Choice==1 means selecting sated odor

# summary of the choice_ss data
summary_choice_ss = choice_dat_ss %>%
  group_by(SubID,StimOrder,StimLoc,Cond,Sess,PrePost,Set) %>%
  dplyr::summarise(Choice=mean(Choice,na.rm = T))

# separted by stimulation location
p1=summary_choice_ss %>%
  ggplot(aes(x=Cond,y=Choice,fill=interaction(PrePost,Set))) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(aes(color=interaction(PrePost,Set)),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
              size = 1, alpha = 0.8) +
  facet_wrap(~StimLoc) +
  labs(x = "", y = "Choice of Sated Odor (Sweet vs. Savory)")

# collapsing aOFC and pOFC stimulations
p2=summary_choice_ss %>%
  ggplot(aes(x=Cond,y=Choice,fill=interaction(PrePost,Set))) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(aes(color=interaction(PrePost,Set)),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
              size = 1, alpha = 0.8) +
  labs(x = "", y = "Choice of Sated Odor (Sweet vs. Savory)")

pdf(file.path(FigDir,'ChoiceSatedOdor.pdf'),8,5)
print(p1)
print(p2)
dev.off()

# focus on Set B
# only on pOFC
use.dat = subset(choice_dat_ss, Set=='B' & StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + (1|SubID), data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ (1|SubID), data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)
# lower choice in sham-sham condition

# only on aOFC
use.dat = subset(choice_dat_ss, Set=='B' & StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + (1|SubID), data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ (1|SubID), data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# no sig difference b/t conditions

# if including both aOFC & pOFC
use.dat = subset(choice_dat_ss, Set=='B' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ Cond + StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond * StimLoc + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_2,model_choice_3)
# interaction efffect between condition and location
# was significant

# Set A
# if only focused on pOFC
use.dat = subset(choice_dat_ss, Set=='A' &
                   StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond * PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# p = 0.090
# with trend but not sig

# if only focused on aOFC
use.dat = subset(choice_dat_ss, Set=='A' &
                   StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond * PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + PrePost + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
# no sig

###############################################
# improved from above by considering session-wise odor bias
# using pre(A) as baseline
# then get difference from pre to post
###############################################

pre_choice = subset(summary_choice_ss,Set=='A' & PrePost=='Pre')$Choice
post_choice_A = subset(summary_choice_ss,Set=='A' & PrePost=='Post')$Choice
post_choice_B = subset(summary_choice_ss,Set=='B' & PrePost=='Post')$Choice
conds = subset(summary_choice_ss,Set=='A' & PrePost=='Pre')$Cond
# check if there is pre-post correlations of choices
df_pre_post_summary = data.frame(
  'conds'=conds,
  'pre_choice'=pre_choice,
  'post_choice_A'=post_choice_A,
  'post_choice_B'=post_choice_B
)

pcorr1=ggplot(df_pre_post_summary,aes(x=pre_choice,y=post_choice_A,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")
pcorr2=ggplot(df_pre_post_summary,aes(x=pre_choice,y=post_choice_B,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")
pcorr3=ggplot(df_pre_post_summary,aes(x=post_choice_A,y=post_choice_B,color=conds)) +
  geom_point() +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Choice of Sated Odor")

avg_changeAB = (post_choice_A + post_choice_B)/2 - pre_choice
avg_changeB = post_choice_B - pre_choice
avg_changeA = post_choice_A - pre_choice
avg_BoverA = post_choice_B - post_choice_A

# create a corrected summary dataset
# containing choice update for each subject and each TMS condition
summary_choice_corrected = summary_choice_ss %>%
  ungroup() %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  select(SubID,StimLoc,StimOrder,Cond,Sess) %>%
  mutate(ChoiceChangeAB=avg_changeAB, # both A and B
         ChoiceChangeA=avg_changeA, # only A
         ChoiceChangeB=avg_changeB, # only B
         ChoiceBoverA=avg_BoverA)

# save this corrected summary for later use
save(summary_choice_corrected,file = '../ProcessedData/Choice_corrected_dat.RData')

# set A and B 
cc1=summary_choice_corrected %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre',
       y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")

# set B only
cc2=summary_choice_corrected %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(B only) - Pre',
       y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")

# set A only
cc3=summary_choice_corrected %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'Post(A only) - Pre',
       y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")

# both set A and B
# separated by two stim locations
levels(summary_choice_corrected$StimLoc)=c('aOFC','pOFC')

c1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = '',
       y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")

# based on c1, focus on pOFC stim
# separated by diff stim orders
c2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  subset(StimLoc %in% c('Posterior')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc + StimOrder) +
  scale_fill_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  scale_color_manual(values = c("sham-sham" = "gray","cTBS-sham" = "green","sham-cTBS" = "blue")) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre', y = "Choice of Sated Odor") + common +
  theme(legend.position = "none")

# or focus on these two conditions, but look at how satiation effect
# changes over sessions
# did not see higher satiation in the first session number
c3=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  subset(StimLoc %in% c('Posterior')) %>%
  mutate(Sess=factor(Sess)) %>%
  ggplot(aes(x=Sess,y=ChoiceChangeAB,fill=Sess)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  scale_fill_manual(values = use.col.sess) +
  scale_color_manual(values = use.col.sess) +
  geom_jitter(aes(color=Sess,group=SubID), position = pd, size = 2, alpha = 0.8) +
  labs(x = "", title = 'Post(A & B) - Pre',
       y = "Choice of Sated Odor") + common + theme(legend.position = "none")

pdf(file.path(FigDir,'ChoiceSatedOdor_Changes.pdf'),8,7)
print(pcorr1)
print(pcorr2)
print(pcorr3)
print(cc1)
print(cc2)
print(cc3)
dev.off()

pdf(file.path(FigDir,'Day2_TMS_ChoiceSatedOdor_Changes.pdf'),7,5)
print(c1)
#print(c2)
#print(c3)
dev.off()


#######################################################
# stat testing on Day 2 TMS effect
# test using Post(A & B) - Pre, summarized choice data
# pOFC
# overall devaluation effect regardless of conditions: yes
# (later realizing this also contains cTBS-sham condition)
value_Y1Z1 = subset(summary_choice_corrected,StimLoc=='Posterior')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.013

value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                    StimLoc=='Posterior')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='Posterior')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? No
wilcox.test(value_Z1,alternative = 'l') # lower than zero? Yes, p=0.026
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s

# stat testing on aOFC
# overall devaluation effect regardless of conditions: yes
value_Y1Z1 = subset(summary_choice_corrected,
                    StimLoc=='Anterior')$ChoiceChangeAB
wilcox.test(value_Y1Z1,alternative = 'l') # lower than zero? Yes, p=0.010

value_Y1 = subset(summary_choice_corrected, Cond=='sham-cTBS'& 
                    StimLoc=='Anterior')$ChoiceChangeAB
value_Z1 = subset(summary_choice_corrected, Cond=='sham-sham'& 
                    StimLoc=='Anterior')$ChoiceChangeAB
wilcox.test(value_Y1,alternative = 'l') # lower than zero? p = 0.065
wilcox.test(value_Z1,alternative = 'l') # lower than zero? n.s.
wilcox.test(value_Y1,value_Z1,paired = T,alternative = 'g') # n.s


#######################################################
# testing Day 2 TMS effect at trial level
#######################################################

# go back to choice_dat_ss to do the testing at trial level
# first calculate a session-wise pre-meal bias for each sub
pre_choice_df = summary_choice_ss %>%
  subset(Set=='A' & PrePost=='Pre') %>%
  group_by(SubID,Cond) %>%
  dplyr::summarise(base=Choice)
# insert this base value to the whole choice dat ss df
choice_dat_ss_w_base = choice_dat_ss %>%
  mutate(base=NA)
for(i in 1:nrow(choice_dat_ss_w_base)){
  choice_dat_ss_w_base$base[i] = pre_choice_df$base[
    pre_choice_df$SubID==choice_dat_ss_w_base$SubID[i] &
      pre_choice_df$Cond==choice_dat_ss_w_base$Cond[i]]
}

# with session-wise baseline odor preference as a covariate
# focus on the two conditions of interest
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' & 
                   StimLoc=='Posterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)
# adding condition sig improved model fit
# and both 'Cond' & 'base' variables are sig

# repeate this analysis on aOFC
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Anterior' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0)
summary(model_choice_1)
# adding condition did not improve model fit

# put aOFC and pOFC together 
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ StimLoc * Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ StimLoc + Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)
summary(model_choice_1)
summary(model_choice_2)
# adding StimLoc itself did not improve model fit, but 
# adding the interaction b/t StimLoc and Cond did
# this nicely aligns with the hypothesis that pOFC 
# selectively disrupted outcome devaluation

###############################
## Day1 TMS effect on probe choice
######## look at cTBS-sham & sham-sham comparison ###
# focus on posterior OFC stimulation
###############################

# set A only
p1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = 'set A', y = "selecting sated odor (post-pre meal)") + 
  common +
  theme(legend.position = "none")

# set B only
p2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = '', y = "Choice of sated odor") + 
  common +
  theme(legend.position = "none")

pp1=ggpubr::ggarrange(p1,p2,ncol = 2, nrow = 1)

# focus on set B
# but separated by stimulation orders
p3=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeB,fill=Cond)) +
  facet_wrap(~StimOrder) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = 'set B', y = "selecting sated odor") + 
  common +
  theme(legend.position = "none")

# focus on set A
# but separated by stimulation orders
p4=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceChangeA,fill=Cond)) +
  facet_wrap(~StimOrder) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  labs(x = "", title = 'set A', y = "selecting sated odor") + 
  common +
  theme(legend.position = "none")

# effect of session time 1,2,3 on set A
p5=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  mutate(Sess=factor(Sess)) %>%
  ggplot(aes(x=Sess,y=ChoiceChangeA,fill=Sess)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(x = "", title = 'set A', y = "selecting sated odor") + 
  common + theme(legend.position = "none")

# effect of session time 1,2,3 on set B
p6=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  mutate(Sess=factor(Sess)) %>%
  ggplot(aes(x=Sess,y=ChoiceChangeB,fill=Sess)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(x = "", title = 'set B', y = "selecting sated odor") + 
  common + theme(legend.position = "none")

# plot B over A to assess day 1 TMS effect on devaluation
# where A serves as the baseline
# (disruption should be choosing more of the sated odor)
# but did not really see the effect
p7=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','cTBS-sham') & StimLoc == 'Posterior') %>%
  ggplot(aes(x=Cond,y=ChoiceBoverA,fill=Cond)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  scale_fill_manual(values = use.col.conds) +
  scale_color_manual(values = use.col.conds) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  labs(x = "", title = 'B - A', y = "Choice of Sated Odor") + 
  common + theme(legend.position = "none")

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_Changes.pdf'),8,7)
print(pp1)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_Changes_setB.pdf'),4,4)
print(p2)
dev.off()

#######################################################
# stat testing on Day 1 TMS effect
# test using summarized choice data
# focus on pOFC stim
#######################################################

# on set A
value_XZ_A = subset(summary_choice_corrected,
                    Cond %in% c('sham-sham','cTBS-sham') &
                      StimLoc=='Posterior')$ChoiceChangeA
wilcox.test(value_XZ_A,alternative = 'l') # lower than zero? Yes, p=0.00037
value_X_A = subset(summary_choice_corrected, Cond=='cTBS-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeA
value_Z_A = subset(summary_choice_corrected, Cond=='sham-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeA
wilcox.test(value_X_A,alternative = 'l') # lower than zero? Yes, p=0.0088
wilcox.test(value_Z_A,alternative = 'l') # lower than zero? Yes, p=0.011
wilcox.test(value_X_A,value_Z_A,paired = T,alternative = 'g') # n.s

# on set B
value_XZ_B = subset(summary_choice_corrected,
                    Cond %in% c('sham-sham','cTBS-sham') &
                      StimLoc=='Posterior')$ChoiceChangeB
wilcox.test(value_XZ_B,alternative = 'l') # lower than zero? p = 0.064
value_X_B = subset(summary_choice_corrected, Cond=='cTBS-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeB
value_Z_B = subset(summary_choice_corrected, Cond=='sham-sham'& 
                     StimLoc=='Posterior')$ChoiceChangeB
wilcox.test(value_X_B,alternative = 'l') # lower than zero? No, p=0.29
wilcox.test(value_Z_B,alternative = 'l') # lower than zero? No, p=0.092
wilcox.test(value_X_B,value_Z_B,paired = T,alternative = 'g') # n.s

# on set B over set A
value_X_BoverA = subset(summary_choice_corrected, Cond=='cTBS-sham'& 
                     StimLoc=='Posterior')$ChoiceBoverA
value_Z_BoverA = subset(summary_choice_corrected, Cond=='sham-sham'& 
                     StimLoc=='Posterior')$ChoiceBoverA
wilcox.test(value_X_BoverA,alternative = 'g') # greater than zero? p=0.054
wilcox.test(value_Z_BoverA,alternative = 't') # p=0.4
wilcox.test(value_X_BoverA,value_Z_BoverA,paired = T,alternative = 'g') # p=0.165


#######################################################
# testing Day 1 TMS effect at trial level
#######################################################

# with considering session-wise baseline odor preference
# pOFC stimulation only
# on set A
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                     StimLoc=='Posterior' & Set=='A' &
                     Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of condition
summary(model_choice_0)
coef_summary <- summary(model_choice_1)
coef_est <- coef_summary$coefficients["CondcTBS-sham", "Estimate"]
coef_se <- coef_summary$coefficients["CondcTBS-sham", "Std. Error"]
z_value <- coef_est / coef_se
p_value <- 1 - pnorm(z_value)
cat("Z-value:", z_value, "\n")
cat("One-sided p-value:", p_value, "\n") # one-sided p: 0.84

# on Set B
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' & Set=='B' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # p=0.075
summary(model_choice_1)
coef_summary <- summary(model_choice_1)
coef_est <- coef_summary$coefficients["CondcTBS-sham", "Estimate"]
coef_se <- coef_summary$coefficients["CondcTBS-sham", "Std. Error"]
z_value <- coef_est / coef_se
p_value <- 1 - pnorm(z_value)
cat("Z-value:", z_value, "\n")
cat("One-sided p-value:", p_value, "\n") # one-sided p: 0.037

# interactions b/t set A and B

use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   StimLoc=='Posterior' & 
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_2 <- glmer(Choice ~ Set * Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Set + Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
summary(model_choice_1) # set B has higher sated choices
anova(model_choice_2,model_choice_1)


#######################################################
# change of choice in the trial course on pre-meal task
# I realized that disruption of identity learning should
# also be shown in the pre-meal choice task
# (1) if we assume people have a session-wise pre-meal odor
# baseline, then people need to learn which stimuli to choose
# (2) if we believe such a learning process, then we
# can't use the whole trial course to serve as 
# pre-meal odor preference

#######################################################

# only include SW/SA choices, 18 trials per run
pre_meal_choice = choice_dat %>%
  subset(PrePost=='Pre' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) # Choice==1 means selecting sated odor
  
# split choices into 3 parts and take each mean
pre_meal_choice = pre_meal_choice %>%
  mutate(TrialCount=rep(1:18,nrow(pre_meal_choice)/18)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 3)) %>%
  group_by(SubID, StimLoc, Cond, part) %>%
  summarise(avg_choice = mean(Choice,na.rm=T)) %>%
  ungroup()

pre1=pre_meal_choice %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Pre-meal task',
       title = 'aOFC stim') +
  scale_color_manual(values = use.col.conds) 

pre2=pre_meal_choice %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Pre-meal task',
       title = 'pOFC stim') + 
  scale_color_manual(values = use.col.conds) 


pdf(file.path(FigDir,'Pre_meal_choices.pdf'),8,7)
print(pre1)
print(pre2)
dev.off()


#######################################################
# post-meal choices
#######################################################

# only include SW/SA choices, 36 trials per run
post_meal_choice = choice_dat %>%
  subset(PrePost=='Post' & ChoiceType==3) %>%
  mutate(ChosenOdor=mapvalues(ChosenOdor,from=c(1,2),to=c('W','A'))) %>%
  mutate(Choice=ifelse(ChosenOdor==Devalued,1,0)) # Choice==1 means selecting sated odor
  
# split choices into 3 parts and take each mean
post_meal_choice = post_meal_choice %>%
  mutate(TrialCount=rep(1:36,nrow(post_meal_choice)/36)) %>%
  group_by(SubID, Cond) %>%
  mutate(part = ntile(TrialCount, 6)) %>%
  group_by(SubID, StimLoc, Cond, part) %>%
  summarise(avg_choice = mean(Choice,na.rm=T)) %>%
  ungroup()

post1=post_meal_choice %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Post-meal task',
       title = 'aOFC stim') +
  scale_color_manual(values = use.col.conds) 

post2=post_meal_choice %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x=part,y=avg_choice,color=Cond)) +
  facet_wrap(~SubID) +
  geom_jitter(height = 0, width = 0.1) +
  geom_smooth(fill = NA) +
  labs(y='Selecting sated odor',x='Post-meal task',
       title = 'pOFC stim') +
  scale_color_manual(values = use.col.conds) 
  
pdf(file.path(FigDir,'Post_meal_choices.pdf'),8,7)
print(post1)
print(post2)
dev.off()
