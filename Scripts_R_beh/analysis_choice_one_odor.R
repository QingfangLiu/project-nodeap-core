
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

summary_choice_one = choice_dat_one %>%
  group_by(SubID,StimLoc,Cond,PrePost,OdorType) %>%
  reframe(Choice=mean(ChosenOdor,na.rm = T))

# collapsing across sessions
c0 = summary_choice_one %>%
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
  annotate('text',x=1.5,y=1.15,label='***') +
  annotate('segment',x=1,y=1.1,xend=2,yend=1.1) +
  annotate('segment',x=1,y=1.05,xend=1,yend=1.1) +
  annotate('segment',x=2,y=1.05,xend=2,yend=1.1) +
  labs(x = NULL, title = NULL, y = "Choice of odor (vs. air)") + common +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir,'Day2_ChoiceOneOdor.pdf'),4,4)
print(c0)
dev.off()


# use trial-level data 
# test: if less odor choices with sated odor
model_choice_0 <- glmer(ChosenOdor ~ (1|SubID), 
                        data = choice_dat_one,family = 'binomial')
model_choice_1 <- glmer(ChosenOdor ~ PrePost + (1|SubID), 
                        data = choice_dat_one,family = 'binomial')
model_choice_2 <- glmer(ChosenOdor ~ PrePost * OdorType + (1|SubID), 
                        data = choice_dat_one,family = 'binomial')

# choice change from pre to post meal?
anova(model_choice_0,model_choice_1)
summary(model_choice_1)

# does odor type (sated vs non-sated) affect choice change?
anova(model_choice_1,model_choice_2)
summary(model_choice_1)
summary(model_choice_2)

