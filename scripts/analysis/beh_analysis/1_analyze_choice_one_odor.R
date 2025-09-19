

############################################################
# One-Odor (vs Air) Choices — Meal Effect & Sated Status
#
# INPUT
# - beh_data_processed/choice_dat.RData  (trial-level choices)
#   Loaded via project setup in scripts/utils/Setup.R
#
# WHAT THIS SCRIPT DOES
# - Filters to one-odor vs air trials (sweet=1, savory=2).
# - Labels each trial as **sated** vs **non-sated** based on whether
#   the odor type matches the devalued outcome for that subject/session.
# - Summarizes P(choose odor vs air) by Pre/Post and OdorType.
# - Draws a box+jitter plot (collapsed across sessions) comparing
#   Pre vs Post and sated vs non-sated.
# - Fits GLMMs (binomial) at the trial level to test:
#   (i) Pre vs Post change; (ii) interaction with OdorType.
#
# OUTPUTS
# - Figure P(choose odor) by OdorType with Pre/Post overlay.
# - not included in the paper

############################################################

# ChoiceType: 1,2 (1: sweet, 2:savory)
# ChosenOdor: 0,1 (1 for selecting odor)

# --- One-Odor vs Air Choices: Effect of Meal & Sated Status -------------------

rm(list = ls())

# Setup & paths
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

processed_dir <- file.path(project_folder, "beh_data_processed")

# Load full trial-level choice data
load(file = file.path(processed_dir, "choice_dat.RData"))

choice_dat_one = choice_dat %>%
  subset(ChoiceType<3) %>%
  mutate(ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('W','A'))) %>%
  mutate(OdorType=ifelse(ChoiceType==Devalued,'sated','non-sated'))

summary_choice_one = choice_dat_one %>%
  group_by(SubID,StimLoc,Cond,PrePost,OdorType) %>%
  reframe(Choice=mean(ChosenOdor,na.rm = T))

# Box + jitter: P(choose odor vs air) by OdorType, Pre vs Post (collapsed across sessions)
summary_choice_one %>%
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



# Trial-level GLMMs
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

