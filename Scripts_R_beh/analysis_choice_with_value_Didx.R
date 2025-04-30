
# this code does mixed-effect models and comparisons on choice data
# only look at choices between sated option and non-sated option
# with added covariates of 
# (1) learned cue weights
# (2) selective satiation index (Didx): should positively scale with sated odor choice

# looked at Day 2 TMS effect first
# then Day 1 TMS effect

# previously looked data for two sets separately but erased for cleanness
# last update: 4/30/2025


rm(list = ls())
source('Setup.R')

# load the choice (ss) with baseline and values
load(file = '../ProcessedData/choice_dat_ss_w_base_values.RData')
load(file = '../ProcessedData/SelectSate_dat.RData')

# create a var ValueDiff from difference of sated and non-sated
choice_dat_ss_w_base = use_choice_dat_ss %>%
  mutate(ValueDiff = ValueSated - ValueNonSated,
         Didx = NA)

# add the session-wise measure Didx to the trial-wise dataset
for(i in 1:nrow(choice_dat_ss_w_base)){
  choice_dat_ss_w_base$Didx[i] = SelectSate_dat$Didx[
    SelectSate_dat$SubID==choice_dat_ss_w_base$SubID[i] &
      SelectSate_dat$Cond==choice_dat_ss_w_base$Cond[i]]
}

################## Day 2 TMS effect ##################

# focus on pOFC
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' & 
                   StimLoc=='pOFC' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # sig TMS p = 0.00034
summary(model_choice_1)
# everything sig: TMS, Valuediff, Didx, base

# focus on aOFC
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' & 
                   StimLoc=='aOFC' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_1,model_choice_0) # not sig TMS
summary(model_choice_0)
# everything else sig: Valuediff, Didx, base


# put aOFC and pOFC together 
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   Cond %in% c('sham-cTBS','sham-sham'))
model_choice_2 <- glmer(Choice ~ StimLoc * Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ StimLoc + Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)
summary(model_choice_2)
# sig: TMS&Loc interaction, Didx, ValueDiff, base


###### then how about Day 1 TMS effect? ##############
# note NODEAP_17 is a posterior cTBS subject
# with only sess 1 (sham-cTBS) and 3 (sham-sham) so remove this
# when analyzing day 1 TMS effect

# pOFC subjects
use_trial_dat_pOFC = choice_dat_ss_w_base %>%
                  subset((!SubID=='NODEAP_17') &
                         PrePost=='Post' &
                         StimLoc=='pOFC' & 
                         Cond %in% c('cTBS-sham','sham-sham')) %>%
  mutate(CuePair=factor(CuePair)) %>%
  filter(complete.cases(.))   # remove all NA trials

model_choice_0 <- glmer(Choice ~ Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of TMS, p=0.2438
summary(model_choice_0)
# sig. valuediff, base, marginal sig. Didx

# adding Sess
model_choice_2 <- glmer(Choice ~ Didx + Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond + Didx + Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
anova(model_choice_2,model_choice_3) # no effect of TMS, p=0.244
summary(model_choice_2)
# sig. value diff, base
# Didx, Sess not sig
anova(model_choice_2,model_choice_0) # sess was not sig, p=0.14

# removing Didx
model_choice_4 <- glmer(Choice ~ ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
model_choice_5 <- glmer(Choice ~ Cond + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_pOFC,family = 'binomial')
anova(model_choice_4,model_choice_5) # no effect of TMS, p=0.129
summary(model_choice_4)

predictions <- ggpredict(model_choice_5, terms = c("Cond"))
# Model-Adjusted Data
use_trial_dat_pOFC$fitted_choice <- fitted(model_choice_5, type = "response")

p_pOFC = ggplot(predictions, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = use_trial_dat_pOFC, 
              aes(x = Cond, y = fitted_choice, color = Cond), 
              width = 0.3, alpha = 0.1, show.legend = F) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = use.col.conds) +
  labs(y = "Choice of sated odor (adjusted)", x = NULL, title = 'pOFC') +
  common


# aOFC subjects
use_trial_dat_aOFC = choice_dat_ss_w_base %>%
  subset(PrePost=='Post' &
           StimLoc=='aOFC' & 
           Cond %in% c('cTBS-sham','sham-sham')) %>%
  mutate(CuePair=factor(CuePair)) %>%
  filter(complete.cases(.))   # remove all NA trials

model_choice_0 <- glmer(Choice ~ Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # marginal effect of TMS, p=0.066
summary(model_choice_0)
summary(model_choice_1)

model_choice_2 <- glmer(Choice ~ Sess + Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond + Didx + Sess + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC,family = 'binomial')
model_choice_4 <- glmer(Choice ~ Cond * Sess + Didx + ValueDiff + base + (1|SubID), 
                        data = use_trial_dat_aOFC,family = 'binomial')
anova(model_choice_2,model_choice_0) # sig. effect of Sess, p=0.0019
anova(model_choice_2,model_choice_3) # no effect of TMS after having sess, p=0.23
anova(model_choice_2,model_choice_4)
anova(model_choice_4,model_choice_3) # interaction of Sess & TMS
anova(model_choice_4,model_choice_0)
anova(model_choice_4,model_choice_0,model_choice_2,model_choice_3)
summary(model_choice_4)

predictions <- ggpredict(model_choice_3, terms = c("Cond"))
# Model-Adjusted Data
use_trial_dat_aOFC$fitted_choice <- fitted(model_choice_3, type = "response")

p_aOFC = ggplot(predictions, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = use_trial_dat_aOFC, 
              aes(x = Cond, y = fitted_choice, colour = Cond), 
              width = 0.3, alpha = 0.1, show.legend = F) +
  scale_color_manual(values = use.col.conds) +
  coord_cartesian(ylim = c(0,1)) +
  labs(y = "Choice of sated odor (adjusted)", x = NULL, title = 'aOFC') +
  common 

pdf(file.path(FigDir,'Day1_TMS_ChoiceSatedOdor_covariates_explained.pdf'),7,4)
ggarrange(p_aOFC,p_pOFC)
dev.off()



# put aOFC and pOFC together 
# to see if there's interaction b/t stim loc & TMS cond
use.dat = subset(choice_dat_ss_w_base, PrePost=='Post' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_2 <- glmer(Choice ~ StimLoc * Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
summary(model_choice_2)
model_choice_1 <- glmer(Choice ~ StimLoc + Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
summary(model_choice_1)
model_choice_0 <- glmer(Choice ~ Cond + Didx +
                          ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
summary(model_choice_0)

anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)


# visualize this based on the model with both aOFC and pOFC
new_data <- expand.grid(
  StimLoc = c("aOFC", "pOFC"),
  Cond = c("sham-sham", "cTBS-sham"),
  Didx = mean(use.dat$Didx, na.rm = TRUE),
  ValueDiff = mean(use.dat$ValueDiff, na.rm = TRUE),
  base = mean(use.dat$base, na.rm = TRUE)
)
new_data$predicted <- predict(model_choice_2, 
                              newdata = new_data, 
                              type = "response", re.form = NA)

q1 = ggplot(new_data, aes(x = StimLoc, y = predicted, fill = Cond)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Predicted Probability of Choice", x = NULL, title = NULL) +
  common +
  theme(legend.title = element_blank())
q2 = ggplot(new_data, aes(x = Cond, 
                     y = predicted, group = StimLoc, color = StimLoc)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = use.col.ap.ofc) +
  labs(y = "Choice of sated odor", x = NULL, title = NULL) +
  common +
  theme(legend.title = element_blank())

predictions <- ggpredict(model_choice_2, terms = c("Cond", "StimLoc"))
ggplot(predictions, aes(x = x, y = predicted, group = group, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Predicted Probability of Choice", x = NULL,
       title = NULL,
       color = "StimLoc") +
  scale_color_manual(values = use.col.ap.ofc) +
  common

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

ggplot(predictions, aes(x = x, y = predicted, color = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.3, linewidth = 1) +
  facet_grid2(cols = vars(group), 
              scales = 'free',
              strip = strip,
              axes = 'all') +
  labs(x = NULL, y = "Predicted Probability of Choice",
    title = NULL) +
  scale_color_manual(values = use.col.conds) +
  common +
  theme(legend.position = "none") 

