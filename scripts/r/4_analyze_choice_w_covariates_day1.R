

# This code performs mixed-effects modeling and group comparisons on choice data, 
# focusing only on trials involving a choice between a sated and a non-sated option.
#
# The models include two key covariates:
# (1) Learned cue weights
# (2) Selective satiation index (Didx), which is expected to positively predict choices of the sated odor.
#
# The analysis first examines the effect of TMS on Day 1
# with a similar code examines day 2


rm(list = ls())
source('Setup.R')

load(file.path(pro_dat_dir,'choice_dat_ss_w_base_values.RData'))

#############  Day 1 TMS effect ##############

# pOFC subjects
# note NODEAP_17 is a posterior cTBS subject
# with only sess 1 (sham-cTBS) and 3 (sham-sham) so remove this
# when analyzing day 1 TMS effect

use.dat.pOFC = use_choice_dat_ss %>%
  subset((!SubID=='NODEAP_17') &
           PrePost=='Post' & StimLoc=='pOFC' & 
           Cond %in% c('cTBS-sham','sham-sham')) %>%
  filter(complete.cases(.))   # remove all NA trials

model_choice_0 <- glmer(Choice ~ Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # no effect of TMS, p=0.2438
summary(model_choice_0)
# sig. valuediff, base, marginal sig. Didx

# adding Sess
model_choice_2 <- glmer(Choice ~ Didx + Sess + ValueDiff + base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond + Didx + Sess + ValueDiff + base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
anova(model_choice_2,model_choice_3) # no effect of TMS, p=0.244
summary(model_choice_2)
# sig. value diff, base
# Didx, Sess not sig
anova(model_choice_2,model_choice_0) # sess was not sig, p=0.14

use.dat.pOFC$fitted_choice <- fitted(model_choice_0, type = "response")

# aOFC subjects
use.dat.aOFC = use_choice_dat_ss %>%
  subset(PrePost=='Post' & StimLoc=='aOFC' &
           Cond %in% c('cTBS-sham','sham-sham')) %>%
  filter(complete.cases(.))   # remove all NA trials

model_choice_0 <- glmer(Choice ~ Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
model_choice_1 <- glmer(Choice ~ Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # marginal effect of TMS, p=0.066
summary(model_choice_0)
summary(model_choice_1)

model_choice_2 <- glmer(Choice ~ Sess + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
model_choice_3 <- glmer(Choice ~ Cond * Sess + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
anova(model_choice_2,model_choice_0) # sig. effect of Sess, p=0.0019
anova(model_choice_2,model_choice_3) # adding Cond*Sess improved model fit compared to model w/ only Sess
summary(model_choice_3)

use.dat.aOFC$fitted_choice <- fitted(model_choice_3, type = "response")

# calculate session-wise summary of choices
summary_choice_ss_fitted_pOFC = use.dat.pOFC %>%
  group_by(SubID,StimLoc,Cond,Sess) %>%
  reframe(avg_fitted_post_meal=mean(fitted_choice,na.rm = T)) 
summary_choice_ss_fitted_aOFC = use.dat.aOFC %>%
  group_by(SubID,StimLoc,Cond,Sess) %>%
  reframe(avg_fitted_post_meal=mean(fitted_choice,na.rm = T)) 
summary_choice_ss_fitted = rbind(summary_choice_ss_fitted_pOFC,
                                 summary_choice_ss_fitted_aOFC)
summary_choice_ss_fitted = summary_choice_ss_fitted %>%
  arrange(SubID)


p_values <- data.frame(
  StimLoc = c("aOFC", "pOFC"),  
  p_value = c('*', 'n.s.'),     
  y_start = c(1, 1),        
  y_end = c(1.03, 1.03), # text y location
  y_low = c(0.98, 0.98)
)

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',face = "bold",size = 16))

c1 = summary_choice_ss_fitted %>%
  ggplot(aes(x=Cond,y=avg_fitted_post_meal)) +
  geom_line(aes(group=SubID), position = pd_dodge, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill=Cond),width = 0.6, fatten = NULL,
               outlier.alpha = 0, alpha = 0.4) +
  stat_summary(fun = mean, geom = "crossbar", 
               fun.min = mean, fun.max = mean,  # makes top = bottom = mean
               width = 0.6,                     # controls horizontal bar width
               color = "black", 
               fatten = 0)  +   
  geom_jitter(aes(color=Cond,group=SubID), position = pd_jitter, 
              size = 2, alpha = 0.8) +
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
    aes(x = 1.5, y = y_end, label = p_value),
    inherit.aes = FALSE,
    size = 3.5) +
  annotate("text",x=1.5,y=1.1,label='Mixed-effects logistic model:') +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common +
  scale_y_continuous(
    breaks = seq(0, 1.1, by = 0.2),  # tick marks at 0, 0.2, 0.4, ..., 1.0
    limits = c(0, 1.1)               # ensure the axis runs from 0 to 1
  ) +
  labs(x = NULL, title = NULL, y = "Predicted probability \n of choosing sated odor") + 
  theme(legend.position = "none")

pdf(file.path(FigPaperDir,'Day1_TMS_ChoiceSatedOdor_fitted.pdf'),7,4)
print(c1)
dev.off()



library(dplyr)
library(ggplot2)

# Filter only the two conditions of interest
plot_data <- use.dat.aOFC %>%
  filter(Cond %in% c("sham-sham", "cTBS-sham")) %>%
  group_by(SubID, Sess, StimOrder_day1, Cond) %>%
  #reframe(mean = mean(Choice,na.rm=T)) %>%
  reframe(mean = mean(fitted_choice,na.rm=T))

# Plot
ggplot(plot_data, aes(x = Sess, y = mean, fill = Cond)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(group = as.factor(Sess))) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # optional: show individual points
  facet_wrap(~ StimOrder_day1) +
  scale_fill_manual(values = use.col.conds) +
  labs(
    title = NULL,
    x = "Session number",
    y = "Choosing sated odors"
  ) +
  common




get_subject_summary <- function(data, stimloc_label) {
  data %>%
    dplyr::group_by(SubID, Cond) %>%
    summarise(
      mean_actual = mean(Choice),
      mean_predicted = mean(fitted_choice),
      .groups = "drop"  # optional: remove grouping
    ) %>%
    mutate(StimLoc = stimloc_label)
}

# Apply to both groups
subject_summary_pOFC <- get_subject_summary(use.dat.pOFC, "pOFC")
subject_summary_aOFC <- get_subject_summary(use.dat.aOFC, "aOFC")

# Combine
combined_summary <- bind_rows(subject_summary_pOFC, subject_summary_aOFC)

# Compute correlations
cor_pOFC <- cor(subject_summary_pOFC$mean_predicted, subject_summary_pOFC$mean_actual)
cor_aOFC <- cor(subject_summary_aOFC$mean_predicted, subject_summary_aOFC$mean_actual)

# Text
cor_text_pOFC <- paste0("r = ", round(cor_pOFC, 2))
cor_text_aOFC <- paste0("r = ", round(cor_aOFC, 2))

label_df <- data.frame(
  StimLoc = c("pOFC", "aOFC"),
  label = c(cor_text_pOFC, cor_text_aOFC),
  x = 0.2,  # adjust as needed
  y = 0.9   # adjust as needed
)

bt_sub_corr = ggplot(combined_summary, aes(x = mean_predicted, y = mean_actual, color = Cond)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = use.col.conds) + 
  labs(
    x = "Mean fitted choice probability",
    y = "Mean actual choice rate",
    title = NULL
  ) +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),  # tick marks at 0, 0.2, 0.4, ..., 1.0
    limits = c(0, 1)               # ensure the axis runs from 0 to 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),  # tick marks at 0, 0.2, 0.4, ..., 1.0
    limits = c(0, 1)               # ensure the axis runs from 0 to 1
  ) +
  common +
  geom_text(data = label_df, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, size = 5) +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir,'Day1_Choice_fitted_actual.pdf'),8,4)
print(bt_sub_corr)
dev.off()


#########################################

# Step 1: Filter out actual values for both conditions
actual_sham <- combined_summary %>%
  filter(Cond == "sham-sham") %>%
  select(SubID, StimLoc, mean_actual) %>%
  rename(mean_actual_sham = mean_actual)

actual_cTBS <- combined_summary %>%
  filter(Cond == "cTBS-sham") %>%
  select(SubID, StimLoc, mean_actual) %>%
  rename(mean_actual_cTBS = mean_actual)

# Step 2: Merge and calculate actual TMS effect
actual_effect_df <- left_join(actual_sham, actual_cTBS, by = c("SubID", "StimLoc")) %>%
  mutate(TMS_effect_actual = mean_actual_cTBS - mean_actual_sham)

# Step 3: Repeat for predicted values
predicted_sham <- combined_summary %>%
  filter(Cond == "sham-sham") %>%
  select(SubID, StimLoc, mean_predicted) %>%
  rename(mean_predicted_sham = mean_predicted)

predicted_cTBS <- combined_summary %>%
  filter(Cond == "cTBS-sham") %>%
  select(SubID, StimLoc, mean_predicted) %>%
  rename(mean_predicted_cTBS = mean_predicted)

predicted_effect_df <- left_join(predicted_sham, predicted_cTBS, by = c("SubID", "StimLoc")) %>%
  mutate(TMS_effect_predicted = mean_predicted_cTBS - mean_predicted_sham)

# Step 4: Merge actual and predicted effect
tms_effect_summary <- left_join(actual_effect_df %>% select(SubID, StimLoc, TMS_effect_actual),
                                predicted_effect_df %>% select(SubID, StimLoc, TMS_effect_predicted),
                                by = c("SubID", "StimLoc"))

# Step 5: Save to CSV
write_csv(tms_effect_summary, file.path(pro_dat_dir,"TMS_effect_summary_day1.csv"))



#########################################

library(pROC)

# Function to compute ROC dataframe for one group
compute_roc_df <- function(data, stimloc_label) {
  conds <- unique(data$Cond)
  
  roc_list <- lapply(conds, function(cond_name) {
    dat <- filter(data, Cond == cond_name)
    roc_obj <- roc(dat$Choice, dat$fitted_choice)
    auc_val <- round(auc(roc_obj), 2)
    data.frame(
      TPR = roc_obj$sensitivities,
      FPR = 1 - roc_obj$specificities,
      Cond = cond_name,
      AUC = auc_val,
      StimLoc = stimloc_label
    )
  })
  
  bind_rows(roc_list)
}

# Compute for pOFC and aOFC
roc_df_pOFC <- compute_roc_df(use.dat.pOFC, "pOFC")
roc_df_aOFC <- compute_roc_df(use.dat.aOFC, "aOFC")

# Combine
roc_df_all <- bind_rows(roc_df_pOFC, roc_df_aOFC)

auc_summary <- roc_df_all %>%
  group_by(StimLoc, Cond) %>%
  summarise(AUC = first(AUC))  # AUC is constant per curve

auc_labels <- auc_summary %>%
  mutate(
    label = paste0("AUC = ", round(AUC, 2)),
    x = case_when(
      Cond == "sham-sham" ~ 0.6,
      Cond == "cTBS-sham" ~ 0.6
    ),
    y = case_when(
      Cond == "sham-sham" ~ 0.2,
      Cond == "cTBS-sham" ~ 0.1
    )
  )


# Plot
roc = ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Cond)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_color_manual(values = use.col.conds) + 
  geom_text(data = auc_labels, aes(x = x, y = y, label = label, color = Cond), inherit.aes = FALSE) +
  labs(
    x = "False positive rate (FPR)",
    y = "True positive rate (TPR)",
    title = NULL
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),  # tick marks at 0, 0.2, 0.4, ..., 1.0
    limits = c(0, 1)               # ensure the axis runs from 0 to 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),  # tick marks at 0, 0.2, 0.4, ..., 1.0
    limits = c(0, 1)               # ensure the axis runs from 0 to 1
  ) +
  common + theme(legend.position = 'none')

pdf(file.path(FigPaperDir,'Day1_Choice_roc.pdf'),8,4)
print(roc)
dev.off()


##################################

# put aOFC and pOFC together) 
# to see if there's interaction b/t stim loc & TMS cond
use.dat = subset(use_choice_dat_ss, PrePost=='Post' &
                   Cond %in% c('cTBS-sham','sham-sham'))
model_choice_2 <- glmer(Choice ~ StimLoc * Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_1 <- glmer(Choice ~ StimLoc + Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')
model_choice_0 <- glmer(Choice ~ Cond + Didx + ValueDiff + base + (1|SubID), 
                        data = use.dat,family = 'binomial')

summary(model_choice_0)
summary(model_choice_1)
summary(model_choice_2)
anova(model_choice_0,model_choice_1)
anova(model_choice_2,model_choice_1)
anova(model_choice_2,model_choice_0)

