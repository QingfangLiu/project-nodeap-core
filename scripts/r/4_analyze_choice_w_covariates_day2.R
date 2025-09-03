
# This code performs mixed-effects modeling and group comparisons on choice data, 
# focusing only on trials involving a choice between a sated and a non-sated option.
#
# The models include two key covariates:
# (1) Learned cue weights
# (2) Selective satiation index (Didx), which is expected to positively predict choices of the sated odor.
#
# The analysis first examines the effect of TMS on Day 2
# with a similar code examines day 1

rm(list = ls())
source('Setup.R')

load(file.path(pro_dat_dir,'choice_dat_ss_w_base_values.RData'))

# pOFC group
use.dat.pOFC = use_choice_dat_ss %>%
                  subset(PrePost=='Post' & 
                         StimLoc=='pOFC' & 
                         Cond %in% c('sham-cTBS','sham-sham')) %>%
                  filter(complete.cases(.)) 

model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + Didx +
                          base + (1|SubID), 
                        data = use.dat.pOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # sig TMS p = 0.00034
summary(model_choice_1)
# everything sig: TMS, Valuediff, Didx, base

use.dat.pOFC$fitted_choice <- fitted(model_choice_1, type = "response")

# add in response to a reviewer

# get unique subject IDs
subs <- unique(use.dat.pOFC$SubID)

# storage
results <- data.frame(SubID = subs, p_value = NA)

for (i in seq_along(subs)) {
  # filter out one subject
  dat_sub <- filter(use.dat.pOFC, SubID != subs[i])
  
  # fit the two models
  m1 <- glmer(Choice ~ Cond + ValueDiff + Didx + base + (1 | SubID),
              data = dat_sub, family = "binomial")
  
  m0 <- glmer(Choice ~ ValueDiff + Didx + base + (1 | SubID),
              data = dat_sub, family = "binomial")
  
  # likelihood ratio test
  lrt <- anova(m1, m0)
  
  # extract p-value (2nd row, "Pr(>Chisq)")
  results$p_value[i] <- lrt$`Pr(>Chisq)`[2]
}

print(results)


results_ordered <- results %>% arrange(p_value)
print(results_ordered)

# Summarize LOSO results
summary_stats <- results %>%
  summarise(
    min_p   = min(p_value, na.rm = TRUE),
    max_p   = max(p_value, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    mean_p   = mean(p_value, na.rm = TRUE)
  )

print(summary_stats)

# Identify subjects with p > 0.05
influential <- results %>%
  filter(p_value > 0.05)

print(influential)





# aOFC group
use.dat.aOFC = use_choice_dat_ss %>%
                subset(PrePost=='Post' & StimLoc=='aOFC' & 
                       Cond %in% c('sham-cTBS','sham-sham')) %>%
                filter(complete.cases(.)) 
model_choice_1 <- glmer(Choice ~ Cond + ValueDiff + Didx + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
model_choice_0 <- glmer(Choice ~ ValueDiff + Didx + base + (1|SubID), 
                        data = use.dat.aOFC,family = 'binomial')
anova(model_choice_1,model_choice_0) # not sig TMS
summary(model_choice_0)
# everything else sig: Valuediff, Didx, base

use.dat.aOFC$fitted_choice <- fitted(model_choice_0, type = "response")

# get unique SubIDs for aOFC dataset
subs <- unique(use.dat.aOFC$SubID)

# storage
results_aOFC <- data.frame(SubID = subs, p_value = NA)

for (i in seq_along(subs)) {
  # filter out one subject
  dat_sub <- filter(use.dat.aOFC, SubID != subs[i])
  
  # fit full and reduced models
  model1 <- glmer(Choice ~ Cond + ValueDiff + Didx + base + (1 | SubID),
                  data = dat_sub, family = "binomial")
  
  model0 <- glmer(Choice ~ ValueDiff + Didx + base + (1 | SubID),
                  data = dat_sub, family = "binomial")
  
  # likelihood ratio test
  lrt <- anova(model1, model0, test = "Chisq")
  
  # store p-value (second row, "Pr(>Chisq)")
  results_aOFC$p_value[i] <- lrt$`Pr(>Chisq)`[2]
}

# order by p-value
results_aOFC <- results_aOFC %>% arrange(p_value)

print(results_aOFC)

# summary stats
summary_stats_aOFC <- results_aOFC %>%
  summarise(
    min_p   = min(p_value, na.rm = TRUE),
    max_p   = max(p_value, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    mean_p   = mean(p_value, na.rm = TRUE)
  )

print(summary_stats_aOFC)



summary_choice_ss_fitted_pOFC = use.dat.pOFC %>%
  group_by(SubID,StimLoc,Cond,Sess) %>%
  reframe(avg_fitted_post_meal=mean(fitted_choice,na.rm = T)) 
summary_choice_ss_fitted_aOFC = use.dat.aOFC %>%
  group_by(SubID,StimLoc,Cond,Sess) %>%
  reframe(avg_fitted_post_meal=mean(fitted_choice,na.rm = T)) 
summary_choice_ss_fitted = rbind(summary_choice_ss_fitted_pOFC,
                                 summary_choice_ss_fitted_aOFC)

p_values <- data.frame(
  StimLoc = c("aOFC", "pOFC"),  
  p_value = c('n.s.', '***'),     
  y_start = c(1, 1),        
  y_end = c(1.03, 1.03), # text y location
  y_low = c(0.98, 0.98)
)

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
         text_x = elem_list_text(color = 'white',face = "bold",size = 16))

c2 = summary_choice_ss_fitted %>%
  ggplot(aes(x=Cond,y=avg_fitted_post_meal)) +
  geom_line(aes(group=SubID), position = pd_dodge, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill = Cond), 
               width = 0.6, fatten = NULL,
               outlier.alpha = 0, 
               alpha = 0.4) +
  stat_summary(fun = mean, geom = "crossbar", 
               fun.min = mean, fun.max = mean,  # makes top = bottom = mean
               width = 0.6,                     # controls horizontal bar width
               color = "black", 
               fatten = 0)  +                     # removes vertical bar ends
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
    aes(
      x = 1.5,         # Position the label between the two groups
      y = y_end,     # Use the pre-calculated y position
      label = p_value),
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

pdf(file.path(FigPaperDir,'Day2_TMS_ChoiceSatedOdor_fitted.pdf'),7,4)
print(c2)
dev.off()


#########################################

# want to get some plots to indicate good model fit of each group
# subject-wise correlation
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

pdf(file.path(FigPaperDir,'Day2_Choice_fitted_actual.pdf'),8,4)
print(bt_sub_corr)
dev.off()

#########################################

# Step 1: Filter out actual values for both conditions
actual_sham <- combined_summary %>%
  filter(Cond == "sham-sham") %>%
  select(SubID, StimLoc, mean_actual) %>%
  rename(mean_actual_sham = mean_actual)

actual_cTBS <- combined_summary %>%
  filter(Cond == "sham-cTBS") %>%
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
  filter(Cond == "sham-cTBS") %>%
  select(SubID, StimLoc, mean_predicted) %>%
  rename(mean_predicted_cTBS = mean_predicted)

predicted_effect_df <- left_join(predicted_sham, predicted_cTBS, by = c("SubID", "StimLoc")) %>%
  mutate(TMS_effect_predicted = mean_predicted_cTBS - mean_predicted_sham)

# Step 4: Merge actual and predicted effect
tms_effect_summary <- left_join(actual_effect_df %>% select(SubID, StimLoc, TMS_effect_actual),
                                predicted_effect_df %>% select(SubID, StimLoc, TMS_effect_predicted),
                                by = c("SubID", "StimLoc"))

# Step 5: Save to CSV
write_csv(tms_effect_summary, file.path(pro_dat_dir,"TMS_effect_summary_day2.csv"))

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
      Cond == "sham-cTBS" ~ 0.6
    ),
    y = case_when(
      Cond == "sham-sham" ~ 0.2,
      Cond == "sham-cTBS" ~ 0.1
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

pdf(file.path(FigPaperDir,'Day2_Choice_roc.pdf'),8,4)
print(roc)
dev.off()

##################################

# put aOFC and pOFC together testing interactions
use.dat = use_choice_dat_ss %>%
                   subset(PrePost=='Post' &
                   Cond %in% c('sham-cTBS','sham-sham')) %>%
                   filter(complete.cases(.)) 
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

# tried plotting use the interaction model
use.dat$fitted_choice <- fitted(model_choice_2, type = "response")
summary_choice_ss_fitted = use.dat %>%
  group_by(SubID,StimLoc,Cond,Sess) %>%
  reframe(avg_fitted_post_meal=mean(fitted_choice,na.rm = T)) 

summary_choice_ss_fitted %>%
  ggplot(aes(x=Cond,y=avg_fitted_post_meal)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(aes(fill=Cond),width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  facet_wrap2(~StimLoc,strip = strip) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) + common 





