

############################################################
# Conditioning Task — Learning Across Runs and Day-1 TMS
#
# INPUT
# - beh_data_processed/Conditioning.RData   (trial-level: OdorChosen, Run, Sess,
#   SubID, StimLoc, Cond_day1, StimOrder_day1, etc.) loaded via Setup.R.
#
# WHAT THIS SCRIPT DOES
# - Plots learning (P[choose odor-predicting stim]) across runs:
#     • overall by StimLoc
#     • by Day-1 condition (cTBS vs sham)
#     • by session (Sess)
# - Fits GLMMs (binomial) testing effects of Run, Sess, Day-1 TMS (Cond_day1),
#   and interactions; repeats models within aOFC / pOFC.
# - Compares cTBS vs sham at Run 5 with paired t-tests within StimLoc.
#
# OUTPUTS (PDFs in FigPaperDir)
# - Conditioning_byStimOrder.pdf     (by StimLoc × StimOrder_day1)
# - Conditioning_day1.pdf            (by StimLoc × Cond_day1)

############################################################

rm(list = ls())

# Setup & paths -----------------------------------------------------------
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

processed_dir <- file.path(project_folder, "beh_data_processed")

# Data -------------------------------------------------------------------
load(file.path(processed_dir, "Conditioning.RData"))
# expected: conditioning_dat with columns:
# SubID, Run, Sess, StimLoc (aOFC/pOFC), Cond_day1 (cTBS/sham),
# StimOrder_day1, OdorChosen (0/1)

# -----------------------------------------------------------------------------
# Day-1 discrimination learning 
# Test: Does choice of the odor-predictive stimulus increase over RUN?
# Model: trial-level GLMM with a fixed effect of Run and a random intercept
#        for subjects. Data include both StimLoc groups (aOFC & pOFC).

model_conditioning_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_0 <- glmer(OdorChosen ~ (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_1,model_conditioning_0)
summary(model_conditioning_1)

# Interpretation: A significant positive Run effect (and LRT m1 vs m0) supports
#                 successful discrimination learning across the five blocks.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Day-1 TMS influence + last-block (Run 5) comparison — matches paper text
#
# “Although this increase was influenced by TMS applied prior to the discrimination
# task (cTBS vs. sham; p = 1.27 × 10^-7), there was no effect of cTBS on choices of
# odor-predicting stimuli in the last block of the discrimination task in either group
# (aOFC: p = 0.605; pOFC: p = 0.967, t-test).”
#
# Below: aggregate Run 5 only, compute subject means per Cond_day1 within each StimLoc,
#        and compare cTBS vs sham separately for aOFC and pOFC with t-tests.
# -----------------------------------------------------------------------------

conditioning_dat_run5 <- subset(conditioning_dat, Run == 5)

library(dplyr)
library(tidyr)

# Step 1 — restrict to Run 5
run5_data <- conditioning_dat_run5 %>%
  filter(Run == 5)

# Step 2 — subject mean choice per Day-1 condition (cTBS / sham)
agg_data <- run5_data %>%
  group_by(SubID, StimLoc, Cond_day1) %>%
  summarise(mean_choice = mean(OdorChosen), .groups = "drop")

# Step 3 — split by StimLoc and widen to cTBS / sham columns (for tests)
agg_aofc <- agg_data %>% filter(StimLoc == "aOFC") %>%
  pivot_wider(names_from = Cond_day1, values_from = mean_choice)
agg_pofc <- agg_data %>% filter(StimLoc == "pOFC") %>%
  pivot_wider(names_from = Cond_day1, values_from = mean_choice)

# Step 4 — last-block comparisons reported in the paper (Figure 3A text)
t_aofc <- t.test(agg_aofc$cTBS, agg_aofc$sham, paired = TRUE)  # aOFC: p ≈ 0.605
t_pofc <- t.test(agg_pofc$cTBS, agg_pofc$sham, paired = TRUE)  # pOFC: p ≈ 0.967

cat("Run-5 comparison (aOFC):\n");  print(t_aofc)
cat("\nRun-5 comparison (pOFC):\n"); print(t_pofc)


# -----------------------------------------------------------------------------
# Day-1 discrimination: overall TMS (cTBS vs sham), session effect, interaction
# Matches paper text:
#  - Significant main effect of TMS
#  - Significant main effect of session number (1st, 2nd, 3rd)
#  - Significant Session × TMS interaction
# Notes:
#  - We control for learning across blocks with Run as a covariate.
#  - All inferential p-values are the Wald tests shown in summary(model_final).
# -----------------------------------------------------------------------------

model_conditioning1_4 <- glmer(OdorChosen ~ Run + Sess + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_5 <- glmer(OdorChosen ~ Run + Sess + Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_6 <- glmer(OdorChosen ~ Run + Sess * Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')

anova(model_conditioning1_1,model_conditioning1_4)
# adding sess helps when there is only run
anova(model_conditioning1_5,model_conditioning1_4)
# adding day 1 TMS helps when there are runs and sess
anova(model_conditioning1_5,model_conditioning1_6)
# interaction b/t sess and day 1 TMS conditions

summary(model_conditioning1_6)


####### plotting below #########


# plot overall learning trend of each aOFC and pOFC group

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

subject_cond_means <- conditioning_dat %>%
  group_by(SubID, Run, StimLoc, Cond_day1) %>%
  summarise(Choice = mean(OdorChosen, na.rm = TRUE), .groups = "drop")

disc_overall_day1 = ggplot(subject_cond_means, aes(x = Run, y = Choice, color = Cond_day1)) +
  facet_wrap2(~StimLoc, scales = 'free_y', strip = strip) +  # Same layout as before
  geom_line(aes(group = interaction(SubID, Cond_day1)), alpha = 0.3, linewidth = 0.3) +  # individual subjects
  stat_summary(aes(group = Cond_day1), fun = mean, geom = "line", linewidth = 1.5) +  # group-level mean per condition
  stat_summary(aes(group = Cond_day1), fun.data = mean_se, geom = "errorbar",
               width = 0.3, linewidth = 1) +  # SE bars
  coord_cartesian(ylim = c(0.4, 1)) +
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "Blocks", 
       y = "Percentage of selecting \n odor-predicting stimulus", color = "Condition", title = NULL) +
  common  

# Figure 3A
pdf(file.path(FigPaperDir,'Conditioning_day1.pdf'),8,4)
print(disc_overall_day1)
dev.off()


# Figure: By StimLoc × StimOrder_day1 (same subjects per panel) --------

custom_labeller <- labeller(
  StimOrder_day1 = label_value
)

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
out2 = conditioning_dat %>%
  group_by(SubID,Run,Sess,StimOrder_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F, width = 0.5,
               geom = "errorbar", linewidth = 0.8) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(StimOrder_day1), 
              scales = 'fixed',
              axes = 'all',
              strip = strip,
              labeller = custom_labeller) +
  scale_color_manual(values = use.col.sess) +
  labs(x = "Blocks", y = "P(Choosing rewarding stim)", fill = 'Session',
       color = 'Session', linetype = 'Session') +
  common +
  theme(legend.position = 'inside',
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.7,'cm'),
        legend.position.inside = c(0.9,0.8))

# Figure S1B
pdf(file.path(FigDir,'Conditioning_byStimOrder.pdf'),9,6)
print(out2)
dev.off()


