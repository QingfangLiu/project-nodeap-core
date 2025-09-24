
############################################################
# Sweet–Savory Choice Analysis
#
# PURPOSE
# Analyze sweet–savory (ChoiceType == 3) odor choices to quantify
# TMS effects on choosing the sated odor, and produce summary
# datasets + figures for downstream modeling and the paper.
#
# INPUTS
# - Setup + globals: scripts/utils/Setup.R
# - Trial-level choices (RData):   beh_data_processed/choice_dat.RData
#
# OUTPUT DATASETS (RData, in beh_data_processed/)
# - Summary_Choice_corrected_dat.RData
#     Session-level summary per SubID/Cond/Sess with
#     ChoiceChangeAB = ((post_A + post_B)/2) − pre_A.
# - choice_dat_ss_w_base.RData
#     Trial-level sweet–savory data with a subject×condition
#     baseline column (pre-meal choice rate).
# - pre_choice_df_sated_pref_ss.RData
#     Subject×condition baseline table used for joins.
#
# FIGURES (PDF, in FigPaperDir)
# - Choice_corr_paper.pdf  (3 panels)
#     A: pre-meal (Set A)  vs post-meal (Set A)  — Pearson r + LM fit
#     B: pre-meal (Set A)  vs post-meal (Set B)  — Pearson r + LM fit
#     C: post-meal (Set A) vs post-meal (Set B)  — Pearson r + LM fit
# - choice_sated_day1.pdf
#     Box + jitter of P(Choosing sated) by Cond (sham–sham, cTBS–sham),
#     faceted by StimLoc; Pre vs Post indicated by linetype/shape.
# - choice_sated_day2.pdf
#     Box + jitter of P(Choosing sated) by Cond (sham–sham, sham–cTBS),
#     faceted by StimLoc; Pre vs Post indicated by linetype/shape.
#
# KEY STEPS (brief)
# 1) Filter to sweet–savory trials; define Choice = 1 if chosen odor
#    equals the devalued (sated) odor.
# 2) Aggregate to session level (mean Choice, mean RT).
# 3) Build correlation dataset (pre_A, post_A, post_B) and plot 3 correlations.
# 4) Compute ChoiceChangeAB and save session-level summary.
# 5) Compute subject×condition pre-meal baseline; join back to trials;
#    save trial-level with baseline and the baseline table.
############################################################

rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# Load full choice dataset
load(file = file.path(processed_dir, 'choice_dat.RData'))

# Step 1: Filter for sweet-savory choices (ChoiceType == 3)
choice_dat_ss <- choice_dat %>%
  filter(ChoiceType == 3) %>%
  mutate(
    ChosenOdor = mapvalues(ChosenOdor, from = c(1, 2), to = c('W', 'A')),
    Choice = ifelse(ChosenOdor == Devalued, 1, 0)  # 1 = chose sated odor
  )

# Step 2: Summarize by subject/session
summary_choice_ss <- choice_dat_ss %>%
  group_by(SubID, StimLoc, StimOrder, StimOrder_day1, StimOrder_day2,
           Cond, Sess, PrePost, Set) %>%
  reframe(
    Choice = mean(Choice, na.rm = TRUE),
    rt = mean(rt, na.rm = TRUE)
  )

# Step 3: Extract data for correlation analysis
pre_choice_A  <- subset(summary_choice_ss, Set == 'A' & PrePost == 'Pre')$Choice
post_choice_A <- subset(summary_choice_ss, Set == 'A' & PrePost == 'Post')$Choice
post_choice_B <- subset(summary_choice_ss, Set == 'B' & PrePost == 'Post')$Choice
conds         <- subset(summary_choice_ss, Set == 'A' & PrePost == 'Pre')$Cond

df_pre_post_summary <- data.frame(
  conds   = conds,
  pre_A   = pre_choice_A,
  post_A  = post_choice_A,
  post_B  = post_choice_B
)

# Step 4: Correlation plots
pcorr1 <- ggplot(df_pre_post_summary, aes(x = pre_A, y = post_A)) +
  geom_point(aes(color = conds), size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = TRUE, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(x = 'P(choose sated), pre-meal (set A)', 
       y = 'P(choose sated), post-meal (set A)', title = NULL) +
  common +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

pcorr2 <- ggplot(df_pre_post_summary, aes(x = pre_A, y = post_B)) +
  geom_point(aes(color = conds), size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = TRUE, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(x = 'P(choose sated), pre-meal (set A)', 
       y = 'P(choose sated), post-meal (set B)', title = NULL) +
  common +
  theme(legend.position = 'none')

pcorr3 <- ggplot(df_pre_post_summary, aes(x = post_A, y = post_B)) +
  geom_point(aes(color = conds), size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = TRUE, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(x = 'P(choose sated), post-meal (set A)', 
       y = 'P(choose sated), post-meal (set B)', title = NULL) +
  common +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir, 'Corr_choice_pre_post_meal.pdf'), width = 14, height = 4)
print(ggarrange(pcorr1, pcorr2, pcorr3, ncol = 3, nrow = 1))
dev.off()

# Step 5: Compute average choice change (from pre-meal to post-meal)
avg_changeAB <- (post_choice_A + post_choice_B) / 2 - pre_choice_A

summary_choice_corrected <- summary_choice_ss %>%
  filter(PrePost == 'Pre') %>%
  select(SubID, StimLoc, StimOrder, Cond, Sess, StimOrder_day1, StimOrder_day2) %>%
  mutate(ChoiceChangeAB = avg_changeAB)

save(summary_choice_corrected, 
     file = file.path(processed_dir,'Summary_Choice_corrected_dat.RData'))

# Step 6: Add subject/session-wise baseline (pre-meal bias) to trial-level data
pre_choice_df <- summary_choice_ss %>%
  filter(PrePost == 'Pre') %>%
  group_by(SubID, Cond) %>%
  summarise(base = mean(Choice, na.rm = TRUE), .groups = 'drop')

choice_dat_ss_w_base <- choice_dat_ss %>%
  left_join(pre_choice_df, by = c('SubID', 'Cond'))

save(choice_dat_ss_w_base,
     file = file.path(processed_dir, 'choice_dat_ss_w_base.RData'))
save(pre_choice_df,
     file = file.path(processed_dir, 'pre_choice_df_sated_pref_ss.RData'))

# draw raw plots of the sated effect without any covariates considered

summary_choice_ss_day1 <- choice_dat_ss %>%
  filter(Cond %in% c('sham-sham','cTBS-sham')) %>%
  group_by(SubID, StimLoc, 
           Cond, PrePost) %>%
  reframe(Choice = mean(Choice, na.rm = TRUE))

# Make sure PrePost is ordered
summary_choice_ss_day1$PrePost <- factor(summary_choice_ss_day1$PrePost, levels = c("Pre", "Post"))

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

p_day1 = ggplot(summary_choice_ss_day1, aes(x = Cond, y = Choice)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = F) +
  geom_jitter(aes(shape=PrePost,colour = Cond),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.6) +
  facet_wrap2(~StimLoc, scales = 'free_y', strip = strip) +
  scale_color_manual(values = use.col.conds) +
  labs(x = NULL, y = "P (Choosing sated)", title = NULL) +
  common +
  theme(legend.position = 'none')



summary_choice_ss_day2 <- choice_dat_ss %>%
  filter(Cond %in% c('sham-sham','sham-cTBS')) %>%
  group_by(SubID, StimLoc, 
           Cond, PrePost) %>%
  reframe(Choice = mean(Choice, na.rm = TRUE))

# Make sure PrePost is ordered
summary_choice_ss_day2$PrePost <- factor(summary_choice_ss_day2$PrePost, levels = c("Pre", "Post"))

p_day2 = ggplot(summary_choice_ss_day2, aes(x = Cond, y = Choice)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = F) +
  geom_jitter(aes(shape=PrePost,colour = Cond),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.6) +
  facet_wrap2(~StimLoc, scales = 'free_y', strip = strip) +
  scale_color_manual(values = use.col.conds) +
  labs(x = NULL, y = "P (Choosing sated)", title = NULL) +
  common +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir,'choice_sated_day1.pdf'),8,4)
print(p_day1)
dev.off()

pdf(file.path(FigPaperDir,'choice_sated_day2.pdf'),8,4)
print(p_day2)
dev.off()


########################################################

# Purpose
# Assess the overall devaluation effect on choices using
# subject-level change scores and visualize the effect by
# TMS target (aOFC vs pOFC).
#
# Inputs
# - Setup.R (sourced): defines paths (e.g., `project_folder`,
#   `processed_dir`) and shared objects used in plotting
#   (e.g., `common`, `use.col.ap.ofc`, `pd`) as well as
#   required libraries.
# - file.path(processed_dir, "Summary_Choice_corrected_dat.RData"):
#   loads an object named `summary_choice_corrected`, a data
#   frame with at least:
#     * SubID : subject identifier
#     * StimLoc : "aOFC" or "pOFC" (TMS target)
#     * ChoiceChangeAB : numeric, (post-meal − pre-meal) choice
#       for sated odor, aggregated per subject
#
# What the script does
# 1) Loads precomputed subject-level change scores.
# 2) Runs one-sided Wilcoxon signed-rank tests testing whether
#    median ChoiceChangeAB < 0 within:
#       - pOFC group
#       - aOFC group
#    and also runs a (two-sided) test on all subjects combined.
# 3) Plots the distribution of ChoiceChangeAB by StimLoc 
#
# Outputs
# - Figure: A ggplot object rendered to the active device showing
#   ChoiceChangeAB (post − pre) by StimLoc (aOFC, pOFC).
#   (this figure is not included in the final manuscript)
#
########################################################


# -----------------------------------------------------
# Statistical testing: devaluation effect within each TMS target
# -----------------------------------------------------

# pOFC group only
value_pOFC <- subset(summary_choice_corrected, StimLoc == 'pOFC')$ChoiceChangeAB
wilcox.test(value_pOFC, alternative = 't')  # Significant: p = 0.0044

# aOFC group only
value_aOFC <- subset(summary_choice_corrected, StimLoc == 'aOFC')$ChoiceChangeAB
wilcox.test(value_aOFC, alternative = 't')  # Significant: p = 0.044

# both groups
value <- summary_choice_corrected$ChoiceChangeAB
wilcox.test(value, alternative = 't')  # Significant: p = 0.0005307

# -----------------------------------------------------
# Plot: Day 2 effect (Post - Pre) across aOFC and pOFC
# -----------------------------------------------------

# Create custom background highlight for each x-axis group
muh_grob <- rectGrob(
  x = 1:2, y = 0,
  gp = gpar(color = 'black', fill = use.col.ap.ofc, alpha = 1)
)

# Violin + boxplot of change scores by stimulation site
summary_choice_corrected %>%
  ggplot(aes(x = StimLoc, y = ChoiceChangeAB)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(group = SubID), position = pd, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = 1, y = 0.45, label = "*", size = 4, color = "black") +     # aOFC sig
  annotate("text", x = 2, y = 0.75, label = "***", size = 4, color = "black") +   # pOFC sig
  labs(
    x = NULL,
    y = "Choice of sated odor\n(post - pre-meal)",
    title = NULL
  ) +
  common +
  coord_cartesian(clip = 'off') +
  theme(
    axis.text.x = element_text(margin = margin(t = 10), color = 'white'),
    legend.position = 'inside',
    legend.position.inside = c(0.3, 0.96),
    legend.key.size = unit(0.28, 'cm')
  ) +
  annotation_custom(grob = muh_grob, xmin = 0, xmax = 1, ymin = -1.08, ymax = -0.88)


