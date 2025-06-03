
# This script analyzes sweet-savory odor choices to assess the effect of TMS on decision-making.
# 
# Main steps:
# 1. Load and filter the raw dataset to include only sweet-savory (ChoiceType == 3) trials.
# 2. Create a binary choice outcome indicating selection of the sated odor.
# 3. Summarize choices at the session level (by subject, condition, and timepoint).
# 4. Visualize pre/post-meal choice correlations across conditions using ggplot.
# 5. Compute an average choice change score (post - pre) for Set A and Set B trials.
# 6. Generate a cleaned summary dataset for downstream statistical modeling.
# 7. Compute and integrate subject-level baseline preference into trial-level data for modeling.
# 
# Outputs:
# - Session-level summary with choice change scores.
# - Trial-level dataset with baseline bias.
# - Correlation plots saved to PDF.


rm(list = ls())
source('Setup.R')  # Load paths, libraries, and any global settings

# Load full choice dataset
load(file = file.path(pro_dat_dir, 'choice_dat.RData'))

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
lab1 <- 'pre-meal (set A)'
lab2 <- 'post-meal (set A)'
lab3 <- 'post-meal (set B)'

pcorr1 <- ggplot(df_pre_post_summary, aes(x = pre_A, y = post_A)) +
  geom_point(aes(color = conds), size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = TRUE, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(x = lab1, y = lab2, title = NULL) +
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
  labs(x = lab1, y = lab3, title = NULL) +
  common +
  theme(legend.position = 'none')

pcorr3 <- ggplot(df_pre_post_summary, aes(x = post_A, y = post_B)) +
  geom_point(aes(color = conds), size = 2) +
  scale_color_manual(values = use.col.conds) +
  geom_smooth(method = "lm", se = TRUE, color = 'black') +
  stat_cor(method = 'pearson', size = 5) +
  labs(x = lab2, y = lab3, title = NULL) +
  common +
  theme(legend.position = 'none')

pdf(file.path(FigPaperDir, 'Choice_corr_paper.pdf'), width = 8, height = 8)
ggarrange(pcorr1, pcorr2, pcorr3, ncol = 2, nrow = 2)
dev.off()

# Step 5: Compute average choice change (from pre-meal to post-meal)
avg_changeAB <- (post_choice_A + post_choice_B) / 2 - pre_choice_A

summary_choice_corrected <- summary_choice_ss %>%
  filter(PrePost == 'Pre') %>%
  select(SubID, StimLoc, StimOrder, Cond, Sess, StimOrder_day1, StimOrder_day2) %>%
  mutate(ChoiceChangeAB = avg_changeAB)

save(summary_choice_corrected, file = '../ProcessedData/Summary_Choice_corrected_dat.RData')

# Step 6: Add subject/session-wise baseline (pre-meal bias) to trial-level data
pre_choice_df <- summary_choice_ss %>%
  filter(PrePost == 'Pre') %>%
  group_by(SubID, Cond) %>%
  summarise(base = mean(Choice, na.rm = TRUE), .groups = 'drop')

choice_dat_ss_w_base <- choice_dat_ss %>%
  left_join(pre_choice_df, by = c('SubID', 'Cond'))

save(choice_dat_ss_w_base, file = '../ProcessedData/choice_dat_ss_w_base.RData')
save(pre_choice_df, file = '../ProcessedData/pre_choice_df_sated_pref_ss.RData')
