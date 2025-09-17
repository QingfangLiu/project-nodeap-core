
#######################################################
# The analysis is based on
# subject-level summary data: average post-meal choice (Set A & B)
# minus pre-meal choice (Set A). 
#
# Main steps:
# - Load summary data with per-subject choice changes
# - Use one-sided Wilcoxon tests to assess overall devaluation effects
# - Plot choice change separately for aOFC and pOFC stimulation sites
#######################################################

rm(list = ls())
source('Setup.R')  # Load paths, libraries, and custom settings

# Load precomputed subject-level choice change data
load(file = file.path(pro_dat_dir, 'Summary_Choice_corrected_dat.RData'))

# -----------------------------------------------------
# Statistical testing: devaluation effect within each TMS target
# -----------------------------------------------------

# pOFC: test if post-pre choice change is less than 0 (stronger devaluation)
value_pOFC <- subset(summary_choice_corrected, StimLoc == 'pOFC')$ChoiceChangeAB
wilcox.test(value_pOFC, alternative = 'less')  # Significant: p = 0.0023

# aOFC: same test
value_aOFC <- subset(summary_choice_corrected, StimLoc == 'aOFC')$ChoiceChangeAB
wilcox.test(value_aOFC, alternative = 'less')  # Significant: p = 0.024

# both groups
value <- summary_choice_corrected$ChoiceChangeAB
wilcox.test(value, alternative = 't')  

# -----------------------------------------------------
# Plot: Day 2 effect (Post - Pre) across aOFC and pOFC
# -----------------------------------------------------

# Create custom background highlight for each x-axis group
muh_grob <- rectGrob(
  x = 1:2, y = 0,
  gp = gpar(color = 'black', fill = use.col.ap.ofc, alpha = 1)
)

# Violin + boxplot of change scores by stimulation site
s1 <- summary_choice_corrected %>%
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

# Save the plot
pdf(file.path(FigPaperDir, 'Day2_overall_ChoiceSatedOdor.pdf'), width = 4, height = 4)
print(s1)
dev.off()
