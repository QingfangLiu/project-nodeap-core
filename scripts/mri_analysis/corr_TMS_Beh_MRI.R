
# =============================================================================
# Script: Corr_TMS_Beh_MRI.R
# Purpose: Correlate behavioral TMS effects with neural TMS effects derived
#          from cVAE embeddings; and visualize/group-level distance differences.
#
# Inputs (from `project-nodeap-core`):
#   - data_mri_processed/merged_day1_tms_effect_0.csv
#       • Columns used: TMS_effect_predicted, day1_tms_effect_mri, StimLoc
#   - data_mri_processed/merged_day2_tms_effect_0.csv
#       • Columns used: TMS_effect_predicted, day2_tms_effect_mri, StimLoc
#   - data_mri_processed/dist_summary_roi_0.csv
#       • Columns used: d_null_real, d_null_sham, subject
#
# Outputs:
#   - figures/Corr_TMS_beh_mri_day1.pdf      (Day 1 correlation by StimLoc)
#   - figures/Corr_TMS_beh_mri_day2.pdf      (Day 2 correlation by StimLoc)
#   - figures/TMS_mri_dist_comp.pdf          (Distance-to-null comparison)
#
# =============================================================================

rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# Load the data
df_day1 = read.csv(file.path(project_folder,"data_mri_processed","merged_day1_tms_effect_0.csv"))
df_day2 = read.csv(file.path(project_folder,"data_mri_processed","merged_day2_tms_effect_0.csv"))
df <- read.csv(file.path(project_folder,"data_mri_processed","dist_summary_roi_0.csv"))  

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',face = "bold",size = 16))

p_day1 = ggscatter(df_day1, 
          x = "TMS_effect_predicted", 
          y = "day1_tms_effect_mri", 
          add = "reg.line", 
          conf.int = TRUE,
          cor.coef = TRUE, 
          cor.method = "spearman") +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_y_continuous(limits = c(-2, 3.2)) +
  scale_x_continuous(limits = c(-0.5, 0.8)) +
  labs(
    title = NULL,
    x = "Behavioral TMS effect",
    y = "Neural TMS effect (aOFC seed)"
  ) +
  common +
  theme(legend.position = "none")

pdf(file.path(FigPaperDir,'Corr_TMS_beh_mri_day1.pdf'),6.5,4)
print(p_day1)
dev.off()


p_day2 = ggscatter(df_day2, 
          x = "TMS_effect_predicted", 
          y = "day2_tms_effect_mri", 
          add = "reg.line", 
          conf.int = TRUE,
          cor.coef = TRUE, 
          cor.method = "spearman") +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip) +
  scale_y_continuous(limits = c(-2, 3.2)) +
  scale_x_continuous(limits = c(-0.4, 0.5)) +
  labs(
    title = NULL,
    x = "Behavioral TMS effect",
    y = "Neural TMS effect (aOFC seed)"
  ) +
  common +
  theme(legend.position = "none")

print(p_day2)

############# plot avg. distance difference 

t.test(df$d_null_real, df$d_null_sham, paired = TRUE)

df_long <- df %>%
  pivot_longer(cols = c(d_null_sham, d_null_real),
               names_to = "Condition",
               values_to = "Distance")

p_values <- data.frame(
  y_start = 5.3,        
  y_end = 5.4,
  y_low = 5.2
)

p_comp_dist = ggplot(df_long, aes(x = Condition, y = Distance)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "gray90", color = "black") +
  geom_line(aes(group = subject), color = "gray60", alpha = 0.5, position = pd, linewidth = 0.5) +
  geom_jitter(aes(group = subject), position = pd, size = 2, alpha = 0.8) +
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
      label = '*'
    ),
    inherit.aes = FALSE,
    size = 8
  ) +
  scale_x_discrete(labels = c("d_null_real" = "cTBS", "d_null_sham" = "sham")) +
  labs(title = NULL, y = "Distance to null", x = NULL) + common

pdf(file.path(FigPaperDir,'TMS_mri_dist_comp.pdf'),3,3.5)
print(p_comp_dist)
dev.off()

