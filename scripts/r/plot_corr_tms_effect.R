

rm(list = ls())
source('Setup.R')

# plot correlation b/ neural tms effect and beh tms effect
# also the conn index derived from roi-based FC (did not find anything)

df_day1 = read.csv(file.path(pro_mri_dir,"merged_day1_tms_effect_0.csv"))
ConnIdx = read.csv("../../data_mri_processed/ConnIdx.csv")

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


df_day2 = read.csv(file.path(pro_mri_dir,"merged_day2_tms_effect_0.csv"))
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

pdf(file.path(FigPaperDir,'Corr_TMS_beh_mri_day2.pdf'),6,4)
print(p_day2)
dev.off()



########## corr with conn idx

df_day1_more <- df_day1 %>%
  left_join(ConnIdx %>% select(SubID, FC), by = "SubID")
ggscatter(df_day1_more, 
                   x = "day1_tms_effect", 
                   y = "FC", 
                   add = "reg.line", 
                   conf.int = TRUE,
                   cor.coef = TRUE, 
                   cor.method = "spearman") +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip)

df_day2_more <- df_day2 %>%
  left_join(ConnIdx %>% select(SubID, FC), by = "SubID")
ggscatter(df_day2_more, 
          x = "TMS_effect_predicted", 
          y = "FC", 
          add = "reg.line", 
          conf.int = TRUE,
          cor.coef = TRUE, 
          cor.method = "spearman") +
  facet_wrap2(~StimLoc,scales = 'free',strip = strip)

# nothing from here

############# plot avg. distance difference 

# Load the data
df <- read.csv(file.path(pro_mri_dir,"dist_summary_roi_0.csv"))  

df <- df %>%
  left_join(
    ConnIdx %>% select(SubID, StimLoc) %>% distinct(),
    by = c("subject" = "SubID")
  )

# Perform paired t-test
t_test_result <- t.test(df$d_null_real, df$d_null_sham, paired = TRUE)

# Print results
print(t_test_result)

# Reshape to long format correctly
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


# Reshape to long format for ANOVA
df_long <- df %>%
  pivot_longer(
    cols = c(d_null_sham, d_null_real),
    names_to = "Condition",
    values_to = "Distance"
  )

library(ez)  
# Run mixed (within-subject) ANOVA
anova_result <- ezANOVA(
  data = df_long,
  dv = Distance,
  wid = subject,
  within = .(Condition),
  between = .(StimLoc),
  detailed = TRUE
)

anova_result

