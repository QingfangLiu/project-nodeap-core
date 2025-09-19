
# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
rm(list = ls())

project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))  # paths, libs, theme

processed_dir <- file.path(project_folder, "beh_data_processed")

# ------------------------------------------------------------------
# Load odor ratings
# ------------------------------------------------------------------
load(file = file.path(processed_dir, "Odor_ratings_dat.RData"))

# standardize those ratings (z-score) within each sub
dat <- Odor_ratings_dat %>%
  group_by(SubID) %>%
  mutate(Pleasant = (Pleasant - mean(Pleasant, na.rm = TRUE)) / 
           sd(Pleasant, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------------------------

# Compute change in pleasantness (Post - Pre) for each odor per subject
# Create a new data frame `change_df` with only the change scores
pre_plea  <- subset(dat, PrePost == 'Pre-meal')$Pleasant
post_plea <- subset(dat, PrePost == 'Post-meal')$Pleasant
change_plea <- post_plea - pre_plea
change_df <- subset(dat, PrePost == 'Pre-meal') %>%
  mutate(
    PrePost = 'Change',
    Pleasant = change_plea
  )

# Fit the mixed-effects model
model0 <- lmer(Pleasant ~ (1|SubID), data = change_df)
model1 <- lmer(Pleasant ~ IfDevalue + (1|SubID), data = change_df)
anova(model1,model0) 
# collapsing across all 3 sessions, there was a devaluation effect

# ------------------------------------------------------------------

# Selective satiation index per subject/session: (sated − non-sated) Δ
# change of pleasantness rating (sated - non-sated)
# post-pre: more negative value means better selective satiation effect

ND = subset(change_df,IfDevalue=='non-sated')$Pleasant
D = subset(change_df,IfDevalue=='sated')$Pleasant
Didx = D - ND
sel_df = subset(change_df,IfDevalue=='non-sated')
# these cols are meanless now after reducing
sel_df$IfDevalue = NULL
sel_df$PrePost = NULL 
sel_df$Pleasant = NULL
sel_df$IfDevalue = NULL
sel_df$OdorType = NULL
sel_df$Didx = Didx

# Compare the mixed-effects model
model0 <- lmer(Didx ~ (1|SubID), data = sel_df)
model1 <- lmer(Didx ~ Cond + (1|SubID), data = sel_df)
model2 <- lmer(Didx ~ Sess + (1|SubID), data = sel_df)
model3 <- lmer(Didx ~ Devalued + (1|SubID), data = sel_df)
model4 <- lmer(Didx ~ StimLoc + (1|SubID), data = sel_df)
summary(model0)
anova(model0,model1)
anova(model0,model2)
anova(model0,model3)
summary(model3)
anova(model0,model4)
# devaluation effect was not modulated by Cond, Sess, DevaluedOdor, StimLoc

# ------------------------------------------------------------------

# Panel-wise test: is the selective-satiation index (Didx) < 0 ?
# (one-sample t-test per StimLoc × Cond panel)

# Group by StimLoc and Cond, then apply t-test on Didx
t_results <- sel_df %>%
  group_by(StimLoc, Cond) %>%
  summarise(
    n = n(),
    t_result = list(t.test(Didx)),
    .groups = "drop"
  )

# Extract t-test stats into readable format
t_results_clean <- t_results %>%
  rowwise() %>%
  mutate(
    t_stat = t_result$statistic,
    df = t_result$parameter,
    p_value = t_result$p.value,
    mean = t_result$estimate,
    conf_low = t_result$conf.int[1],
    conf_high = t_result$conf.int[2]
  ) %>%
  select(-t_result)

print(t_results_clean)

# ------------------------------------------------------------------

# save this selective satiation index for later use
SelectSate_dat = sel_df
Sate_dat = change_df
save(SelectSate_dat, Sate_dat, 
     file = file.path(processed_dir, "SelectSate_dat.RData"))

# ------------------------------------------------------------------

# 📦 Create customized facet strip styling
strip <- strip_themed(
  background_y = elem_list_rect(fill = use.col.ap.ofc),
  text_x = elem_list_text(color = use.col.conds, face = "bold", size = 16),
  text_y = elem_list_text(color = "white", face = "bold", size = 16)
)

p_values <- data.frame(
  y_start = 10.8,     
  y_text = 11,
  y_low = 10
)

p1=ggplot(Odor_ratings_dat,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = F) +
  geom_jitter(aes(shape=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.6) +
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
      y = y_text,     # Use the pre-calculated y position
      label = '***'
    ),
    inherit.aes = FALSE,
    size = 8
  ) +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(Cond), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  theme(legend.position = 'inside',
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.35,'cm'),
        legend.position.inside = c(0.1,0.1)) +
  theme(
    strip.text.x.top = element_text(face = "bold"),
    strip.text.y.left = element_text(face = "bold")
  )

pdf(file.path(FigPaperDir,'PleasantRatings_by_cond.pdf'),8,6)
print(p1)
dev.off()
