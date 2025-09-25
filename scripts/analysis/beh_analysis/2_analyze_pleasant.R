
# ======================================================================
# Purpose:
#   - Standardize odor pleasantness ratings within subject
#   - Compute change scores (Post − Pre)
#   - Estimate devaluation effect via mixed-effects models
#   - Derive selective satiation index (Didx = Δ_sated − Δ_non-sated)
#   - Test Didx < 0 per (StimLoc × Cond) panel and convert p to stars
#   - Save processed data and generate facet plot with significance marks
#
# Inputs:
#   - processed_dir/Odor_ratings_dat.RData  (expects Odor_ratings_dat)
#
# Outputs:
#   - processed_dir/SelectSate_dat.RData  (SelectSate_dat, Sate_dat)
#   - Figs_paper/PleasantRatings_by_cond.pdf
#
# Notes:
#   - Plot uses raw ratings for visualization (Odor_ratings_dat),
#     while modeling uses standardized values (z-scored within SubID).
# ======================================================================

rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# ------------------------------------------------------------------
# Load odor ratings
# ------------------------------------------------------------------
load(file = file.path(processed_dir, "Odor_ratings_dat.RData"))

# standardize those ratings (z-score) within each sub
dat <- Odor_ratings_dat %>%
  group_by(SubID) %>%
  mutate(Pleasant = as.numeric(scale(Pleasant))) %>%  # z-score within SubID
  ungroup()

# ------------------------------------------------------------------

# Compute change in pleasantness (Post - Pre) for each odor per subject
# Create a new data frame `change_df` with only the change scores
pre_plea  <- subset(dat, PrePost == 'Pre-meal')$Pleasant
post_plea <- subset(dat, PrePost == 'Post-meal')$Pleasant
change_plea <- post_plea - pre_plea
change_df <- subset(dat, PrePost == 'Pre-meal') %>%
  mutate(PrePost = 'Change', Pleasant = change_plea)

# ----------------------------------------------------------------------
# Devaluation effect (collapsed across sessions)
# Mixed-effects model: Pleasant ~ IfDevalue + (1|SubID)
# ----------------------------------------------------------------------

model0 <- lmer(Pleasant ~ (1|SubID), data = change_df)
model1 <- lmer(Pleasant ~ IfDevalue + (1|SubID), data = change_df)
anova(model1,model0) 
# collapsing across all 3 sessions, there was a devaluation effect

# ------------------------------------------------------------------

# ----------------------------------------------------------------------
#    Selective satiation index (Didx = Δ_sated − Δ_non-sated)
#    more negative value means better selective satiation effect
#    Then test modulation by Cond / Sess / Devalued / StimLoc
# ----------------------------------------------------------------------

ND = subset(change_df,IfDevalue=='non-sated')$Pleasant
D = subset(change_df,IfDevalue=='sated')$Pleasant
Didx = D - ND
sel_df = subset(change_df,IfDevalue=='non-sated')
# drop columns no longer meaningful after reduction
sel_df$IfDevalue = NULL
sel_df$PrePost = NULL 
sel_df$Pleasant = NULL
sel_df$IfDevalue = NULL
sel_df$OdorType = NULL
sel_df$Didx = Didx

# Compare mixed-effects models
model0 <- lmer(Didx ~ (1|SubID), data = sel_df)
model1 <- lmer(Didx ~ Cond + (1|SubID), data = sel_df)
model2 <- lmer(Didx ~ Sess + (1|SubID), data = sel_df)
model3 <- lmer(Didx ~ Devalued + (1|SubID), data = sel_df)
model4 <- lmer(Didx ~ StimLoc + (1|SubID), data = sel_df)

anova(model0,model1)
anova(model0,model2)
anova(model0,model3)
anova(model0,model4)
# Comment: no significant modulation by Cond, Sess, DevaluedOdor, or StimLoc (all p > .05)

# ----------------------------------------------------------------------
#    Panel-wise tests: Is Didx < 0 ? (one-sample t-tests per StimLoc×Cond)
#    + p → stars helper + panel-wise annotations table (p_annots)
# ----------------------------------------------------------------------

t_results <- sel_df %>%
  group_by(StimLoc, Cond) %>%
  reframe(
    n = n(),
    t_result = list(t.test(Didx, mu = 0, alternative = "less")))

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

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}

p_annots <- t_results_clean %>%
  mutate(stars = vapply(p_value, p_to_stars, character(1))) %>%
  select(StimLoc, Cond, stars)


# ------------------------------------------------------------------

# save this selective satiation index for later use
SelectSate_dat = sel_df
Sate_dat = change_df
save(SelectSate_dat, Sate_dat, 
     file = file.path(processed_dir, "SelectSate_dat.RData"))

# ------------------------------------------------------------------

# Create customized facet strip styling
strip <- strip_themed(
  background_y = elem_list_rect(fill = use.col.ap.ofc),
  text_x = elem_list_text(color = use.col.conds, size = 16),
  text_y = elem_list_text(color = "white", face = "bold", size = 16)
)

# fixed heights
p_values <- tibble::tibble(
  y_start = 10.8,
  y_text  = 11,
  y_low   = 10
)

# combine: add fixed heights to every StimLoc×Cond row
p_values_df <- tidyr::crossing(p_annots, p_values) %>%
  mutate(
    # make ns bigger/smaller if needed
    size_label = ifelse(stars == "n.s.", 8, 8),
    # shift n.s. up a bit
    y_text     = ifelse(stars == "n.s.", 11.8, y_text)
  )

# ----------------------------------------------------------------------
#    Plot: Odor pleasantness by IfDevalue × PrePost, faceted by StimLoc × Cond
#    (Uses raw Odor_ratings_dat for visualization; stats from standardized data)
# ----------------------------------------------------------------------

p1 <- ggplot(Odor_ratings_dat, aes(x = IfDevalue, y = Pleasant)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(aes(shape = PrePost),
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2),
              size = 1.5, alpha = 0.6) +
  geom_segment(
    data = p_values_df, aes(x = 1, xend = 2, y = y_start, yend = y_start),
    inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(
    data = p_values_df, aes(x = 1, xend = 1, y = y_low, yend = y_start),
    inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(
    data = p_values_df, aes(x = 2, xend = 2, y = y_low, yend = y_start),
    inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_text(
    data = p_values_df,
    aes(x = 1.5, y = y_text, label = stars, size = size_label),
    inherit.aes = FALSE, show.legend = F) +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  facet_grid2(rows = vars(StimLoc),
              cols = vars(Cond),
              scales = "fixed",
              strip = strip,
              axes = "all") +
  theme(legend.position = "inside",
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.position.inside = c(0.1, 0.1)) +
  theme(
    strip.text.x.top  = element_text(face = "bold"),
    strip.text.y.left = element_text(face = "bold")
  )

pdf(file.path(FigPaperDir,'PleasantRatings_by_cond.pdf'),8,6)
print(p1)
dev.off()


