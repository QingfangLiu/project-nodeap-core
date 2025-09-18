

# plot tms effect
# stats come from mixed effect models including value diff and Didx 
# in another code


rm(list = ls())
source('Setup.R')  # Load paths, libraries, and custom settings

# Load full choice dataset
load(file = file.path(pro_dat_dir, 'Summary_Choice_corrected_dat.RData'))


p_values <- data.frame(
  StimLoc = c("aOFC", "pOFC"),  
  p_value = c(0.655, 0.00034),     
  y_start = c(0.74, 0.74),        
  y_end = c(0.81, 0.81),
  y_low = c(0.72, 0.72)
)

strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
c1=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, 
              size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  facet_wrap2(~StimLoc,scales = 'free',
              strip = strip,
              axes = 'all') +
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
      label = paste0("p = ", p_value)
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  annotate("text",x=1.5,y=0.9,label='Mixed-effects logistic model:') +
  scale_fill_manual(values = use.col.conds) +
  scale_color_manual(values = use.col.conds) +
  coord_cartesian(ylim = c(-0.8,0.95)) +
  labs(x = NULL, title = NULL,
       y = "Choice of sated odor\n (post - pre-meal)") + common +
  theme(legend.position = "none")

custom_labeller <- labeller(
  StimOrder_day2 = label_value
)

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
# based on c1
# separated by diff stim orders of Day 2 TMS
c2=summary_choice_corrected %>%
  subset(Cond %in% c('sham-sham','sham-cTBS')) %>%
  ggplot(aes(x=Cond,y=ChoiceChangeAB,fill=Cond)) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Cond,group=SubID), position = pd, size = 2, alpha = 0.9) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(StimOrder_day2), 
              scales = 'fixed',
              axes = 'all',
              strip = strip,
              labeller = custom_labeller) +
  scale_fill_manual(values = use.col.conds) +
  scale_color_manual(values = use.col.conds) +
  labs(x = NULL, title = NULL, y = "Choice of sated odor\n (post - pre-meal)") + 
  common + theme(legend.position = "none")

pdf(file.path(FigPaperDir,'Day2_TMS_ChoiceSatedOdor_Changes.pdf'),7,4)
print(c1)
dev.off()

pdf(file.path(FigPaperDir,'Day2_TMS_ChoiceSatedOdor_Changes_by_order.pdf'),11,8)
print(c2)
dev.off()

