
# look at session-wise odor preference baseline
# from ss type of choices
# and from one_odor type of choices
# preference has been saved as subject & condition-wise measure

rm(list = ls())
source('Setup.R')
load(file = '../ProcessedData/pre_choice_df_sated_pref_one_odor.RData')
load(file = '../ProcessedData/pre_choice_df_sated_pref_ss.RData')

pre_choice_df = pre_choice_df %>%
  mutate(base_one_odor = pre_choice_df_sated_pref$pref) %>%
  mutate(base_one_odor = (base_one_odor+1)/2) # rescale from (-1,1) to (0,1)

p1=pre_choice_df %>%
  ggplot(aes(x=base,y=base_one_odor,group=Cond,color=Cond,fill=Cond)) +
  geom_point(aes(color=Cond,fill=Cond),size=1.5,alpha=0.8) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  geom_smooth(method = rlm, se = T, 
              alpha = 0.2, linetype = 3) + common +
  stat_cor(method = 'pearson') +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
  labs(title = "Odor baseline preference", 
       x = 'from savory-sweet choices',
       y = 'from odor-air choices') +
  theme(legend.position = 'topright')

pdf(file.path(FigDir,'OdorBaselinePreferenceCorr.pdf'),4,4)
print(p1)
dev.off()

# based on high correlation of the odor preference
# from different types of odors, compute the weighted average of the two
# ratio of ss choices vs. one-odor choices: 9 vs. 6
# but one-odor choices need to be cut half
pre_choice_df$base_overall = 1/4*pre_choice_df$base_one_odor +
                                3/4*pre_choice_df$base
save(pre_choice_df,
     file = '../ProcessedData/pre_choice_df_sated_pref_overall.RData')
# the saved overall baseline preference is on 0-1 scale


