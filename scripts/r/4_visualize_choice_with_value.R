
# this code visualize choice with different value differences

rm(list = ls())
source('Setup.R')

load(file.path(pro_dat_dir,'choice_dat_ss_w_base_values.RData'))

# to round the learned values
target_values <- seq(0.5, 1, by = 0.05)
target_values
use_choice_dat_ss = use_choice_dat_ss %>%
  subset(SubID != 'NODEAP_17') %>%
  mutate(Rounded_ValueSated = as.numeric(sapply(ValueSated, function(x) 
    target_values[which.min(abs(target_values - x))]))) %>%
  mutate(Rounded_ValueNonSated = as.numeric(lapply(ValueNonSated, function(x) 
    target_values[which.min(abs(target_values - x))]))) %>%
  filter(complete.cases(Choice))
  
# get a summary across subjects and conditions
summary_use_choice_dat_ss = use_choice_dat_ss %>%
  dplyr::group_by(Rounded_ValueSated,Rounded_ValueNonSated) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),
          nsample=n()) 

summary_use_choice_dat_ss_vdiff = use_choice_dat_ss %>%
  mutate(value_diff = Rounded_ValueSated - Rounded_ValueNonSated) %>%
  mutate(value_diff = factor(value_diff)) %>%
  dplyr::group_by(value_diff) %>%
  reframe(MeanChoice=mean(Choice,na.rm=T),
          nsample=n()) %>%
  mutate(value_diff = varhandle::unfactor(value_diff))

common_setting = 
  theme(axis.text.x = element_text(colour = "black"), # axis color
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5), # center title
        panel.background = element_rect(fill = "gray90"),  # Lighten the panel background
        panel.grid = element_line(color = "gray85"), # Adjust grid line colors
        #plot.margin = margin(t = 10, r = 5.5, b = 5.5, l = 5.5, unit = 'pt'),
        text = element_text(size=16)) # text size & font

use_xlab = expression(w[SA])
use_ylab = expression(w[NS])

p_choice_1 = ggplot(summary_use_choice_dat_ss,
       aes(x = Rounded_ValueSated, 
           y = Rounded_ValueNonSated, 
           size = log(nsample),  # Circle size indicates log sample size
           color = MeanChoice)) +  # Circle color indicates MeanChoice
  geom_point(shape = 16, alpha = 0.8) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = expression(w[SA]), y = expression(w[NS]),
       size = expression(log(n[trial]))) +
  scale_color_gradient2(low = "red", high = "green", 
                        mid = 'white', midpoint = 0.5,
                        name = "P(SA)") +
  theme(legend.position = "right") + common_setting

p_choice_2 = summary_use_choice_dat_ss_vdiff %>%
  ggplot(aes(x=value_diff,
             y=MeanChoice)) +
  geom_point(aes(size = log(nsample)),
             shape = 16, alpha = 0.6) +
  geom_smooth(method = lm, se = TRUE, color = 'black',
              linetype = 2) + 
  scale_size_continuous(range = c(2, 12)) +  # Adjust circle size range
  labs(x = expression(w[SA]-w[NS]), 
       size = expression(log(n[trial])),
       y = 'P(SA)') +
  common_setting + 
  theme_classic(base_size = 16)

cor.test(summary_use_choice_dat_ss_vdiff$value_diff, 
    summary_use_choice_dat_ss_vdiff$MeanChoice)

pdf(file.path(FigPaperDir,'Choice_w_values_agg_gradient.pdf'),5,4)
print(p_choice_1)
dev.off()

pdf(file.path(FigPaperDir,'Choice_w_values_agg.pdf'),5,4)
print(p_choice_2)
dev.off()



