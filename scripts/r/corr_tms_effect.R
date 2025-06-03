

rm(list = ls())
source('Setup.R')


# plot correlation b/ neural tms effect and beh tms effect

df_day1 = read.csv(file.path(pro_mri_dir,"merged_day1_tms_effect_0.csv"))
strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = 'white',face = "bold",size = 16))

p_day1 = ggscatter(df_day1, 
          x = "TMS_effect_predicted", 
          y = "day1_tms_effect", 
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

pdf(file.path(FigPaperDir,'Corr_TMS_beh_mri_day1.pdf'),7,4)
print(p_day1)
dev.off()


df_day2 = read.csv(file.path(pro_mri_dir,"merged_day2_tms_effect_0.csv"))
p_day2 = ggscatter(df_day2, 
          x = "TMS_effect_predicted", 
          y = "day2_tms_effect", 
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

pdf(file.path(FigPaperDir,'Corr_TMS_beh_mri_day2.pdf'),7,4)
print(p_day2)
dev.off()
