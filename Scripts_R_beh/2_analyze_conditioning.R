

rm(list = ls())
source('Setup.R')

###############################
load('../ProcessedData/Conditioning.RData')

# this one, each panel has the same group of people
# separated by stimulation order

custom_labeller <- labeller(
  StimOrder_day1 = label_value
)

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))
out2 = conditioning_dat %>%
  group_by(SubID,Run,Sess,StimOrder_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F, width = 0.5,
               geom = "errorbar", linewidth = 0.8) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(StimOrder_day1), 
              scales = 'fixed',
              axes = 'all',
              strip = strip,
              labeller = custom_labeller) +
  scale_color_manual(values = use.col.sess) +
  labs(x = "Runs", y = "P(Choosing rewarding stim)", fill = 'Session',
       color = 'Session', linetype = 'Session') +
  common +
  theme(legend.position = 'inside',
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.7,'cm'),
        legend.position.inside = c(0.9,0.8))

pdf(file.path(FigDir,'Conditioning_paper1.pdf'),9,6)
print(out2)
dev.off()



model_conditioning_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_0 <- glmer(OdorChosen ~ (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_1,model_conditioning_0)
summary(model_conditioning_1)
# yes conditioning improves over runs


####################################
# focus on cond_day1: cTBS vs sham
# so cTBS only includes cTBS-sham; while
# sham include both sham-sham and sham-cTBS
####################################

model_conditioning1_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_2 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_3 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')

anova(model_conditioning1_2,model_conditioning1_3)
anova(model_conditioning1_1,model_conditioning1_2)
summary(model_conditioning1_2)
summary(model_conditioning1_3) 

model_conditioning1_4 <- glmer(OdorChosen ~ Run + Sess + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_5 <- glmer(OdorChosen ~ Run + Sess + Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
model_conditioning1_6 <- glmer(OdorChosen ~ Run + Sess * Cond_day1 + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')

anova(model_conditioning1_1,model_conditioning1_4)
# adding sess helps when there is only run
anova(model_conditioning1_5,model_conditioning1_4)
# adding day 1 TMS helps when there are runs and sess
anova(model_conditioning1_5,model_conditioning1_6)
# interaction b/t sess and day 1 TMS conditions
summary(model_conditioning1_6)



model_conditioning1_4 <- glmer(OdorChosen ~ Run * Cond_day1 + StimLoc + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_4,model_conditioning1_3) # need StimLoc? p=0.055



# aOFC only
use.dat = subset(conditioning_dat,StimLoc=='aOFC')
use.dat = use.dat[complete.cases(use.dat$OdorChosen), ] # remove NA values

model_conditioning_aOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_aOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_aOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')


anova(model_conditioning_aOFC_0,model_conditioning_aOFC_1) # day1 TMS main effect
anova(model_conditioning_aOFC_1,model_conditioning_aOFC_2) # no interaction
summary(model_conditioning_aOFC_1)
# describe: faster learning with runs under sham for aOFC stimulation group




# pOFC only
use.dat = subset(conditioning_dat,StimLoc=='pOFC')
use.dat = use.dat[complete.cases(use.dat$OdorChosen), ] # remove NA values

model_conditioning_pOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')

anova(model_conditioning_pOFC_0,model_conditioning_pOFC_1) # day1 TMS main effect
anova(model_conditioning_pOFC_1,model_conditioning_pOFC_2) # interaction b/t TMS & runs
summary(model_conditioning_pOFC_2)





# do some testing here
use_dat = subset(conditioning_dat,StimOrder_day1=='cTBS on 1st')
model_cond_sep_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                          data = use_dat,family = 'binomial')
model_cond_sep_2 <- glmer(OdorChosen ~ Run + Sess + (1|SubID), 
                          data = use_dat,family = 'binomial')
anova(model_cond_sep_1,model_cond_sep_2)
summary(model_cond_sep_2)
# sess only matters for 'cTBS on 1st' group of ppts



# plot overall learning trend of each aOFC and pOFC group

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

# Summarize within-subject means
subject_means <- conditioning_dat %>%
  group_by(SubID, Run, StimLoc) %>%
  summarise(Choice = mean(OdorChosen, na.rm = TRUE), .groups = "drop")

# Then plot with mean + SE across subjects
disc_overall = ggplot(subject_means, aes(x = Run, y = Choice)) +
  facet_wrap2(~StimLoc, scales = 'free_y', strip = strip) +
  geom_line(aes(group = SubID), show.legend = FALSE, alpha = 0.5, linewidth = 0.2) +  # individual subject lines
  stat_summary(aes(group = 1), fun = mean, geom = "line", 
               linewidth = 1.5, color = "black") +  # group-level mean
  stat_summary(aes(group = 1), fun.data = mean_se, geom = "errorbar", 
               width = 0.3, linewidth = 1, color = "black") +  # SE bars
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = "Runs", y = "P (Choosing rewarding stim)",
       title = NULL) + common


pdf(file.path(FigPaperDir,'Conditioning_overall.pdf'),8,4)
print(disc_overall)
dev.off()
