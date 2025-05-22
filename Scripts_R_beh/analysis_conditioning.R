
rm(list = ls())
source('Setup.R')

###############################
#### analyze conditioning data ####
###############################

conditioning_dat1 = read.xlsx('../ProcessedData/Conditioning.xlsx')
#conditioning_dat2 = read.xlsx('../ProcessedData/Conditioning_NODEAP_17.xlsx')
# not include NODEAP_17 due to lack of cTBS-sham session

conditioning_dat = conditioning_dat1 %>%
  mutate(SelectedCue = if_else(SelectedCue == 0, NA_integer_, SelectedCue)) %>% # replace zero with NA
  mutate(Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('Sweet','Savory')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                        to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                             'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('C-S-S','C-S-S','S-C-S','S-S-C','S-C-S','S-S-C')),
         StimOrder_day2 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('S-C-S','S-S-C','C-S-S','C-S-S','S-S-C','S-C-S'))) %>%
  mutate(SubID=factor(SubID),
         StimOrder=factor(StimOrder),
         StimOrder_day1=factor(StimOrder_day1),
         StimOrder_day2=factor(StimOrder_day2),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Set=factor(Set),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         ChoiceType=factor(ChoiceType),
         StimLoc=factor(StimLoc))

levels(conditioning_dat$StimLoc)=c('aOFC','pOFC')

###############################
save(conditioning_dat,file = '../ProcessedData/Conditioning.RData')
###############################

# does conditioning improve over runs?

p01=conditioning_dat %>%
  subset(StimLoc=='aOFC') %>%
  mutate(Sess=factor(Sess)) %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Sess,group=Sess)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.sess) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Anterior Stim)") 
p02=conditioning_dat %>%
  subset(StimLoc=='pOFC') %>%
  mutate(Sess=factor(Sess)) %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Sess,group=Sess)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.sess) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Posterior Stim)") 

p03=conditioning_dat %>%
  subset(StimLoc=='aOFC') %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Cond_day1,group=Cond_day1)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Anterior Stim)") 

p04=conditioning_dat %>%
  subset(StimLoc=='pOFC') %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=Cond_day1,group=Cond_day1)) +
  facet_wrap(~SubID) +
  geom_line() +
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "", y = "P(Choosing Odors)",
       title = "Conditioning (Posterior Stim)") 

pdf(file.path(FigDir,'Conditioning_subs.pdf'),10,8)
print(p01)
print(p02)
print(p03)
print(p04)
dev.off()

# collapsing across 3 conditions/sessions
p1=conditioning_dat %>%
  group_by(SubID,Run) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=SubID)) +
  geom_line(show.legend = F,alpha=0.5) +
  geom_line(stat = "summary", fun = "mean", linewidth = 2, color = "black") +
  labs(x = "Runs", y = "P(Choosing Odors)",
       title = "Conditioning (all subjects)") +
  common

model_conditioning_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_0 <- glmer(OdorChosen ~ (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_1,model_conditioning_0)
summary(model_conditioning_1)
# yes conditioning improves over runs

# does different TMS expt condition affect learning?
p3=conditioning_dat %>%
  group_by(SubID,Run,Cond) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common

# same as p3 but using error bar
p3a=conditioning_dat %>%
  group_by(SubID,Run,Cond) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +  
  stat_summary(fun.data = "mean_se", na.rm = T, show.legend = F,
               geom = "errorbar", width = 0.2) +  
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common

model_conditioning_2 <- glmer(OdorChosen ~ Run + Cond + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
model_conditioning_3 <- glmer(OdorChosen ~ Run * Cond + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_2,model_conditioning_1)
anova(model_conditioning_2,model_conditioning_3)
summary(model_conditioning_3)
# interaction effect between cTBS-sham/sham-cTBS conditions and runs


# is this modulated by different TMS locations?
p4=conditioning_dat %>%
  group_by(SubID,Run,Cond,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond, fill = Cond)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) + 
  labs(x = "Runs", y = "Choosing Odors", color = "") +
  common +
  theme(legend.position = c(0.85,0.3))

model_conditioning_4 <- glmer(OdorChosen ~ Run * Cond + StimLoc + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_4,model_conditioning_3) # need StimLoc? p=0.056
model_conditioning_5 <- glmer(OdorChosen ~ Run * Cond * StimLoc + (1|SubID), 
                              data = conditioning_dat,family = 'binomial')
anova(model_conditioning_4,model_conditioning_5)
anova(model_conditioning_3,model_conditioning_5)
summary(model_conditioning_5)


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

# conditioning performance across runs, by day1 TMS
pp1=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", linewidth = 0, show.legend = F,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P (Selecting odor-paired stim)", 
       color = "Day1 TMS") +
  common +
  theme(legend.position = c(0.75, 0.2))



# conditioning performance across runs, by diff day1 TMS & stim locs
pp2=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "P(Selecting odor-paired stim)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  common +
  theme(legend.position = c(0.85,0.3))
pp2

# same to pp2, but using error bar
conditioning_dat_summary = conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T))

conditioning_dat_summary %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1)) +
  facet_wrap(~StimLoc,scales = 'free_y') +
  scale_color_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "errorbar", linewidth = 1, width = 0.3) +  
  labs(x = "Runs", y = "P (Selecting odor-paired stim)", 
       color = "Day1 TMS") +
  common +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(legend.position = c(0.85,0.3),
        strip.background = element_blank(),
        panel.border = element_blank())


model_conditioning1_4 <- glmer(OdorChosen ~ Run * Cond_day1 + StimLoc + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_4,model_conditioning1_3) # need StimLoc? p=0.055

model_conditioning1_5 <- glmer(OdorChosen ~ Run * Cond_day1 + Run * StimLoc + (1|SubID), 
                               data = conditioning_dat,family = 'binomial')
anova(model_conditioning1_4,model_conditioning1_5)
anova(model_conditioning1_3,model_conditioning1_5)
summary(model_conditioning1_5)
# from the models, it seems to suggest an interaction b/t runs & conditions
# and an interaction b.t runs & stim locations

# aOFC only
use.dat = subset(conditioning_dat,StimLoc=='aOFC')
use.dat = use.dat[complete.cases(use.dat$OdorChosen), ] # remove NA values

model_conditioning_aOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_aOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_aOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_aOFC_3 <- glmer(OdorChosen ~ I(Run^2) + Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')  # quadratic term of run
model_conditioning_aOFC_4 <- glmer(OdorChosen ~ I(Run^2) + Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')  

anova(model_conditioning_aOFC_0,model_conditioning_aOFC_1) # day1 TMS main effect
anova(model_conditioning_aOFC_1,model_conditioning_aOFC_2) # no interaction
summary(model_conditioning_aOFC_1)
# describe: faster learning with runs under sham for aOFC stimulation group
anova(model_conditioning_aOFC_1,model_conditioning_aOFC_3) # add quadratic term helps
summary(model_conditioning_aOFC_3)
anova(model_conditioning_aOFC_4,model_conditioning_aOFC_3) 


conditioning_dat$Cond_day1_StimLoc <- interaction(conditioning_dat$Cond_day1, conditioning_dat$StimLoc)
conditioning_dat = conditioning_dat[complete.cases(conditioning_dat$OdorChosen), ] # remove NA values

summary_data <- conditioning_dat %>%
  group_by(Run, Cond_day1_StimLoc) %>%
  summarize(actual_choice = mean(OdorChosen),
            se_choice = sd(OdorChosen) / sqrt(n()))
ggplot(summary_data,aes(x = Run, color = Cond_day1_StimLoc)) +
  #scale_color_manual(values = use.col.conds.day1) +
  geom_line(aes(y = actual_choice, color = Cond_day1_StimLoc)) +
  geom_point(aes(y = actual_choice, color = Cond_day1_StimLoc), size = 3) +
  geom_errorbar(aes(ymin = actual_choice - se_choice, 
                    ymax = actual_choice + se_choice, color = Cond_day1_StimLoc), width = 0.2) +
  labs(x = "Runs", y = "Choosing stimuli predicting odors", 
       color = "Day1 TMS") +
  common +
  theme(legend.position = c(0.85,0.3))

model_conditioning_both_3 <- glmer(OdorChosen ~ I(Run^2) + Run + Cond_day1_StimLoc + (1|SubID), 
                                   data = conditioning_dat,family = 'binomial')
summary(model_conditioning_both_3)

# pOFC only
use.dat = subset(conditioning_dat,StimLoc=='pOFC')
use.dat = use.dat[complete.cases(use.dat$OdorChosen), ] # remove NA values

model_conditioning_pOFC_0 <- glmer(OdorChosen ~ Run + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_1 <- glmer(OdorChosen ~ Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_2 <- glmer(OdorChosen ~ Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_3 <- glmer(OdorChosen ~ I(Run^2) + Run + Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
model_conditioning_pOFC_4 <- glmer(OdorChosen ~ I(Run^2) + Run * Cond_day1 + (1|SubID), 
                                   data = use.dat,family = 'binomial')
anova(model_conditioning_pOFC_0,model_conditioning_pOFC_1) # day1 TMS main effect
anova(model_conditioning_pOFC_1,model_conditioning_pOFC_2) # interaction b/t TMS & runs
summary(model_conditioning_pOFC_2)
anova(model_conditioning_pOFC_1,model_conditioning_pOFC_3)
anova(model_conditioning_pOFC_4,model_conditioning_pOFC_3)
summary(model_conditioning_pOFC_4)

# look at Sess, Cond_day1, runs together
# note that subjects in each panel are not the same 
conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, linetype = factor(Sess),
             color = factor(Sess), fill = factor(Sess))) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  facet_wrap(~StimLoc + Cond_day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       title = "Conditioning by stim loc & Day1 TMS",
       color = "Session", fill = "Session", linetype = "Session") +
  common

# flip TMS cond & sess number
conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, fill = Cond_day1,
             color = Cond_day1)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  facet_wrap(~StimLoc + factor(Sess)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  
  stat_summary(fun.data = "mean_se", linewidth = 0,
               geom = "ribbon", alpha = 0.3) +  
  labs(x = "Runs", y = "P(Choosing Odors)", 
       title = "Conditioning by stim loc & Session #",
       color = "", fill = "") +
  common

# overall change of conditioning across sessions 1,2,3
q_sess = conditioning_dat %>%
  group_by(SubID,Run,Sess) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess),linetype=factor(Sess),
             fill=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", show.legend = F,
               linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "ribbon", alpha = 0.3, linewidth = 0) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  labs(x = "Runs", y = "P(Selecting odor-paired stim)",
       color = 'Session', title = '') +
  common

# overall change of conditioning across sessions 1,2,3, by stim loc
conditioning_dat %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess),fill=factor(Sess))) +
  facet_wrap(~StimLoc,scales = 'free_y') +
  stat_summary(fun = mean, geom = "line", show.legend = F,
               linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F, width = 0.5,
               geom = "errorbar", linewidth = 0.8) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(x = "Runs", y = "P(Selecting odor-paired stim)",color = 'Session') +
  coord_cartesian(ylim = c(0.5, 1)) +
  common

# this one, each panel has the same group of people
# separated by stimulation order
levels(conditioning_dat$StimOrder_day1)=
  c('cTBS on 1st','cTBS on 2nd','cTBS on 3rd')

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

conditioning_dat %>%
  group_by(SubID,Run,Sess,StimOrder_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,linetype=factor(Sess),
             color=factor(Sess),fill=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +
  facet_wrap(~StimOrder_day1) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  scale_linetype_manual(values=c("dotted","twodash", "solid"))+
  labs(x = "Runs", y = "P(Selecting Odor-paired stim)", fill = 'Session',
       color = 'Session', linetype = 'Session',
       title = '') +
  common

# do some testing here
use_dat = subset(conditioning_dat,StimOrder_day1=='cTBS on 1st')
model_cond_sep_1 <- glmer(OdorChosen ~ Run + (1|SubID), 
                               data = use_dat,family = 'binomial')
model_cond_sep_2 <- glmer(OdorChosen ~ Run + Sess + (1|SubID), 
                          data = use_dat,family = 'binomial')
anova(model_cond_sep_1,model_cond_sep_2)
summary(model_cond_sep_2)
# sess only matters for 'cTBS on 1st' group of ppts

pdf(file.path(FigDir,'Conditioning.pdf'),10,8)
print(p1)
print(p3)
print(p4)
print(pp1)
print(pp2)
dev.off()

pdf(file.path(FigDir,'Conditioning_error_bar.pdf'),7,4)
print(pp2a)
dev.off()


# plot overall learning trend of each aOFC and pOFC group

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

disc_overall = conditioning_dat %>%
  group_by(SubID,Run,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", linewidth = 1, 
               show.legend = F, width = 0.3,
               geom = "errorbar") +  
  facet_wrap2(~StimLoc,scales = 'free_y',strip = strip) +
  labs(x = "Runs", y = "P (Choosing rewarding stim)", 
       title = NULL,
       color = "Day1 TMS") +
  coord_cartesian(ylim = c(0.5, 1)) +
  common + theme(legend.position = c(0.85, 0.3)) 


pdf(file.path(FigPaperDir,'Conditioning_overall.pdf'),8,4)
print(disc_overall)
dev.off()

############ check RT 

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))
prt1 = conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc,Sess) %>%
  reframe(meanRT=mean(rt,na.rm = T)) %>%
  ggplot(aes(x = Run, y = meanRT, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap2(~StimLoc,strip = strip,
             scales = 'free_y') +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1,show.legend = T) + 
  stat_summary(fun.data = "mean_se", show.legend = F, width = 0.5,
               geom = "errorbar", linewidth = 0.8) +  
  labs(x = "Runs", y = "Response times (ms)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  coord_cartesian(ylim = c(1000,1650)) +
  common +
  theme(legend.position = c(0.85, 0.8)) 

prt2 = conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc,Sess) %>%
  reframe(meanRT=mean(rt,na.rm = T)) %>%
  ggplot(aes(x = Run, y = meanRT, color = Cond_day1, fill = Cond_day1)) +
  facet_wrap(~StimLoc + factor(Sess)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = "mean_se", 
               geom = "ribbon", alpha = 0.3, linewidth = 0) +  
  labs(x = "Runs", y = "RT (ms)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  common

# rt change with sess 
prt3=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc,Sess,StimOrder_day1) %>%
  mutate(Sess=factor(Sess,levels=c(1,2,3))) %>%
  reframe(meanRT=mean(rt,na.rm = T)) %>%
  ggplot(aes(x = Run, y = meanRT, color = Sess, fill = Sess)) +
  facet_wrap2(~StimLoc,strip = strip,
             scales = 'free_y') +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  stat_summary(fun = mean, geom = "line", linewidth = 1,show.legend = T) + 
  stat_summary(fun.data = "mean_se", show.legend = F, width = 0.5,
               geom = "errorbar", linewidth = 0.8) +  
  labs(x = "Runs", y = "Response times (ms)", 
       color = "Day1 TMS", fill = "Day1 TMS") +
  coord_cartesian(ylim = c(1000,1650)) +
  common +
  theme(legend.position = c(0.85, 0.8)) 


pdf(file.path(FigDir,'Conditioning_rt.pdf'),10,8)
print(prt1)
print(prt2)
dev.off()

# no TMS effect at all
model_cond_rt0 <- lmer(rt ~ (1|SubID), data = conditioning_dat)
model_cond_rt1 <- lmer(rt ~ Run + (1|SubID), data = conditioning_dat)
summary(model_cond_rt1)
anova(model_cond_rt1,model_cond_rt0) # rt decreases over runs

model_cond_rt2 <- lmer(rt ~ Run + Cond_day1 + (1|SubID), data = conditioning_dat)
model_cond_rt3 <- lmer(rt ~ Run * Cond_day1 + (1|SubID), data = conditioning_dat)
anova(model_cond_rt1,model_cond_rt2)
anova(model_cond_rt3,model_cond_rt1)

model_cond_rt4 <- lmer(rt ~ Run + Sess + (1|SubID), data = conditioning_dat)
anova(model_cond_rt1,model_cond_rt4) # adding sess improves model fit
summary(model_cond_rt4)

model_cond_rt5 <- lmer(rt ~ Run + Sess + Cond_day1 + (1|SubID), data = conditioning_dat)
anova(model_cond_rt5,model_cond_rt4)

model_cond_rt6 <- lmer(rt ~ Run + Sess * Cond_day1 + (1|SubID), data = conditioning_dat)
anova(model_cond_rt6,model_cond_rt4)

##########################################

pdf(file.path(FigDir,'Conditioning_paper.pdf'),9,8)
ggarrange(prt1,prt3,
          ncol = 1,nrow = 2,
          common.legend = F)
dev.off()




