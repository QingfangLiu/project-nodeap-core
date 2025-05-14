
rm(list = ls())
source('Setup.R')

# This code looks at devaluation effect
# 9.27.24: major changes:
# add the abnormal NODEAP_17 into this dataset
# standarize ratings before stat testing & more comparisons

###############################
# analyze pleasantness ratings data
###############################

dat1 = read.xlsx('../ProcessedData/Pleasant.xlsx')
dat2 = read.xlsx('../ProcessedData/Pleasant_NODEAP_17.xlsx')

dat = rbind(dat1,dat2) %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre-meal','Post-meal')),
         IfDevalue=mapvalues(IfDevalue,from=c('0','1'),to=c('non-sated','sated')),
         StimLoc=mapvalues(StimLoc,from=c('Anterior','Posterior'),to=c('aOFC','pOFC')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre-meal','Post-meal')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         IfDevalue=factor(IfDevalue,levels=c('non-sated','sated'))) %>%
  arrange(SubID,Cond,IfDevalue)

dat = dat %>% 
  mutate(pair=rep(c(1:(nrow(dat)/2)),each=2))

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = use.col.conds,
                                             face = "bold",
                                             size = 16),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

p1=ggplot(dat,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = F) +
  geom_jitter(aes(shape=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.6) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(Cond), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  theme(legend.position = 'inside',
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.35,'cm'),
        legend.position.inside = c(0.1,0.1))

strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.ofc),
                     text_x = elem_list_text(color = use.col.sess,
                                             face = "bold",
                                             size = 16),
                     text_y = elem_list_text(color = 'white',
                                             face = "bold",
                                             size = 16))

p2=ggplot(dat,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(aes(linetype = PrePost),
               width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0, show.legend = F) +
  geom_jitter(aes(shape=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.6) +
  facet_grid2(rows = vars(StimLoc), 
              cols = vars(Sess), 
              scales = 'fixed',
              strip = strip,
              axes = 'all') +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  theme(legend.position = 'inside',
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.35,'cm'),
        legend.position.inside = c(0.1,0.1))

p3=ggplot(dat,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(aes(linetype = PrePost),
               outlier.alpha = 0, alpha=0.4, show.legend = F) +
  geom_jitter(aes(shape=PrePost),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.2), 
              size = 1.5, alpha = 0.4) +
  scale_shape_manual(values = c(17, 1)) +
  annotate("segment", x = 1, xend = 1, y = 9.8, yend = 10.5) +
  annotate("segment", x = 2, xend = 2, y = 9.8, yend = 10.5) +
  annotate("segment", x = 1, xend = 2, y = 10.5, yend = 10.5) +
  annotate("text", x = 1.5, y = 11.3, label = "p = 2.75e-13", size = 4) +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  theme(legend.position = c(0.25, 0.15)) 

pdf(file.path(FigPaperDir,'PleasantRatings_comp.pdf'),8,12)
ggarrange(p1,p2,nrow = 2)
dev.off()

pdf(file.path(FigPaperDir,'PleasantRatings_main.pdf'),4,4)
print(p3)
dev.off()




# standardize those ratings (z-score) within each sub
dat <- dat %>%
  group_by(SubID) %>%
  mutate(Pleasant = (Pleasant - mean(Pleasant, na.rm = TRUE)) / 
           sd(Pleasant, na.rm = TRUE)) %>%
  ungroup()

# to plot the change of pleasantness ratings (post-pre meal)
# create a new df 'reduced' from 'dat'
pre_plea = subset(dat,PrePost=='Pre-meal')$Pleasant
post_plea = subset(dat,PrePost=='Post-meal')$Pleasant
change_plea = post_plea - pre_plea
reduced = subset(dat,PrePost=='Pre-meal')
reduced$PrePost = 'Change'
reduced$Pleasant = change_plea

reduced = reduced %>% 
  mutate(pair=rep(c(1:(nrow(reduced)/2)),each=2))


# Fit the mixed-effects model
model0 <- lmer(Pleasant ~ (1|SubID), data = reduced)
model1 <- lmer(Pleasant ~ IfDevalue + (1|SubID), data = reduced)
anova(model1,model0) # collapsing across all 3 sessions, there was a devaluation effect
summary(model1)

model2 <- lmer(Pleasant ~ Cond + IfDevalue + (1|SubID), data = reduced)
model3 <- lmer(Pleasant ~ Cond * IfDevalue + (1|SubID), data = reduced)

anova(model2,model1)
summary(model2)
summary(model3)
anova(model2,model3)
anova(model3,model1)
# the devaluation effect was not modulated by conditions


# further create a measure of selective satiation effect
# per session per subject
# change of pleasantness rating (sated - non-sated)
# post-pre: more negative value means better selective satiation effect

ND = subset(reduced,IfDevalue=='non-sated')$Pleasant
D = subset(reduced,IfDevalue=='sated')$Pleasant
Didx = D - ND
rereduced = subset(reduced,IfDevalue=='non-sated')
# these cols are meanless now after reducing
rereduced$IfDevalue = NULL
rereduced$PrePost = NULL 
rereduced$Pleasant = NULL
rereduced$IfDevalue = NULL
rereduced$OdorType = NULL
rereduced$pair = NULL
rereduced$Didx = Didx

# Compare the mixed-effects model
model0 <- lmer(Didx ~ (1|SubID), data = rereduced)
model1 <- lmer(Didx ~ Cond + (1|SubID), data = rereduced)
model2 <- lmer(Didx ~ Sess + (1|SubID), data = rereduced)
model3 <- lmer(Didx ~ Devalued + (1|SubID), data = rereduced)
model4 <- lmer(Didx ~ StimLoc + (1|SubID), data = rereduced)
summary(model0)
anova(model0,model1)
anova(model0,model2)
anova(model0,model3)
summary(model3)
anova(model0,model4)
# devaluation effect was not modulated by Cond, Sess, DevaluedOdor, StimLoc


# save the whole rating data for later use
Odor_ratings_dat = dat
save(Odor_ratings_dat,file = '../ProcessedData/Odor_ratings_dat.RData')

# save this selective satiation index for later use
SelectSate_dat = rereduced
Sate_dat = reduced
save(SelectSate_dat,Sate_dat,
     file = '../ProcessedData/SelectSate_dat.RData')


################## below is not included in the paper #######################



p1=ggplot(rereduced,aes(x=Cond,y=Didx,color=Cond,fill=Cond)) +
  geom_jitter(aes(shape=StimLoc),size = 2, alpha = 0.8,show.legend = F) +
  geom_violin(alpha = 0.1) +
  geom_boxplot(alpha = 0, width = .1, show.legend = F) +
  facet_wrap(~StimLoc) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'black') +
  labs(x = "", y = "Selective satiation index (SSI)") + common +
  scale_color_manual(values = use.col.conds) +
  scale_fill_manual(values = use.col.conds) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 0.5))

p2=ggplot(rereduced,aes(x=factor(Sess),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Selective satiation effect") + common +
  theme(strip.background = element_blank())

p3=ggplot(rereduced,aes(x=factor(Devalued),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Selective satiation effect") + common +
  theme(strip.background = element_blank())

p4=ggplot(rereduced,aes(x=factor(StimLoc),y=Didx)) +
  geom_jitter(size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, outlier.alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  labs(x = "", y = "Selective satiation effect") + common +
  theme(strip.background = element_blank())

p5=ggplot(rereduced,aes(x=Cond_day2,y=Didx,color=Cond_day2,fill=Cond_day2)) +
  geom_jitter(aes(shape=StimLoc),size = 2, alpha = 0.8,show.legend = F) +
  geom_violin(alpha = 0.1) +
  geom_boxplot(alpha = 0, width = .1, show.legend = F) +
  facet_wrap(~StimLoc) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'black') +
  labs(x = "", y = "Selective satiation index (SSI)") + common +
  scale_color_manual(values = use.col.conds.day2) +
  scale_fill_manual(values = use.col.conds.day2) +
  theme(strip.background = element_blank())

pdf(file.path(FigDir,'Selective_satiation_effect.pdf'),4.5,4)
print(p2)
print(p3)
print(p4)
dev.off()

pdf(file.path(FigDir,'Selective_satiation_effect_by_cond.pdf'),9,4)
print(p1)
dev.off()

pdf(file.path(FigDir,'Selective_satiation_effect_by_cond_day2.pdf'),7,4)
print(p5)
dev.off()


# a diff way to plot to show lines pre & post meal
p=ggplot(dat,aes(x=PrePost,y=Pleasant,fill=PrePost)) +
  geom_boxplot(width = 0.6, alpha = 0.4, outlier.alpha = 0) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_jitter(aes(color=PrePost,group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  facet_wrap(~IfDevalue+Cond) +
  labs(x = "", y = "Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        legend.position = 'none')

pdf(file.path(FigDir,'PleasantRatings_w_lines.pdf'),8,8)
print(p)
dev.off()


q1=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~Cond) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q1a=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~StimLoc + Cond) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q2=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.5, color = 'darkgray') +
  facet_wrap(~Sess) +
  labs(x = "", y = "Change of Pleasantness Ratings") + common +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

q3=ggplot(reduced,aes(x=IfDevalue,y=Pleasant)) +
  geom_line(aes(group=pair), position = pd, linewidth = 0.2, color = 'darkgray') +
  geom_jitter(aes(group=pair), position = pd, 
              size = 1, alpha = 0.8) +
  geom_boxplot(width = 0.6, alpha = 0, position = position_dodge(0.75),
               outlier.alpha = 0) +
  labs(x = "", y = "Pleasantness Ratings (Post-Pre)") + common +
  theme(strip.background = element_blank())

pdf(file.path(FigDir,'PleasantChanges.pdf'),6,4)
print(q1)
print(q2)
print(q1a)
dev.off()

# showing standardized values
pdf(file.path(FigDir,'PleasantChanges_all_sessions.pdf'),4,5)
print(q3)
dev.off()




