
rm(list = ls())
source('Setup.R')

# standarize ratings before stat testing & more comparisons
load('../../data_beh_processed/Odor_ratings_dat.RData')

p0=ggplot(Odor_ratings_dat,aes(x=IfDevalue,y=Pleasant)) +
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
  annotate("text", x = 1.5, y = 11.3, label = "***", size = 4) +
  labs(x = NULL, y = "Odor pleasantness", title = NULL) + common +
  theme(legend.position = c(0.25, 0.15)) 

pdf(file.path(FigPaperDir,'PleasantRatings_main.pdf'),4,4)
print(p0)
dev.off()


# 📦 Create customized facet strip styling
strip <- strip_themed(
  background_y = elem_list_rect(fill = use.col.ap.ofc),
  text_x = elem_list_text(color = use.col.conds, face = "bold", size = 16),
  text_y = elem_list_text(color = "white", face = "bold", size = 16)
)

p1=ggplot(Odor_ratings_dat,aes(x=IfDevalue,y=Pleasant)) +
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

p2=ggplot(Odor_ratings_dat,aes(x=IfDevalue,y=Pleasant)) +
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

pdf(file.path(FigPaperDir,'PleasantRatings_comp.pdf'),8,12)
ggarrange(p1,p2,nrow = 2)
dev.off()

pdf(file.path(FigPaperDir,'PleasantRatings_by_cond.pdf'),8,6)
print(p1)
dev.off()


# standardize those ratings (z-score) within each sub
dat <- Odor_ratings_dat %>%
  group_by(SubID) %>%
  mutate(Pleasant = (Pleasant - mean(Pleasant, na.rm = TRUE)) / 
           sd(Pleasant, na.rm = TRUE)) %>%
  ungroup()

# 🔄 Compute change in pleasantness (Post - Pre) for each odor per subject
# Create a new data frame `reduced` with only the change scores
pre_plea  <- subset(dat, PrePost == 'Pre-meal')$Pleasant
post_plea <- subset(dat, PrePost == 'Post-meal')$Pleasant
change_plea <- post_plea - pre_plea

reduced <- subset(dat, PrePost == 'Pre-meal') %>%
  mutate(
    PrePost = 'Change',
    Pleasant = change_plea
  )

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

# save this selective satiation index for later use
SelectSate_dat = rereduced
Sate_dat = reduced
save(SelectSate_dat,Sate_dat,file = '../ProcessedData/SelectSate_dat.RData')
