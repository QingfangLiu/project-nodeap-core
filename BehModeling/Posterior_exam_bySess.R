
rm(list = ls())
source('../Scripts_R/Setup.R')

opt = 'bySess'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'.rds'))
dic_value <- samples$BUGSoutput$DIC
pd_value <- samples$BUGSoutput$pD
cat("DIC:", dic_value, "\n")
cat("pD (Effective Number of Parameters):", pd_value, "\n")

#################################

load(file = '../ProcessedData/Conditioning.RData') 
# this contains some NA trials
use_dat = conditioning_dat  
sub = use_dat$SubID
nsubs = length(unique(sub))

# posterior values
alpha <- samples$BUGSoutput$sims.list$alpha
mu <- samples$BUGSoutput$sims.list$mu
kappa <- samples$BUGSoutput$sims.list$kappa

apply(mu,2,mean)
mean(kappa)

alpha_map = apply(alpha,c(2,3),mean)
write.csv(alpha_map,file = 'alpha_per_sub_sess.csv')
# save alpha MAP estimates

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
SubInfo = subset(SubInfo,SubID %in% unique(use_dat$SubID)) %>%
  mutate(StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS'))) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,
                                    from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                           'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('C-S-S','C-S-S','S-C-S','S-S-C','S-C-S','S-S-C')))

TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder_day1[j]
  if(tmp_order=='S-S-C') TMS_types[[j]] = c('sham','sham','cTBS')
  if(tmp_order=='S-C-S') TMS_types[[j]] = c('sham','cTBS','sham')
  if(tmp_order=='C-S-S') TMS_types[[j]] = c('cTBS','sham','sham')
}


df_alpha = data.frame(TMS_types = unlist(TMS_types),
                      Sess = rep(1:3,nsubs),
                      alpha = as.numeric(t(alpha_map)),
                      StimLoc = rep(SubInfo$StimLoc,each=3),
                      StimOrder_day1 = rep(SubInfo$StimOrder_day1,each=3),
                      SubID = rep(SubInfo$SubID,each=3)) %>%
  mutate(StimOrder_day1=factor(StimOrder_day1))


levels(df_alpha$StimOrder_day1)=
  c('cTBS on 1st','cTBS on 2nd','cTBS on 3rd')
ggplot(df_alpha,aes(x=TMS_types,y=alpha)) +
  facet_wrap(~StimOrder_day1) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1.5, alpha = 0.8) +
  labs(y='Learning rate',x='') +
  common

# this is consistent with the figure using just the choice
# the difference seems to be driven from the 1st session

# model testing: does adding TMS better explain the learning rates?
model0 <- lmer(alpha ~ Sess + (1|SubID), data = df_alpha)
summary(model0)
model1 <- lmer(alpha ~ Sess + TMS_types + (1|SubID), data = df_alpha)
anova(model1,model0) 
summary(model1)
model2 <- lmer(alpha ~ Sess * TMS_types + (1|SubID), data = df_alpha)
anova(model0,model2)

survey_dat = read.xlsx('../ProcessedData/Survey_uncomf_strong.xlsx') %>%
    subset(SubID %in% unique(use_dat$SubID)) %>%
    subset(TMSSess %in% c(1,3,5)) %>%
    mutate(rating = (uncom + strong)/2)
  
# attach some cols from df_alpha to survey_dat
survey_dat = cbind(survey_dat,df_alpha$alpha,df_alpha$StimOrder_day1)
names(survey_dat)[c(10,11)] = c('alpha','StimOrder_day1')

p1=ggplot(survey_dat,aes(x=rating,y=alpha,color=TMS_type)) +
  facet_wrap(~StimOrder_day1,scales = 'free_y') +
  geom_point(size=2.5,alpha=0.8) +
  #geom_smooth(method = "lm", se = FALSE) +  
  stat_cor(method = "pearson") +  
  labs(y='Learning rate',x='Perceived TMS discomfort & intensity') +
  scale_color_manual(values=use.col.conds.day1) +
  common +
  coord_cartesian(ylim = c(0,0.75),xlim = c(-2,2)) +
  theme(legend.position = '')

# avg two sham sessions for line plotting
df_alpha_summary = df_alpha %>%
  dplyr::group_by(TMS_types,SubID,StimLoc,StimOrder_day1) %>%
  reframe(alpha=mean(alpha))
p2=ggplot(df_alpha_summary,aes(x=TMS_types,y=alpha,fill = TMS_types)) +
  facet_wrap(~StimOrder_day1,scales = 'free') +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID,color=TMS_types), position = pd, 
              size = 2, alpha = 0.8) +
  labs(y='Learning rate',x=NULL) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_compare_means(method = "wilcox.test",label.y = 0.8,
                     paired = T) +
  common +
  coord_cartesian(ylim = c(-0.1,0.8)) +
  theme(legend.position = '')

pdf('learning_rate_uncomf_strong.pdf',9,8)
ggarrange(p2,p1,nrow = 2,heights = c(1,1))
dev.off()


