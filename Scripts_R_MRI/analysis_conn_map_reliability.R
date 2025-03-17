

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

tmp_path = '/Volumes/X9Pro/NODEAP/ConnMap_PAID_global/Agg'
subs_corr = read.csv(file.path(tmp_path,'subs_correlation.csv'),header = F)
subs_corr = as.numeric(unlist(subs_corr))

subs_corr_df = data.frame('corr'=subs_corr,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48))
q0=subs_corr_df %>%
  ggplot(aes(x = roi, y = corr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  labs(title = 'All subjects',y = 'Conn Pattern Reliability', x = '') +
  common

subs_corr_cTBS = read.csv(file.path(tmp_path,'subs_correlation_cTBS.csv'),header = F)
subs_corr_cTBS = as.numeric(unlist(subs_corr_cTBS))
subs_corr_cTBS_df = data.frame('corr'=subs_corr_cTBS,
                          'SubID'=rep(Subs,4),
                          'StimLoc'=rep(SubInfo$StimLoc,4),
                          'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                          'TMS'='cTBS')

subs_corr_sham = read.csv(file.path(tmp_path,'subs_correlation_sham.csv'),header = F)
subs_corr_sham = as.numeric(unlist(subs_corr_sham))
subs_corr_sham_df = data.frame('corr'=subs_corr_sham,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='sham')
subs_corr_both_df = rbind(subs_corr_cTBS_df,subs_corr_sham_df)

subs_corr_diff = subs_corr_sham - subs_corr_cTBS
subs_corr_diff_df = data.frame('corr_diff'=subs_corr_diff,
                               'corr'=subs_corr,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48))
q1=subs_corr_diff_df %>%
  ggplot(aes(x=corr,y=corr_diff)) +
  geom_point() +
  geom_smooth(method = "lm", col = "purple") +
  stat_cor(method = "pearson") +
  facet_wrap(~roi) + 
  labs(x = 'Overall Reliability', y = 'Reliability (Sham - cTBS)') +
  common

q2=subs_corr_diff_df %>%
  ggplot(aes(x=corr,y=corr_diff,group = StimLoc)) +
  geom_point(aes(shape=StimLoc)) +
  geom_smooth(method = "lm", col = "purple") +
  stat_cor(method = "pearson") +
  facet_wrap(~roi) + 
  labs(x = 'Overall Reliability', y = 'Reliability (Sham - cTBS)') +
  common


p0=subs_corr_both_df %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'All subjects',y = 'Conn Pattern Reliability', x = '') +
  common

x = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='sham')
y = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='cTBS')
wilcox.test(x$corr,y$corr,paired = T,alternative = 'g') # p=0.1552

# Function to perform Fisher z-transformation
fisher_z <- function(r) {
  return(0.5 * log((1 + r) / (1 - r)))
}

# test with overall corr as a covariate
# this model doesn't work well: covergence issue
subs_corr_both_df$overall_corr = rep(subs_corr,2)
use.dat = subs_corr_both_df %>%
  subset(roi %in% c('aOFC-stim')) %>%
  mutate(SubID=factor(SubID),
         fisher_z_corr = fisher_z(corr),
         TMS=factor(TMS))
model_1 <- lmer(fisher_z_corr ~ TMS + overall_corr + (1|SubID), data = use.dat)
model_0 <- lmer(fisher_z_corr ~ overall_corr + (1|SubID), data = use.dat)
anova(model_1,model_0)
summary(model_0)


p1=subs_corr_both_df %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'aOFC subjects only',y = 'Conn Pattern Reliability', x = '') +
  common

p2=subs_corr_both_df %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'pOFC subjects only',y = 'Conn Pattern Reliability', x = '') +
  common

pdf(file.path(FigDir,'ConnPatternReliability.pdf'),8,5)
print(p0)
print(p1)
print(p2)
dev.off()

pdf(file.path(FigDir,'ConnPatternReliability_overall_diff.pdf'),7,6)
print(q0)
print(q1)
print(q2)
dev.off()


ConnIdx = read.csv('../ProcessedData/ConnIdx.csv')
plot(ConnIdx$ACI,tmp$change)
cor.test(ConnIdx$ACI,tmp$change)

df = ConnIdx[,c('ACI','PCI','OCI','RCI','ARCI')]
df$aOFC_seed=subset(subs_corr_diff_df,roi=='aOFC-seed')$corr_diff
df$pOFC_seed=subset(subs_corr_diff_df,roi=='pOFC-seed')$corr_diff
df$aOFC_stim=subset(subs_corr_diff_df,roi=='aOFC-stim')$corr_diff
df$pOFC_stim=subset(subs_corr_diff_df,roi=='pOFC-stim')$corr_diff

library(GGally)
ggpairs(df)

########################################################
## bootstrap approach to do it

tmp_path = '/Volumes/X9Pro/NODEAP/ConnMap_PAID_global/Agg'
subs_corr_cTBS = read.csv(file.path(tmp_path,'subs_correlation_cTBS_bs.csv'),header = F)
subs_corr_cTBS = as.numeric(unlist(subs_corr_cTBS))
subs_corr_cTBS_df = data.frame('corr'=subs_corr_cTBS,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='cTBS')

subs_corr_sham = read.csv(file.path(tmp_path,'subs_correlation_sham_bs.csv'),header = F)
subs_corr_sham = as.numeric(unlist(subs_corr_sham))
subs_corr_sham_df = data.frame('corr'=subs_corr_sham,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='sham')
subs_corr_both_df = rbind(subs_corr_cTBS_df,subs_corr_sham_df)

p0_bs=subs_corr_both_df %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'All subjects',y = 'Conn Pattern Reliability', x = '') +
  common

x = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='sham')
y = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='cTBS')
wilcox.test(x$corr,y$corr,paired = T,alternative = 'g')

p1_bs=subs_corr_both_df %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'aOFC subjects only',y = 'Conn Pattern Reliability', x = '') +
  common

p2_bs=subs_corr_both_df %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'pOFC subjects only',y = 'Conn Pattern Reliability', x = '') +
  common

x = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='sham' & StimLoc=='Anterior')
y = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='cTBS' & StimLoc=='Anterior')
wilcox.test(x$corr,y$corr,paired = T,alternative = 'g')

pdf(file.path(FigDir,'ConnPatternReliability_bs.pdf'),8,5)
print(p0_bs)
print(p1_bs)
print(p2_bs)
dev.off()


