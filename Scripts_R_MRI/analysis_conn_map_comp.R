

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

tmp_path = '/Volumes/X9Pro/NODEAP/ConnMap_PAID_global/Agg'

subs_corr_cTBS = read.csv(file.path(tmp_path,'subs_corr_map_cTBS.csv'),header = F)
subs_corr_cTBS = as.numeric(unlist(subs_corr_cTBS))
subs_corr_cTBS_df = data.frame('corr'=subs_corr_cTBS,
                          'SubID'=rep(Subs,4),
                          'StimLoc'=rep(SubInfo$StimLoc,4),
                          'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                          'TMS'='cTBS')

subs_corr_sham = read.csv(file.path(tmp_path,'subs_corr_map_sham.csv'),header = F)
subs_corr_sham = as.numeric(unlist(subs_corr_sham))
subs_corr_sham_df = data.frame('corr'=subs_corr_sham,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='sham')
subs_corr_both_df = rbind(subs_corr_cTBS_df,subs_corr_sham_df)

p0=subs_corr_both_df %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'All subjects',y = 'Global Conn (raw)', x = '') +
  common

p1=subs_corr_both_df %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'aOFC subjects only',y = 'Global Conn (raw)', x = '') +
  common

x = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='sham' & StimLoc=='Anterior')
y = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='cTBS' & StimLoc=='Anterior')
wilcox.test(x$corr,y$corr,paired = T,alternative = 'g')

p2=subs_corr_both_df %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'pOFC subjects only',y = 'Global Conn (raw)', x = '') +
  common

pdf(file.path(FigDir,'MeanConn.pdf'),8,5)
print(p0)
print(p1)
print(p2)
dev.off()


##################

subs_corr_cTBS = read.csv(file.path(tmp_path,'abs_subs_corr_map_cTBS.csv'),header = F)
subs_corr_cTBS = as.numeric(unlist(subs_corr_cTBS))
subs_corr_cTBS_df = data.frame('corr'=subs_corr_cTBS,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='cTBS')

subs_corr_sham = read.csv(file.path(tmp_path,'abs_subs_corr_map_sham.csv'),header = F)
subs_corr_sham = as.numeric(unlist(subs_corr_sham))
subs_corr_sham_df = data.frame('corr'=subs_corr_sham,
                               'SubID'=rep(Subs,4),
                               'StimLoc'=rep(SubInfo$StimLoc,4),
                               'roi'=rep(c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'),each=48),
                               'TMS'='sham')
subs_corr_both_df = rbind(subs_corr_cTBS_df,subs_corr_sham_df)

p0=subs_corr_both_df %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'All subjects',y = 'Global Conn (abs)', x = '') +
  common

p1=subs_corr_both_df %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'aOFC subjects only',y = 'Global Conn (abs)', x = '') +
  common

x = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='sham' & StimLoc=='Anterior')
y = subset(subs_corr_both_df,roi=='aOFC-stim' & TMS=='cTBS' & StimLoc=='Anterior')
wilcox.test(x$corr,y$corr,paired = T,alternative = 'g')

p2=subs_corr_both_df %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x = TMS, y = corr, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'pOFC subjects only',y = 'Global Conn (abs)', x = '') +
  common

pdf(file.path(FigDir,'MeanConn_abs.pdf'),8,5)
print(p0)
print(p1)
print(p2)
dev.off()

