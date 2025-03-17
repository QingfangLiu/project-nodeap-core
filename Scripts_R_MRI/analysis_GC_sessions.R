
# this aggregates roi-based global conn results 

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load conn vals from matlab
Option = 'Signed' # Unsigned or Signed or Squared
matfile = paste0('../ProcessedData/GC_subs_sessions_',Option,'.mat')
matdat = readMat(matfile)
convals = matdat$MeanGC

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

roinames = c('aOFC_seed','aOFC_stim','pOFC_seed','pOFC_stim');
n_roi = length(roinames)

# convert this array to df
dfconvals = structure(convals,
                       .Dim = c(n_rest_names, n_roi, length(Subs)), 
                       .Dimnames = structure(list(Sess = rest_names, 
                                                  roi = roinames,
                                                  SubID = as.vector(Subs)), 
                                             .Names = c("Sess", "roi" ,"SubID")))
dfconvals = adply(dfconvals, c(1,2,3))
names(dfconvals) = c(names(dfconvals)[1:3],'GC')
dfconvals$StimLoc = NA
dfconvals$StimOrder = NA

for(subj in Subs){
dfconvals$StimLoc[dfconvals$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
dfconvals$StimOrder[dfconvals$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfconvals$StimOrder=mapvalues(dfconvals$StimOrder,
                    from=c('123','132','213','231','312','321'),
                    to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                         'SC-SS-CS','SS-CS-SC','SS-SC-CS'))

# add a column called TMS_types indicating the TMS type of each session

TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order==123) TMS_types[[j]] = rep(c('N','C','S','S','C','S','S'),4)
  if(tmp_order==132) TMS_types[[j]] = rep(c('N','C','S','S','S','S','C'),4)
  if(tmp_order==213) TMS_types[[j]] = rep(c('N','S','C','C','S','S','S'),4)
  if(tmp_order==231) TMS_types[[j]] = rep(c('N','S','C','S','S','C','S'),4)
  if(tmp_order==312) TMS_types[[j]] = rep(c('N','S','S','C','S','S','C'),4)
  if(tmp_order==321) TMS_types[[j]] = rep(c('N','S','S','S','C','C','S'),4)
}

dfconvals$TMS_types = factor(unlist(TMS_types))
levels(dfconvals$TMS_types) = c('cTBS','null','sham')


Sess_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order==123) Sess_types[[j]] = rep(c('N','cTBS-sham','cTBS-sham','sham-cTBS','sham-cTBS','sham-sham','sham-sham'),4)
  if(tmp_order==132) Sess_types[[j]] = rep(c('N','cTBS-sham','cTBS-sham','sham-sham','sham-sham','sham-cTBS','sham-cTBS'),4)
  if(tmp_order==213) Sess_types[[j]] = rep(c('N','sham-cTBS','sham-cTBS','cTBS-sham','cTBS-sham','sham-sham','sham-sham'),4)
  if(tmp_order==231) Sess_types[[j]] = rep(c('N','sham-cTBS','sham-cTBS','sham-sham','sham-sham','cTBS-sham','cTBS-sham'),4)
  if(tmp_order==312) Sess_types[[j]] = rep(c('N','sham-sham','sham-sham','cTBS-sham','cTBS-sham','sham-cTBS','sham-cTBS'),4)
  if(tmp_order==321) Sess_types[[j]] = rep(c('N','sham-sham','sham-sham','sham-cTBS','sham-cTBS','cTBS-sham','cTBS-sham'),4)
}

dfconvals$Sess_types = factor(unlist(Sess_types))

## summary for reporting to ClinicalTrial purpose
# used cTBS session for sham-cTBS, or cTBS-sham
# but the average for sham-sham
summary_dfconvals = dfconvals %>%
  subset(roi %in% c('aOFC_seed','pOFC_seed')) %>%
  subset(!Sess %in% 'D0') %>%
  subset((Sess_types=='cTBS-sham' & TMS_types=='C') | 
           (Sess_types=='sham-cTBS' & TMS_types=='C') |
           (Sess_types=='sham-sham' & TMS_types=='S')
           ) %>%
  group_by(SubID,Sess_types,StimLoc) %>%
  summarize(meanGC=mean(GC,na.rm=T)) %>%
  ungroup()

summary_summary_dfconvals = summary_dfconvals %>%
  group_by(Sess_types,StimLoc) %>%
  summarize(mean=mean(meanGC,na.rm=T),
            sd=sd(meanGC))

p0=dfconvals %>%
  ggplot(aes(x=roi,y=GC,fill=Sess)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=Sess), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = 'Global conn', x = '') +
  common

p1=dfconvals %>%
  ggplot(aes(x=roi,y=GC,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc,nrow = 2) +
  labs(title = 'All sessions',y = 'Global conn', x = '') +
  common
  
p1a=dfconvals %>%
  ggplot(aes(x=roi,y=GC,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = 'All sessions',y = 'Global conn', x = '') +
  common

# if only using Day 1 TMS data
p2=dfconvals %>%
  subset(Sess %in% c('S1D1','S2D1','S3D1')) %>%
  ggplot(aes(x=roi,y=GC,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc,nrow = 2) +
  labs(title = 'Day 1 sessions',y = 'Global conn', x = '') +
  common

p2a=dfconvals %>%
  subset(Sess %in% c('S1D1','S2D1','S3D1')) %>%
  ggplot(aes(x=roi,y=GC,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = 'Day 1 sessions',y = 'Global conn', x = '') +
  common

# if only using Day 2 TMS data
p3=dfconvals %>%
  subset(Sess %in% c('S1D2','S2D2','S3D2')) %>%
  ggplot(aes(x=roi,y=GC,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc,nrow = 2) +
  labs(title = 'Day 2 sessions',y = 'Global conn', x = '') +
  common

# average across S1,S2,S3 for Day 1
summary_dfconvals_D1 = dfconvals %>%
  subset(Sess %in% c('S1D1','S2D1','S3D1')) %>%
  group_by(SubID,roi,TMS_types,StimLoc) %>%
  summarize(meanGC=mean(GC,na.rm=T))

s1=summary_dfconvals_D1 %>%
ggplot(aes(x=TMS_types,y=meanGC)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,nrow = 1) +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "wilcox.test", 
                     alternative = "t",
                     paired = T,
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.format", 
                     hide.ns = F) + 
  labs(title = 'Day 1 sessions',y = 'Global conn', x = '') +
  common

summary_dfconvals_D2 = dfconvals %>%
  subset(Sess %in% c('S1D2','S2D2','S3D2')) %>%
  group_by(SubID,roi,TMS_types,StimLoc) %>%
  summarize(meanGC=mean(GC,na.rm=T))

s2=summary_dfconvals_D2 %>%
  ggplot(aes(x=TMS_types,y=meanGC)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,nrow = 1) +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "wilcox.test", 
                     alternative = "t",
                     paired = T,
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.format", 
                     hide.ns = F) + 
  labs(title = 'Day 2 sessions',y = 'Global conn', x = '') +
  common


# do some correlations
# if overall there is no sig diff b/t sham and cTBS
# would this depend on the Conn indices?

cTBS_vals_df = summary_dfconvals_D2 %>% 
  subset(TMS_types=='cTBS' & roi=='pOFC_seed')
cTBS_vals = cTBS_vals_df$meanGC
sham_vals_df = summary_dfconvals %>% 
  subset(TMS_types=='sham'& roi=='pOFC_seed')
sham_vals = sham_vals_df$meanGC
TMS_effect_df = subset(sham_vals_df, select = -c(TMS_types,meanGC))
TMS_effect_df$change = sham_vals - cTBS_vals
tmp = subset(TMS_effect_df,roi=='aOFC_seed' & StimLoc=='Anterior')
hist(tmp$change)
median(tmp$change)

ConnIdx = read.csv('../ProcessedData/ConnIdx.csv')
ConnIdx = subset(ConnIdx,StimLoc=='Anterior')
plot(ConnIdx$ACI,tmp$change)
cor.test(ConnIdx$ACI,tmp$change)
plot(ConnIdx$PCI,tmp$change)
cor.test(ConnIdx$PCI,tmp$change)
plot(ConnIdx$OCI,tmp$change)
cor.test(ConnIdx$OCI,tmp$change)
plot(ConnIdx$RCI,tmp$change)
cor.test(ConnIdx$RCI,tmp$change)

filename = paste0('Global_conn_sessions_',Option,'.pdf')
pdf(file.path(FigDir,filename),8,7)
print(p0)
print(p1)
print(p2)
print(p3)
print(s1)
print(s2)
dev.off()


