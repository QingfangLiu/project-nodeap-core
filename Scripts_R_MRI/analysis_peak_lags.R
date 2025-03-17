
# this aggregates peak_lags results across sessions
# see matlab code for details

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load from matlab
matfile = '../ProcessedData/peak_lags_subs_sessions.mat'
matdat = readMat(matfile)
peak_lags = matdat$Mean.peak.lags
sections = unlist(matdat$sections) # 4 sections
n_sections = length(sections)

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

# convert this array to df
dfpeak_lags = structure(peak_lags,
                       .Dim = c(n_rest_names, n_sections, length(Subs)), 
                       .Dimnames = structure(list(Sess = rest_names, 
                                                  section = sections,
                                                  SubID = as.vector(Subs)), 
                                             .Names = c("Sess", "section" ,"SubID")))
dfpeak_lags = adply(dfpeak_lags, c(1,2,3))
names(dfpeak_lags) = c(names(dfpeak_lags)[1:3],'peak_lags')
dfpeak_lags$StimLoc = NA
dfpeak_lags$StimOrder = NA

for(subj in Subs){
  dfpeak_lags$StimLoc[dfpeak_lags$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
  dfpeak_lags$StimOrder[dfpeak_lags$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfpeak_lags$StimOrder=mapvalues(dfpeak_lags$StimOrder,
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

dfpeak_lags$TMS_types = factor(unlist(TMS_types))
levels(dfpeak_lags$TMS_types) = c('cTBS','null','sham')

# add time variable
dfpeak_lags$time <- as.numeric(factor(dfpeak_lags$Sess,
                                      levels = unique(dfpeak_lags$Sess)))

# change of peak lags across sessions
p1=dfpeak_lags %>%
  ggplot(aes(x=section,y=peak_lags,fill=Sess)) +
  geom_boxplot(width = 0.8, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=Sess), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = '', x = '') +
  common

p2=dfpeak_lags %>%
  ggplot(aes(x=section,y=peak_lags,fill=TMS_types)) +
  geom_boxplot(width = 0.8, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = '', x = '') +
  common

p3=dfpeak_lags %>%
  subset(!TMS_types %in% 'null') %>%
  ggplot(aes(x=section,y=peak_lags,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc,nrow = 2) +
  labs(title = '',y = '', x = '') +
  common
  
# reduce the data, so not including Day 0 data
# also only keep one sham and one cTBS data point for each subj
summary_dfpeak_lags = dfpeak_lags %>%
  subset(!TMS_types %in% 'null') %>%
  group_by(SubID,section,TMS_types) %>%
  summarize(peak_lags=mean(peak_lags,na.rm=T))

p4=summary_dfpeak_lags %>%
  ggplot(aes(x=TMS_types,y=peak_lags)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~section,nrow = 1) +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "wilcox.test", 
                     alternative = "t",
                     paired = T,
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.format", 
                     hide.ns = F) + 
  labs(title = '',y = '', x = '') +
  common

filename = 'func_conn_peak_lags.pdf'
pdf(file.path(FigDir,filename),8,5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

