
# this aggregates ALFF results across sessions
# options: ALFF, fALFF, or zALFF
# see matlab code for details

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load from matlab
use_measure = 'ALFF'
matfile = paste0('../ProcessedData/ALFF_subs_sessions_',use_measure,'.mat')
matdat = readMat(matfile)
ALFF = matdat$MeanALFF

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

# this needs to match the matlab code
roinames = c('aOFC_seed','aOFC_stim','pOFC_seed','pOFC_stim');
n_roi = length(roinames)

# convert this array to df
dfALFF = structure(ALFF,
                       .Dim = c(n_rest_names, n_roi, length(Subs)), 
                       .Dimnames = structure(list(Sess = rest_names, 
                                                  roi = roinames,
                                                  SubID = as.vector(Subs)), 
                                             .Names = c("Sess", "roi" ,"SubID")))
dfALFF = adply(dfALFF, c(1,2,3))
names(dfALFF) = c(names(dfALFF)[1:3],'ALFF')
dfALFF$StimLoc = NA
dfALFF$StimOrder = NA

for(subj in Subs){
  dfALFF$StimLoc[dfALFF$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
  dfALFF$StimOrder[dfALFF$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfALFF$StimOrder=mapvalues(dfALFF$StimOrder,
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

dfALFF$TMS_types = factor(unlist(TMS_types))
levels(dfALFF$TMS_types) = c('cTBS','null','sham')

# add time variable
dfALFF$time <- as.numeric(factor(dfALFF$Sess,levels = unique(dfALFF$Sess)))

# exclude data outside of 3sd, for each roi
# only do this for ALFF value (other two are scaled ok)
if (use_measure=='ALFF'){
  print(paste('# data points before excluding is',nrow(dfALFF)))
  dfALFF = dfALFF %>%
   group_by(roi) %>%
   filter(abs(ALFF - mean(ALFF,na.rm=T)) <= 3 * sd(ALFF,na.rm = T))
  print(paste('# data points after excluding is',nrow(dfALFF)))
}

# change of ALFF across sessions
p1=dfALFF %>%
  ggplot(aes(x=roi,y=ALFF,fill=Sess)) +
  geom_boxplot(width = 0.8, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=Sess), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = use_measure, x = '') +
  common

p2=dfALFF %>%
  subset(!TMS_types %in% 'null') %>%
  ggplot(aes(x=roi,y=ALFF,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = use_measure, x = '') +
  common

p3=dfALFF %>%
  subset(!TMS_types %in% 'null') %>%
  ggplot(aes(x=roi,y=ALFF,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc,nrow = 2) +
  labs(title = '',y = use_measure, x = '') +
  common
  
# reduce the data, so not including Day 0 data
# also only keep one sham and one cTBS data point for each subj
summary_dfALFF = dfALFF %>%
  subset(!TMS_types %in% 'null') %>%
  group_by(SubID,roi,TMS_types) %>%
  summarize(ALFF=mean(ALFF,na.rm=T))

p4=summary_dfALFF %>%
  ggplot(aes(x=TMS_types,y=ALFF)) +
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
  labs(title = '',y = use_measure, x = '') +
  common

filename = paste0(use_measure,'.pdf')
pdf(file.path(FigDir,filename),8,5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

