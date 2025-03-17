
# this aggregates tSNRs PAID results across sessions

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load SNRs vals from matlab
matfile = '../ProcessedData/tSNRs_PAID_subs_sessions.mat'
matdat = readMat(matfile)
MeantSNRs = matdat$MeantSNRs

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

roinames = c('aOFC_seed','pOFC_seed','aOFC_stim','pOFC_stim','gm_0.1','wm_0.9','csf_0.9');
n_roi = length(roinames)

# convert this array to df
dfMeantSNRs = structure(MeantSNRs,
                       .Dim = c(n_rest_names, n_roi, length(Subs)), 
                       .Dimnames = structure(list(Sess = rest_names, 
                                                  roi = roinames,
                                                  SubID = as.vector(Subs)), 
                                             .Names = c("Sess", "roi" ,"SubID")))
dfMeantSNRs = adply(dfMeantSNRs, c(1,2,3))
names(dfMeantSNRs) = c(names(dfMeantSNRs)[1:3],'SNRs')
dfMeantSNRs$StimLoc = NA
dfMeantSNRs$StimOrder = NA

for(subj in Subs){
  dfMeantSNRs$StimLoc[dfMeantSNRs$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
  dfMeantSNRs$StimOrder[dfMeantSNRs$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfMeantSNRs$StimOrder=mapvalues(dfMeantSNRs$StimOrder,
                    from=c('123','132','213','231','312','321'),
                    to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                         'SC-SS-CS','SS-CS-SC','SS-SC-CS'))

# add a column called TMS_types indicating the TMS type of each session

TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order==123) TMS_types[[j]] = rep(c('N','C','S','S','C','S','S'),n_roi)
  if(tmp_order==132) TMS_types[[j]] = rep(c('N','C','S','S','S','S','C'),n_roi)
  if(tmp_order==213) TMS_types[[j]] = rep(c('N','S','C','C','S','S','S'),n_roi)
  if(tmp_order==231) TMS_types[[j]] = rep(c('N','S','C','S','S','C','S'),n_roi)
  if(tmp_order==312) TMS_types[[j]] = rep(c('N','S','S','C','S','S','C'),n_roi)
  if(tmp_order==321) TMS_types[[j]] = rep(c('N','S','S','S','C','C','S'),n_roi)
}

dfMeantSNRs$TMS_types = factor(unlist(TMS_types))
levels(dfMeantSNRs$TMS_types) = c('cTBS','null','sham')


p1=dfMeantSNRs %>%
  ggplot(aes(x=roi,y=SNRs,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = 'tSNRs', x = '') +
  common
  
p2=dfMeantSNRs %>%
  ggplot(aes(x=roi,y=SNRs,fill=Sess)) +
  geom_boxplot(width = 0.8, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=Sess), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = 'tSNRs', x = '') +
  common

# add time variable
dfMeantSNRs$time <- as.numeric(factor(dfMeantSNRs$Sess, 
                                     levels = unique(dfMeantSNRs$Sess)))
use.dat = subset(dfMeantSNRs,roi=='pOFC_seed') # test on each roi
model_1 <- lmer(SNRs ~ time + (1|SubID), data = use.dat)
model_0 <- lmer(SNRs ~ (1|SubID), data = use.dat)
anova(model_0,model_1)
summary(model_1)


pdf(file.path(FigDir,'tSNR_sessions.pdf'),12,5)
print(p1)
print(p2)
dev.off()

