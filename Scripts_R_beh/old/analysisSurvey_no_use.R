
# this code analyzes data from the post-experiment survey
# so we have some self-report measure on subjects' experiences
# from the TMS

rm(list = ls())
source('Setup.R')

library(openxlsx)
survey_file = '/Users/liuq13/NODEAP/BonusQuestions.xlsx'
survey_dat = read.xlsx(survey_file)
# Convert columns to numeric
survey_dat[, 2:13] <- lapply(survey_dat[, 2:13], as.numeric)

uncom_dat = t(as.matrix(survey_dat[,seq(2,12,2)]))
strong_dat = t(as.matrix(survey_dat[,seq(3,13,2)]))

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]
nsubs = length(Subs)

SubInfo = SubInfo %>%
  mutate(StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS'))) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,
                                    from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                           'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('C-S-S','C-S-S','S-C-S','S-S-C','S-C-S','S-S-C')))


# add a column called TMS_types indicating the TMS type of each session
TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order=='CS-SC-SS') TMS_types[[j]] = rep(c('C','S','S','C','S','S'),1)
  if(tmp_order=='CS-SS-SC') TMS_types[[j]] = rep(c('C','S','S','S','S','C'),1)
  if(tmp_order=='SC-CS-SS') TMS_types[[j]] = rep(c('S','C','C','S','S','S'),1)
  if(tmp_order=='SC-SS-CS') TMS_types[[j]] = rep(c('S','C','S','S','C','S'),1)
  if(tmp_order=='SS-CS-SC') TMS_types[[j]] = rep(c('S','S','C','S','S','C'),1)
  if(tmp_order=='SS-SC-CS') TMS_types[[j]] = rep(c('S','S','S','C','C','S'),1)
}

df_survey = data.frame(TMS_types = factor(unlist(TMS_types)),
                      TMSSess = rep(1:6,nsubs), # TMS Sess from 1 to 6
                      Sess = rep(c('S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'),nsubs),
                      uncom = as.numeric(uncom_dat),
                      strong = as.numeric(strong_dat),
                      StimLoc = rep(SubInfo$StimLoc,each=6),
                      SubID = rep(SubInfo$SubID,each=6)) %>%
mutate(SubID=factor(SubID),
       StimLoc=factor(StimLoc))

# standardize those ratings (z-score)
df_survey <- df_survey %>%
  group_by(SubID) %>%
  mutate(uncom = (uncom - mean(uncom, na.rm = TRUE)) / 
           sd(uncom, na.rm = TRUE))  %>%
  mutate(strong = (strong - mean(strong, na.rm = TRUE)) / 
           sd(strong, na.rm = TRUE))
levels(df_survey$StimLoc)=c('aOFC','pOFC')
levels(df_survey$TMS_types)=c('cTBS','sham')

# people tend to rate cTBS as more uncomfortable and stronger
# than sham sessions, very clear pattern

p1=ggplot(df_survey,aes(x=TMS_types,y=uncom)) +
  facet_wrap(~StimLoc) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1.5, alpha = 0.8) +
  labs(y='Uncomfortable',x='') +
  common 

p2=ggplot(df_survey,aes(x=TMS_types,y=strong)) +
  facet_wrap(~StimLoc) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1.5, alpha = 0.8) +
  labs(y='Strong',x='') +
  common 

# whether ratings of cTBS were influenced by which session they received it
# not really

p3=df_survey %>%
  subset(TMS_types == 'cTBS') %>%
  ggplot(aes(x=TMSSess,y=uncom)) +
  geom_point() + 
  facet_wrap(~StimLoc,scales = 'free_y') +
  labs(y='Uncomfortable',x='') +
  common

p4=df_survey %>%
  subset(TMS_types == 'cTBS') %>%
  ggplot(aes(x=TMSSess,y=uncom)) +
  geom_point() + 
  facet_wrap(~StimLoc,scales = 'free_y') +
  labs(y='strong',x='') +
  common

p5=ggplot(df_survey,aes(x=strong,y=uncom)) +
  facet_wrap(~StimLoc) +
  geom_point(aes(color=TMS_types)) +
  labs(y='strong',x='Uncomfortable') +
  common

pdf(file.path(FigDir,'Survey_uncomf_strong.pdf'),6,4)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()

# save survey data frame
write.xlsx(df_survey,file = '../ProcessedData/Survey_uncomf_strong.xlsx')




