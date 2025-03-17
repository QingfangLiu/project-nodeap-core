
# this aggregates motion results across sessions

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load SNRs vals from matlab
matfile = '../ProcessedData/motion_subs_sessions.mat'
matdat = readMat(matfile)
all.ctr.badvol = matdat$all.ctr.badvol
all.fd = matdat$all.fd

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

df_motion = structure(all.fd,
                       .Dim = c(n_rest_names, length(Subs)), 
                       .Dimnames = structure(list(Sess = rest_names, 
                                                  SubID = as.vector(Subs)), 
                                             .Names = c("Sess","SubID")))
df_motion = adply(df_motion, c(1,2))
names(df_motion) = c(names(df_motion)[1:2],'motion')
df_motion$StimLoc = NA
df_motion$StimOrder = NA

for(subj in Subs){
  df_motion$StimLoc[df_motion$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
  df_motion$StimOrder[df_motion$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

df_motion$StimOrder=mapvalues(df_motion$StimOrder,
                    from=c('123','132','213','231','312','321'),
                    to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                         'SC-SS-CS','SS-CS-SC','SS-SC-CS'))

# add a column called TMS_types indicating the TMS type of each session
TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order==123) TMS_types[[j]] = c('N','C','S','S','C','S','S')
  if(tmp_order==132) TMS_types[[j]] = c('N','C','S','S','S','S','C')
  if(tmp_order==213) TMS_types[[j]] = c('N','S','C','C','S','S','S')
  if(tmp_order==231) TMS_types[[j]] = c('N','S','C','S','S','C','S')
  if(tmp_order==312) TMS_types[[j]] = c('N','S','S','C','S','S','C')
  if(tmp_order==321) TMS_types[[j]] = c('N','S','S','S','C','C','S')
}

df_motion$TMS_types = factor(unlist(TMS_types))
levels(df_motion$TMS_types) = c('cTBS','null','sham')

use_df = df_motion %>%
  subset(TMS_types %in% c('cTBS','sham')) %>% # only include cTBS & sham
  mutate(time = as.numeric(factor(Sess,levels = unique(Sess)))) # add time variable

# how motion varies with TMS types?
p1=use_df %>%
  ggplot(aes(x=TMS_types,y=motion,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = 'Framewise Displacement (mm)', x = '') +
  common

use_df %>%
  group_by(TMS_types,StimLoc) %>%
  summarise(mean=mean(motion,na.rm=T),
            sd=sd(motion,na.rm=T))

use_df %>%
  group_by(TMS_types) %>%
  summarise(motion=mean(motion,na.rm=T))

model0 <- lmer(motion ~ (1|SubID), data = use_df)
model1 <- lmer(motion ~ TMS_types + (1|SubID), data = use_df)
anova(model1,model0) 
summary(model1)

model2 <- lmer(motion ~ TMS_types + StimLoc + (1|SubID), data = use_df)
anova(model2,model0) 

p1a=use_df %>%
  ggplot(aes(x=TMS_types,y=motion,fill=TMS_types)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=TMS_types,group=TMS_types), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  facet_wrap(~StimLoc) +
  labs(title = '',y = 'Framewise Displacement (mm)', x = '') +
  common

# how motion changes with sessions
p2=use_df %>%
  ggplot(aes(x=Sess,y=motion,fill=Sess)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color=Sess,group=Sess), 
              position = position_jitterdodge(dodge.width = 0.6, 
                                              jitter.width = 0.3),
              size = 0.5, alpha = 0.4) +
  labs(title = '',y = 'Framewise Displacement (mm)', x = '') +
  common
  
pdf(file.path(FigDir,'MRI_motion_sessions.pdf'),8,7)
print(p1)
print(p1a)
print(p2)
dev.off()

