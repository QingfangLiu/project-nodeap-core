

rm(list = ls())
source('Setup.R')
HomeDir = '/Users/liuq13/NODEAP'

## for calculating calories
# each subject, each D2 session
sess_names = c('S1D2','S2D2','S3D2')
sub_dat = NULL
sess_ctr = 0
Calories = matrix(NA,48,length(sess_names))
for(sess in sess_names){
  sess_ctr = sess_ctr + 1
tmp_dat = read.xlsx(file.path(HomeDir,'NODEAP_DataCollectionSheet_QL.xlsx'),sheet = sess)
Calories[,sess_ctr] = as.numeric(tmp_dat$Calories)

}

# put the calories with the sub cond together
sub_cond_dat = read.xlsx('../ProcessedData/SubConds.xlsx') %>%
  mutate(StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                    to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                         'SC-SS-CS','SS-CS-SC','SS-SC-CS')))

cbind(sub_cond_dat$SubID,Calories)

Calories[10,2] = NA # remove this session not included in the beh analysis

mean(Calories,na.rm=T)
sd(Calories,na.rm=T)/sqrt(48)

# overall of assignment
table(sub_cond_dat$StimLoc) # 23 aOFC, 25 pOFC
table(sub_cond_dat$Sex) # 32 females, 16 males
table(sub_cond_dat$Sex,sub_cond_dat$StimLoc) 

table(sub_cond_dat$StimOrder)
table(sub_cond_dat$Sex,sub_cond_dat$StimOrder)

# for each session
sess_names = c('S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
sub_dat = NULL
sess_ctr = 0
uncomfortable = matrix(NA,48,length(sess_names))
strong = matrix(NA,48,length(sess_names))
for(sess in sess_names){
  sess_ctr = sess_ctr + 1
  tmp_dat = read.xlsx(file.path(HomeDir,'NODEAP_DataCollectionSheet_QL.xlsx'),sheet = sess)
  uncomfortable[,sess_ctr] = as.numeric(tmp_dat$`Uncomfortable.TMS.(0-10)`)
  strong[,sess_ctr] = as.numeric(tmp_dat$`Strong.TMS.(0-10)`)
}

TMS_types = NULL
for(j in 1:nrow(sub_cond_dat)){
  tmp_order = sub_cond_dat$StimOrder[j]
  if(tmp_order==123) TMS_types[[j]] = c('C','S','S','C','S','S')
  if(tmp_order==132) TMS_types[[j]] = c('C','S','S','S','S','C')
  if(tmp_order==213) TMS_types[[j]] = c('S','C','C','S','S','S')
  if(tmp_order==231) TMS_types[[j]] = c('S','C','S','S','C','S')
  if(tmp_order==312) TMS_types[[j]] = c('S','S','C','S','S','C')
  if(tmp_order==321) TMS_types[[j]] = c('S','S','S','C','C','S')
}

TMS_rating_df = data.frame(SubID=rep(sub_cond_dat$SubID,each=6),
                           StimLoc = rep(sub_cond_dat$StimLoc,each=6),
                           sess_name = rep(sess_names,48),
                           TMSSess = rep(1:6,48), # TMS Sess from 1 to 6
                           TMS_type = unlist(TMS_types),
                           uncomfortable=as.numeric(t(uncomfortable)),
                           strong=as.numeric(t(strong)))
TMS_rating_df = TMS_rating_df %>%
  mutate(SubID=factor(SubID),
         StimLoc=factor(StimLoc),
         TMS_type=factor(TMS_type))

TMS_rating_df %>%
  group_by(TMS_type) %>%
  dplyr::summarise(uncomfortable=mean(uncomfortable,na.rm = T),
                   strong=mean(strong,na.rm = T))

model_1 <- lmer(uncomfortable ~ TMS_type + (1|SubID), data = TMS_rating_df)
model_0 <- lmer(uncomfortable ~ (1|SubID), data = TMS_rating_df)
anova(model_0,model_1)
summary(model_1)

model_1 <- lmer(strong ~ TMS_type + (1|SubID), data = TMS_rating_df)
model_0 <- lmer(strong ~ (1|SubID), data = TMS_rating_df)
anova(model_0,model_1)
summary(model_1)

TMS_rating_df %>%
  group_by(TMS_type,StimLoc) %>%
  dplyr::summarise(uncomfortable=mean(uncomfortable,na.rm = T),
                   strong=mean(strong,na.rm = T))

# cTBS only, does aOFC or pOFC matter?
use_dat = subset(TMS_rating_df,TMS_type=='C')
model_1 <- lmer(uncomfortable ~ StimLoc + (1|SubID), data = use_dat)
model_0 <- lmer(uncomfortable ~ (1|SubID), data = use_dat)
anova(model_0,model_1)

model_1 <- lmer(strong ~ StimLoc + (1|SubID), data = use_dat)
model_0 <- lmer(strong ~ (1|SubID), data = use_dat)
anova(model_0,model_1)

# standardize those ratings (z-score)
TMS_rating_df <- TMS_rating_df %>%
  group_by(SubID) %>%
  mutate(uncom = (uncomfortable - mean(uncomfortable, na.rm = TRUE)) / 
           sd(uncomfortable, na.rm = TRUE))  %>%
  mutate(strong = (strong - mean(strong, na.rm = TRUE)) / 
           sd(strong, na.rm = TRUE))
levels(TMS_rating_df$StimLoc)=c('aOFC','pOFC')
levels(TMS_rating_df$TMS_type)=c('cTBS','sham')

# people tend to rate cTBS as more uncomfortable and stronger
# than sham sessions, very clear pattern

p1=ggplot(TMS_rating_df,aes(x=TMS_type,y=uncom)) +
  facet_wrap(~StimLoc) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1.5, alpha = 0.8) +
  labs(y='Uncomfortable',x='') +
  common 

p2=ggplot(TMS_rating_df,aes(x=TMS_type,y=strong)) +
  facet_wrap(~StimLoc) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1.5, alpha = 0.8) +
  labs(y='Strong',x='') +
  common 

# are ratings of cTBS were influenced by which session they received it?
# not really

p3=TMS_rating_df %>%
  subset(TMS_type == 'cTBS') %>%
  ggplot(aes(x=TMSSess,y=uncom)) +
  geom_point() + 
  facet_wrap(~StimLoc,scales = 'free_y') +
  labs(y='Uncomfortable',x='') +
  common

p4=TMS_rating_df %>%
  subset(TMS_type == 'cTBS') %>%
  ggplot(aes(x=TMSSess,y=uncom)) +
  geom_point() + 
  facet_wrap(~StimLoc,scales = 'free_y') +
  labs(y='strong',x='') +
  common

p5=ggplot(TMS_rating_df,aes(x=strong,y=uncom)) +
  facet_wrap(~StimLoc) +
  geom_point(aes(color=TMS_type)) +
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
write.xlsx(TMS_rating_df,
           file = '../ProcessedData/Survey_uncomf_strong.xlsx')




