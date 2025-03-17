
rm(list = ls())
source('../Scripts_R/Setup.R')

opt = 'bySess'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'rerun_track_w_hier.rds'))

load(file = '../ProcessedData/Conditioning.RData') 
# this contains some NA trials
use_dat = conditioning_dat  
sub = use_dat$SubID
nsubs = length(unique(sub))
sess = use_dat$Sess
nsess = max(sess)

ncuepair = 12

# posterior values
alpha <- samples$BUGSoutput$sims.list$alpha
mu <- samples$BUGSoutput$sims.list$mu
kappa <- samples$BUGSoutput$sims.list$kappa
a <- samples$BUGSoutput$sims.list$a
k <- samples$BUGSoutput$sims.list$k

########## find the learned strength value at the end of learning

w <- samples$BUGSoutput$sims.list$w
dim(w)
use_w = w[,,,121,]    # the learned value in the end (121:120+1)
use_w_mean = apply(use_w,c(2,3,4),mean)


df_w = structure(use_w_mean,
                      .Dim = c(nsubs, nsess, ncuepair), 
                      .Dimnames = structure(list(SubID = unique(as.vector(sub)),
                                                 Sess = 1:nsess,
                                                 CuePair = 1:ncuepair
                                                 ), 
                                            .Names = c("SubID","Sess","CuePair")))
df_w = adply(df_w, c(1,2,3))
names(df_w) = c(names(df_w)[1:3],'w')

df_w = df_w %>%
  arrange(SubID,Sess)

s1=ggplot(df_w,aes(x=Sess,y=w,color=Sess,fill=Sess)) +
  geom_boxplot(alpha=0.4) +
  geom_jitter(alpha=0.5) +
    facet_wrap(~SubID) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(y='Learned stimulus value',x='Session')

pdf(file.path(FigDir,'Conditioning_w_subs_hier.pdf'),10,8)
print(s1)
dev.off()

# save df_w for later use during the choice phase on Day 2
write.csv(df_w,row.names = F,
          file = 'df_w_per_sub_sess_cuepair.csv')

#############
# examine the a parameter
alpha_map = apply(a,c(2,3),mean)

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


# avg two sham sessions for line plotting
df_alpha_summary = df_alpha %>%
  dplyr::group_by(TMS_types,SubID,StimLoc,StimOrder_day1) %>%
  reframe(alpha=mean(alpha))
ggplot(df_alpha_summary,aes(x=TMS_types,y=alpha,fill = TMS_types)) +
  facet_wrap(~StimOrder_day1,scales = 'free') +
  #facet_wrap(~StimLoc+StimOrder_day1,scales = 'free') +
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

