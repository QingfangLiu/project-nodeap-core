

opt = 'bySess'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'.rds'))
dic_value <- samples$BUGSoutput$DIC
pd_value <- samples$BUGSoutput$pD
cat("DIC:", dic_value, "\n")
cat("pD (Effective Number of Parameters):", pd_value, "\n")

opt = 'same'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'.rds'))
dic_value <- samples$BUGSoutput$DIC
pd_value <- samples$BUGSoutput$pD
cat("DIC:", dic_value, "\n")
cat("pD (Effective Number of Parameters):", pd_value, "\n")

opt = 'byTMS'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'.rds'))
dic_value <- samples$BUGSoutput$DIC
pd_value <- samples$BUGSoutput$pD
cat("DIC:", dic_value, "\n")
cat("pD (Effective Number of Parameters):", pd_value, "\n")
# continue with this byTMS model only

#################################

# posterior values
alpha <- samples$BUGSoutput$sims.list$alpha
mu <- samples$BUGSoutput$sims.list$mu
kappa <- samples$BUGSoutput$sims.list$kappa


apply(mu,2,mean)

alpha_map = apply(alpha,c(2,3),mean)
write.csv(alpha_map,file = 'alpha.csv')

boxplot(alpha_map)

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
                      alpha = as.numeric(t(apply(alpha,c(2,3),mean))),
                      SubID = rep(SubInfo$SubID,each=3))

model0 <- lmer(alpha ~ Sess + (1|SubID), data = df_alpha)
summary(model0)

model1 <- lmer(alpha ~ Sess + TMS_types + (1|SubID), data = df_alpha)
anova(model1,model0) 
summary(model1)


df_alpha_summary = df_alpha %>%
  dplyr::group_by(TMS_types,SubID) %>%
  summarise(alpha=mean(alpha))

par(ask=F)
ggplot(df_alpha_summary,aes(x=TMS_types,y=alpha)) +
  geom_boxplot(width = 0.6, alpha = 0.4, position = position_dodge(0.75),
               outlier.alpha = 0) +
  geom_jitter(aes(group=SubID), position = pd, 
              size = 1, alpha = 0.8) +
  geom_line(aes(group=SubID), position = pd, linewidth = 0.5, color = 'darkgray') +
  common


