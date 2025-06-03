
rm(list = ls())
source('../Scripts_R/Setup.R')

opt = 'bySess'
samples = readRDS(file = paste0('./JagsPosteriors/',opt,'rerun_track_w.rds'))

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

pdf(file.path(FigDir,'Conditioning_w_subs.pdf'),10,8)
print(s1)
dev.off()



