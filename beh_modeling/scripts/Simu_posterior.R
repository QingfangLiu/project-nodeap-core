
# this code tries to simulate the discrimination task across runs
# using the posterior estimates of alpha from rjags
# alpha is for each sess and each subj

rm(list=ls()) 

project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

library(rjags)
library(R2jags)

# load conditioning data to use
load(file = file.path(processed_dir,'Conditioning.RData'))

# this contains some NA trials
use_dat = conditioning_dat   

# prepare model inputs as vectors
sub = use_dat$SubID
sub = mapvalues(sub,from = unique(sub),to = 1:length(unique(sub))) # recode sub into consective numbers
nsubs = length(unique(sub))
sess = use_dat$Sess
nsess = max(sess)

ntrials_per_run = 24
nruns = 5
ntrials = ntrials_per_run * nruns

cuepair = use_dat$CuePair
ncuepair = length(unique(cuepair)) # 12 unique cue pairs

CuePair = Resp = array(NA,c(nsubs,nsess,ntrials))
SimuResp = array(NA,c(nsubs,nsess,ntrials))

for(j in 1:nsubs){
  for(c in 1:nsess){
    CuePair[j,c,] = use_dat$CuePair[sub==j & sess==c]
    Resp[j,c,] = use_dat$OdorChosen[sub==j & sess==c]
  }
}

## simulate the model

# load the MAP estimates of alpha
# from the session-specific learning model
alpha = read.csv(file = 'alpha_per_sub_sess.csv',row.names = NULL)
alpha = alpha[,c(2,3,4)]

nsim = 1e4
all_acc = array(NA,c(nsim,nsubs,nsess,nruns))

for(kk in 1:nsim){
  
  cat("Current iter:", kk, "\n")

w = array(NA,c(nsubs,nsess,(ntrials+1),ncuepair))
w[,,1,] = 1/2

for(j in 1:nsubs){
  for(c in 1:nsess){
    for(i in 1:ntrials){
      # update w of current cue pair or pass w of absent cue to next trial
      for(p in 1:ncuepair){
        w[j,c,i+1,p] <- w[j,c,i,p] + ifelse(CuePair[j,c,i]==p,
                                            alpha[j,c] * (1 - w[j,c,i,p]),0)
      }   
    } # end of trials loop
  } # end of sess loop
} # end of subs loop

# simulate response data
acc = array(NA,c(nsubs,nsess,5))
for(j in 1:nsubs){
  for(c in 1:nsess){
    for(i in 1:ntrials){
      SimuResp[j,c,i] = rbinom(1,1,w[j,c,i,CuePair[j,c,i]])
    }
    # separate trials into 5 blocks
    use_simu_resp = SimuResp[j,c,]
    use_simu_resp = matrix(use_simu_resp,24,5)
    acc[j,c,] = apply(use_simu_resp,2,mean)
    }}

all_acc[kk,,,] = acc

}

use_acc = apply(all_acc,c(2,3,4),mean)

pdf('Simu.pdf',8,7)
par(mfrow=c(2,2))
for(j in 1:nsubs){
matplot(t(use_acc[j,,]),type = 'l',ylab = 'Avg acc',
        lty = 1,lwd = 2,ylim = c(0.5,1),
        main = unique(use_dat$SubID)[j])
}
dev.off()

save.image(file = 'Simu_posterior.RData')

# load the data back & do some plotting

load(file = 'Simu_posterior.RData')

# convert array to df
all_acc = all_acc[1:100,,,]
all_acc = apply(all_acc,c(1,3,4),mean) # avg across subs

df_acc = structure(all_acc,
                        .Dim = c(100,nsess, nruns), 
                        .Dimnames = structure(list(sim = 1:100,
                                                   sess = 1:3, 
                                                   run = 1:5), 
                                              .Names = c("sim", "sess" ,"run")))
df_acc = adply(df_acc, c(1,2,3))
names(df_acc) = c(names(df_acc)[1:3],'acc')

# Summarize data: compute mean, lower/upper bounds for the ribbon
# quantify the variance across iterations, not across subjects
df_summary_by_sess <- df_acc %>%
  group_by(run, sess) %>%
  summarise(
    mean_acc = mean(acc),
    lower = quantile(acc, 0.025),  # 2.5th percentile for lower bound (95% CI)
    upper = quantile(acc, 0.975)
    ) 
    
# Now use the summarized data for plotting
ggplot(df_summary_by_sess, aes(x = run, y = mean_acc, fill = as.factor(sess),
                       group = sess, color = as.factor(sess))) +
  #geom_line(size = 1) +  # Add line for the mean accuracy
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(x = "Runs", y = "Simulated accuracy") +
  common
  

# to have TMS info
load(file = 'Simu_posterior.RData')
# convert array to df
use.iter = 10
all_acc = all_acc[1:use.iter,,,]

df_acc = structure(all_acc,
                   .Dim = c(use.iter,nsubs,nsess, nruns), 
                   .Dimnames = structure(list(sim = 1:use.iter,
                                              sub = 1:nsubs,
                                              sess = 1:3, 
                                              run = 1:5), 
                                         .Names = c("sim", "sub", "sess" ,"run")))
df_acc = adply(df_acc, c(1,2,3,4))
names(df_acc) = c(names(df_acc)[1:4],'acc')


SubInfo = read.xlsx('../../ProcessedData/SubConds.xlsx')
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
  if(tmp_order=='S-S-C') TMS_types[[j]] = rep(c('sham','sham','cTBS'),each=5*use.iter)
  if(tmp_order=='S-C-S') TMS_types[[j]] = rep(c('sham','cTBS','sham'),each=5*use.iter)
  if(tmp_order=='C-S-S') TMS_types[[j]] = rep(c('cTBS','sham','sham'),each=5*use.iter)
}


df_acc = df_acc %>%
  arrange(sub,sess) %>%
  mutate(StimLoc = factor(rep(SubInfo$StimLoc,each=5*3*use.iter)),
    TMS_types = factor(unlist(TMS_types)))
levels(df_acc$StimLoc) = c('aOFC','pOFC')

df_acc_avg_across_sub = df_acc %>%
  group_by(run,TMS_types,sim,StimLoc) %>%
  reframe(acc = mean(acc)) 

df_summary_by_TMS <- df_acc_avg_across_sub %>%
  group_by(run, TMS_types, StimLoc) %>%
  reframe(
    mean_acc = mean(acc),
    lower = quantile(acc, 0.025),  # 2.5th percentile for lower bound (95% CI)
    upper = quantile(acc, 0.975)
  ) 


ggplot(df_summary_by_TMS, aes(x = run, y = mean_acc, fill = as.factor(TMS_types),
                       group = TMS_types, color = as.factor(TMS_types))) +
  geom_line(size = 1) +  # Add line for the mean accuracy
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  scale_color_manual(values = use.col.conds.day1) +
  labs(x = "Runs", y = "Simulated accuracy") +
  common


strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

p1=conditioning_dat %>%
  group_by(SubID,Run,Cond_day1,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x = Run, y = Choice, color = Cond_day1, fill = Cond_day1)) +
  scale_color_manual(values = use.col.conds.day1) +
  scale_fill_manual(values = use.col.conds.day1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = "mean_se", linewidth = 1, 
               show.legend = F, width = 0.3,
               geom = "errorbar") +  
  facet_wrap2(~StimLoc,scales = 'free_y',strip = strip) +
  geom_ribbon(data = df_summary_by_TMS,
    aes(x = as.numeric(run), y = mean_acc,
        color = TMS_types,
        fill = TMS_types,
        group = TMS_types,
        ymin = lower, ymax = upper), alpha = 0.2, show.legend = F) +
  labs(x = "Runs", y = "P (Choosing rewarding stim)", 
       title = NULL,
       color = "Day1 TMS") +
  coord_cartesian(ylim = c(0.5, 1)) +
  common + theme(legend.position = c(0.85, 0.3)) 

# overall change of conditioning across sessions 1,2,3
p2=conditioning_dat %>%
  group_by(SubID,Run,Sess,StimLoc) %>%
  reframe(Choice=mean(OdorChosen,na.rm = T)) %>%
  ggplot(aes(x=Run,y=Choice,color=factor(Sess),
             fill=factor(Sess))) +
  stat_summary(fun = mean, geom = "line", show.legend = T,
               linewidth = 1, na.rm = T) +
  stat_summary(fun.data = "mean_se", show.legend = F,
               geom = "errorbar", width = 0.3, linewidth = 1) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  facet_wrap2(~StimLoc,scales = 'free_y',strip = strip) +
  geom_ribbon(data = df_summary_by_sess,
              aes(x = as.numeric(run), y = mean_acc,
                  color = sess,
                  fill = sess,
                  group = sess,
                  ymin = lower, ymax = upper), alpha = 0.2, show.legend = F) +
  labs(x = "Runs", y = "P(Choosing rewarding stim)",
       title = NULL) +
  coord_cartesian(ylim = c(0.5, 1)) +
  common + theme(legend.position = c(0.85, 0.3))


pdf(file.path(FigPaperDir,'Conditioning_acc_obs_sim.pdf'),8,8)
ggarrange(p1,p2,nrow = 2)
dev.off()





