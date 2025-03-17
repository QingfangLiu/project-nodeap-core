
rm(list = ls())
source('Setup.R')
sub_dat = read.xlsx('../NODEAP_DataCollectionSheet.xlsx')
# clean those spaces
# need to do this once and save a new sheet
sub_dat$`Anterior/Posterior`=
  str_replace_all(sub_dat$`Anterior/Posterior`,
                  " ","")

conn_dat = read.xlsx('../ProcessedData/SeedStimConn_subs.xlsx') %>%
  mutate(stim=mapvalues(stim,from=c('aOFC stim','pOFC stim'),
                        to=c('aOFC stimulation region','pOFC stimulation region'))) %>%
  mutate(stim=factor(stim)) %>%
  mutate(StimLoc=rep(sub_dat$`Anterior/Posterior`,4))

aSeedaStim = subset(conn_dat,stim=='aOFC stimulation region'&
                      seed=='aOFC seed')$conn
pSeedaStim = subset(conn_dat,stim=='aOFC stimulation region'&
                      seed=='pOFC seed')$conn
sum(aSeedaStim > pSeedaStim)
wilcox.test(aSeedaStim,pSeedaStim,paired = T,alternative = 'g')

pSeedpStim = subset(conn_dat,stim=='pOFC stimulation region'&
                      seed=='pOFC seed')$conn
aSeedpStim = subset(conn_dat,stim=='pOFC stimulation region'&
                      seed=='aOFC seed')$conn
sum(pSeedpStim > aSeedpStim)
wilcox.test(pSeedpStim,aSeedpStim,paired = T,alternative = 'g')


# this shows the theoretical separation of the two networks
# across all 48 subjects
p=ggplot(conn_dat, aes(x = seed, y = conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~stim) + 
  labs(y='Functional connectivity (a.u.)',x='') +
  common

# use their group assignment
use_dat = subset(conn_dat,stim=='aOFC stimulation region'&
                   StimLoc=='Anterior'|
                   stim=='pOFC stimulation region'&
                   StimLoc=='Posterior')

p1=ggplot(use_dat, aes(x = seed, y = conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~stim) + 
  labs(y='Functional connectivity (a.u.)',x='') +
  common

pdf(file.path(FigDir,'Fig_conn_seed_stim.pdf'),7,5)
print(p)
dev.off()

pdf(file.path(FigDir,'Fig_conn_seed_stim_1.pdf'),7,5)
print(p1)
dev.off()


aSeedaStim = subset(use_dat,stim=='aOFC stimulation region'&
                      seed=='aOFC seed')$conn
pSeedaStim = subset(use_dat,stim=='aOFC stimulation region'&
                      seed=='pOFC seed')$conn
sum(aSeedaStim > pSeedaStim)
wilcox.test(aSeedaStim,pSeedaStim,paired = T,alternative = 'g')

pSeedpStim = subset(use_dat,stim=='pOFC stimulation region'&
                      seed=='pOFC seed')$conn
aSeedpStim = subset(use_dat,stim=='pOFC stimulation region'&
                      seed=='aOFC seed')$conn
sum(pSeedpStim > aSeedpStim)
wilcox.test(pSeedpStim,aSeedpStim,paired = T,alternative = 'g')



