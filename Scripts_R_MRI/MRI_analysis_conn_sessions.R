
# this aggregates conn results across sessions
# b/t seed OFC region and target LPFC region
# based on the stimulation location

rm(list = ls())
source('../Scripts_R_beh/Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load conn vals from matlab
matdat = readMat('../ProcessedData/SeedStimConn_subs_sessions.mat')
convals = matdat$mean.convals

rest_names = c('D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
n_rest_names = length(rest_names)

# convert this array to df
dfconvals = structure(convals,
                       .Dim = c(2L, 2L, length(Subs), 7), 
                       .Dimnames = structure(list(Seed = c('aOFC','pOFC'), 
                                                  Stim = c('aOFC','pOFC'), 
                                                  SubID = as.vector(Subs),
                                                  Sess = rest_names), 
                                             .Names = c("Seed", "Stim" ,"SubID", "Sess")))
dfconvals = adply(dfconvals, c(1,2,3,4))
names(dfconvals) = c(names(dfconvals)[1:4],'Conn')
dfconvals$StimLoc = NA
dfconvals$StimOrder = NA

for(subj in Subs){
dfconvals$StimLoc[dfconvals$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
dfconvals$StimOrder[dfconvals$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfconvals$StimOrder=mapvalues(dfconvals$StimOrder,
                    from=c('123','132','213','231','312','321'),
                    to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                         'SC-SS-CS','SS-CS-SC','SS-SC-CS'))
dfconvals = dfconvals %>%
  arrange(SubID) # reorder the data frame for convenience of next step

levels(dfconvals$Seed)=c('Seed-aOFC','Seed-pOFC')
levels(dfconvals$Stim)=c('aOFC-conn-LPFC','pOFC-conn-LPFC')

# add a column called TMS_types indicating the TMS type of each session
# Day 0 is labeled as null
TMS_types = NULL
for(j in 1:nrow(SubInfo)){
  tmp_order = SubInfo$StimOrder[j]
  if(tmp_order==123) TMS_types[[j]] = rep(c('N','C','S','S','C','S','S'),each=4)
  if(tmp_order==132) TMS_types[[j]] = rep(c('N','C','S','S','S','S','C'),each=4)
  if(tmp_order==213) TMS_types[[j]] = rep(c('N','S','C','C','S','S','S'),each=4)
  if(tmp_order==231) TMS_types[[j]] = rep(c('N','S','C','S','S','C','S'),each=4)
  if(tmp_order==312) TMS_types[[j]] = rep(c('N','S','S','C','S','S','C'),each=4)
  if(tmp_order==321) TMS_types[[j]] = rep(c('N','S','S','S','C','C','S'),each=4)
}

dfconvals$TMS_types = factor(unlist(TMS_types))
levels(dfconvals$TMS_types) = c('cTBS','null','sham')


p0=dfconvals %>%
  ggplot(aes(x = Seed, y = Conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = interaction(SubID,Sess)), color = "gray", alpha = 0.2,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.1, 
              position = position_dodge(0.2)) +
  facet_wrap(~Stim) + 
  labs(y='Functional connectivity (a.u.)',x='',
       title = 'All sessions') +
  common

# use their group assignment
dfconvals_expt = subset(dfconvals,Stim=='aOFC-conn-LPFC'&
                   StimLoc=='Anterior'|
                   Stim=='pOFC-conn-LPFC'&
                   StimLoc=='Posterior')

p1=dfconvals_expt %>%
  ggplot(aes(x = Seed, y = Conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = interaction(SubID,Sess)), color = "gray", alpha = 0.2,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.1, 
              position = position_dodge(0.2)) +
  facet_wrap(~Stim) + 
  labs(y='Functional connectivity (a.u.)',x='',
       title = 'All sessions') +
  common

p1_D0=dfconvals_expt %>%
  subset(Sess=='D0') %>%
  ggplot(aes(x = Seed, y = Conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = interaction(SubID,Sess)), color = "gray", alpha = 0.5,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~Stim) + 
  labs(y='Functional connectivity (a.u.)',x='',
       title = 'Day 0') +
  common


use_dat = dfconvals_expt %>%
  subset(Sess=='D0')
aSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-aOFC')$Conn
pSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-pOFC')$Conn
aSeedaStim = aSeedaStim[!is.na(aSeedaStim)]
pSeedaStim = pSeedaStim[!is.na(pSeedaStim)]
sum(aSeedaStim > pSeedaStim)
length(aSeedaStim)
wilcox.test(aSeedaStim,pSeedaStim,paired = T,alternative = 'g')

pSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-pOFC')$Conn
aSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-aOFC')$Conn
pSeedpStim = pSeedpStim[!is.na(pSeedpStim)]
aSeedpStim = aSeedpStim[!is.na(aSeedpStim)]
sum(pSeedpStim > aSeedpStim)
length(pSeedpStim)
wilcox.test(pSeedpStim,aSeedpStim,paired = T,alternative = 'g')

# look at the change across sham, null and cTBS
p2=dfconvals_expt %>%
  ggplot(aes(x = Seed, y = Conn, fill = TMS_types)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = TMS_types), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), 
              alpha = 0.2) +
  facet_wrap(~Stim) + 
  labs(y='Functional connectivity (a.u.)',x='') +
  common

summary_convals_expt = dfconvals_expt %>%
  group_by(Seed,Stim,SubID,StimLoc) %>%
  summarize(Conn=mean(Conn,na.rm=T)) %>%
  ungroup() %>%
  arrange(SubID)

p3=summary_convals_expt %>%
  ggplot(aes(x = Seed, y = Conn)) +
  geom_violinhalf(fill = 'blue', alpha = 0.5, flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~Stim) + 
  labs(y='Functional connectivity (a.u.)',x='',
       title = 'Avg across sessions') +
  common


use_dat = summary_convals_expt
aSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-aOFC')$Conn
pSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-pOFC')$Conn
sum(aSeedaStim > pSeedaStim)
length(aSeedaStim)
wilcox.test(aSeedaStim,pSeedaStim,paired = T,alternative = 'g')

pSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-pOFC')$Conn
aSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-aOFC')$Conn
sum(pSeedpStim > aSeedpStim)
length(pSeedpStim)
wilcox.test(pSeedpStim,aSeedpStim,paired = T,alternative = 'g')


summary_convals = dfconvals %>%
  group_by(Seed,Stim,SubID,StimLoc) %>%
  summarize(Conn=mean(Conn,na.rm=T)) %>%
  ungroup() %>%
  arrange(SubID)


strip = strip_themed(background_x = elem_list_rect(fill = use.col.ap.lpfc),
                     text_x = elem_list_text(color = 'black',
                                             face = "italic",
                                             size = 16))

library(grid)
muh_grob <- grid::rectGrob(
  x=1:2, y=0, gp=gpar(face='bold',
    color='black', fill=use.col.ap.ofc, alpha=1))

# includes all ppts in each panel
pp=summary_convals %>%
  ggplot(aes(x = Seed, y = Conn)) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_violinhalf(flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap2(~Stim,axes = 'all',strip = strip) + 
  labs(y='Functional connectivity (a.u.)',x=NULL,title=NULL) +
  annotate("text", x=1.5, y=15.8, label='***') +
  annotate("segment",x = 1, xend = 1, y = 15.5, yend = 15) +
  annotate("segment",x = 2, xend = 2, y = 15.5, yend = 15) +
  annotate("segment",x = 1, xend = 2, y = 15.5, yend = 15.5) +
  common +
  coord_cartesian(clip='off') +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  theme(axis.text.x = element_text(margin=margin(t=10),color = 'white')) +
  annotation_custom(
    grob=muh_grob, xmin = 0, xmax = 1, ymin = -4.4, ymax= -2.7)


use_dat = summary_convals
aSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-aOFC')$Conn
pSeedaStim = subset(use_dat,Stim=='Stim-aOFC'&
                      Seed=='Seed-pOFC')$Conn
sum(aSeedaStim > pSeedaStim)
length(aSeedaStim)
wilcox.test(aSeedaStim,pSeedaStim,paired = T,alternative = 't')

pSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-pOFC')$Conn
aSeedpStim = subset(use_dat,Stim=='Stim-pOFC'&
                      Seed=='Seed-aOFC')$Conn
sum(pSeedpStim > aSeedpStim)
length(pSeedpStim)
wilcox.test(pSeedpStim,aSeedpStim,paired = T,alternative = 't')

pdf(file.path(FigPaperDir,'Fig_conn_seed_stim_sessions_paper.pdf'),7,4)
print(pp)
dev.off()
  

pdf(file.path(FigDir,'Fig_conn_seed_stim_sessions.pdf'),8,5)
print(p0)
print(p1)
print(p1_D0)
print(p2)
print(p3)
dev.off()

# calculate indexes
df1 = summary_convals_expt %>% subset(Seed=='Seed-aOFC')
ACI = df1$Conn
df2 = summary_convals_expt %>% subset(Seed=='Seed-pOFC')
PCI = df2$Conn

df_idx = df1 %>%
  select(SubID,StimLoc) %>%
  mutate(ACI=ACI,PCI=PCI,OCI=ACI+PCI,RCI=ACI-PCI,ARCI=abs(RCI))

write.csv(df_idx, "../ProcessedData/ConnIdx.csv", 
          row.names = FALSE)

