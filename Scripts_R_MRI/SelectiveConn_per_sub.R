

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

tmp_path = '/Volumes/X9Pro/NODEAP/ConnMap_PAID_global/Agg'
n_voxels_vec = seq(50,500,50)
opts = c('raw','abs')

roinames = c('aOFC_seed','aOFC_stim', 'pOFC_seed','pOFC_stim');
n_roi = length(roinames)

useFigDir = paste0(FigDir,'SelectedConn')
if (!dir.exists(useFigDir)) {
  dir.create(useFigDir, recursive = TRUE)
} 

for (n_top_voxels in n_voxels_vec){
  for (opt in opts){
  
  mat_name = paste0('FC_corr_sham_cTBS_top',n_top_voxels,'_per_sub.mat')
  matfile = file.path(tmp_path,mat_name)
  matdat = readMat(matfile)
  if(opt=='raw') SelectedConn = matdat$mean.corr.map.TMS
  if(opt=='abs') SelectedConn = matdat$abs.mean.corr.map.TMS
    
# convert this array to df
dfSelectedConn = structure(SelectedConn,
                   .Dim = c(length(Subs), n_roi, 2), 
                   .Dimnames = structure(list(SubID = as.vector(Subs),
                                              roi = roinames,
                                              TMS = c('sham','cTBS')), 
                                         .Names = c("SubID","roi","TMS")))
dfSelectedConn = adply(dfSelectedConn, c(1,2,3))
names(dfSelectedConn) = c(names(dfSelectedConn)[1:3],'Conn')
dfSelectedConn$StimLoc = NA
dfSelectedConn$StimOrder = NA

for(subj in Subs){
  dfSelectedConn$StimLoc[dfSelectedConn$SubID==subj] = SubInfo$StimLoc[SubInfo$SubID==subj]
  dfSelectedConn$StimOrder[dfSelectedConn$SubID==subj] = SubInfo$StimOrder[SubInfo$SubID==subj]
}

dfSelectedConn$StimOrder=mapvalues(dfSelectedConn$StimOrder,
                           from=c('123','132','213','231','312','321'),
                           to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                'SC-SS-CS','SS-CS-SC','SS-SC-CS'))

p0=dfSelectedConn %>%
  ggplot(aes(x = TMS, y = Conn, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = '',y = 'Conn (raw)', x = '') +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test", 
                     alternative = "greater",
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.signif", 
                     hide.ns = F) + 
  common

p1=dfSelectedConn %>%
  subset(StimLoc=='Anterior') %>%
  ggplot(aes(x = TMS, y = Conn, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'aOFC stim only',y = 'Conn (raw)', x = '') +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test", 
                     alternative = "greater",
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.signif", 
                     hide.ns = F) + 
  common

p2=dfSelectedConn %>%
  subset(StimLoc=='Posterior') %>%
  ggplot(aes(x = TMS, y = Conn, color = TMS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap(~roi,ncol = 4) +
  labs(title = 'pOFC stim only',y = 'Conn (raw)', x = '') +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test", 
                     alternative = "greater",
                     comparisons = list(c("sham", "cTBS")), 
                     label = "p.signif", 
                     hide.ns = F) + 
  common

fig_name = paste0('Selected_conn_top',n_top_voxels,'_',opt,'_per_sub.pdf')
pdf(file.path(useFigDir,fig_name),8,5)
print(p0)
print(p1)
print(p2)
dev.off()

}
}


