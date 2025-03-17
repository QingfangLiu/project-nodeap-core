

# this tries to correlate behavioral effects of Day 2 with FC measures
# from resting-state data

rm(list = ls())
source('Setup.R')

SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx')
Subs = SubInfo$SubID[SubInfo$Include==1]

# load FC_Day2 vals from matlab
matfile = '../ProcessedData/FC_Day2_subs_sessions_w_shift_new.mat'
matdat = readMat(matfile)
FCdat_Day2 = matdat$FCdat.Day2

rest_names = c('S1D2','S2D2','S3D2')
n_rest_names = length(rest_names)

aal_name_file = '../Scripts_matlab/FuncConn_AAL/TKsent/aal.mat'
aal_name = unlist(readMat(aal_name_file)) 
num_rois = length(aal_name)

mask_labels = c('aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim')
n_masks = length(mask_labels)

# read in behavioral choice data
load('../ProcessedData/Summary_Choice_corrected_dat.RData')
summary_choice_corrected = summary_choice_corrected %>%
  arrange(SubID,Sess)

# convert this array to df
dfFC_Day2 = structure(FCdat_Day2,
                        .Dim = c(num_rois, n_masks, n_rest_names, length(Subs)), 
                        .Dimnames = structure(list(roi = aal_name,
                                                   mask = mask_labels,
                                                   Sess = rest_names, 
                                                   SubID = as.vector(Subs)), 
                                              .Names = c("roi", "mask", "Sess" ,"SubID")))
dfFC_Day2 = adply(dfFC_Day2, c(1,2,3,4))
names(dfFC_Day2) = c(names(dfFC_Day2)[1:4],'FC')

dfFC_Day2 = dfFC_Day2 %>%
  subset(!(SubID=='NODEAP_17' & Sess=='S2D2')) %>%
  arrange(SubID,Sess) %>%
  mutate(Cond = rep(summary_choice_corrected$Cond,
                    each=num_rois*n_masks)) %>%
  mutate(beheffect = rep(summary_choice_corrected$ChoiceChangeAB,
                         each=num_rois*n_masks)) %>%
  mutate(StimLoc=rep(summary_choice_corrected$StimLoc,
                     each=num_rois*n_masks)) %>%
  subset(StimLoc=='pOFC') # tried focusing on pOFC subjects, or both

# correlate behavioral effect and FCs
# to see if there's any potential correlation
# fully exploratory
for (r in num_rois){
  for(m in n_masks){
    tmp_dat = subset(dfFC_Day2,roi==aal_name[r] & mask==mask_labels[m]) 
    tmp_corr_test = cor.test(tmp_dat$FC,tmp_dat$beheffect)
    if(tmp_corr_test$p.value < 0.1){
      print(paste('Found',aal_name[r],mask_labels[m],
                  tmp_corr_test$p.value))
    }
  }
}

dfFC_Day2 %>%
  subset(roi==aal_name[2] & mask==mask_labels[3]) %>%
  ggplot(aes(x=beheffect,y=FC)) +
  geom_point() +
  #facet_wrap(~Cond) +
  geom_smooth(method = "lm", se = T, color = 'black') +
  stat_cor(method = 'pearson',size = 5) +
  common


# look at sham-sham vs. sham-cTBS difference, indicating
# TMS effect on behavioral and neural

dfFC_Day2_diff = dfFC_Day2 %>%
  subset(Cond %in% c('sham-sham','sham-cTBS'))

df_sham_sham = subset(dfFC_Day2_diff,Cond=='sham-sham')
df_sham_cTBS = subset(dfFC_Day2_diff,Cond=='sham-cTBS')

dfFC_Day2_diff_new = dfFC_Day2 %>%
  subset(Cond == 'sham-sham') %>%
  mutate(beheffect_diff = df_sham_cTBS$beheffect - df_sham_sham$beheffect) %>%
  mutate(FC_diff = df_sham_cTBS$FC - df_sham_sham$FC) 
dfFC_Day2_diff_new$FC = NULL
dfFC_Day2_diff_new$beheffect = NULL
dfFC_Day2_diff_new$Cond = NULL
dfFC_Day2_diff_new$Sess = NULL


for (r in num_rois){
  for(m in n_masks){
    tmp_dat = subset(dfFC_Day2_diff_new,roi==aal_name[r] & mask==mask_labels[m]) 
    tmp_corr_test = cor.test(tmp_dat$FC_diff,tmp_dat$beheffect_diff)
    if(tmp_corr_test$p.value < 0.2){
      print(paste('Found',aal_name[r],mask_labels[m],
                  tmp_corr_test$p.value))
    }
  }
}


dfFC_Day2_diff_new %>%
  subset(roi==aal_name[2] & mask==mask_labels[3]) %>%
  ggplot(aes(x=beheffect_diff,y=FC_diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = 'black') +
  stat_cor(method = 'pearson',size = 5) +
  common

