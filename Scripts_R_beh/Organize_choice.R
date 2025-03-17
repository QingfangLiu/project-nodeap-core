

rm(list = ls())
source('Setup.R')

###############################
# analyze choice data in the pre-meal and post-meal on Day 2
# update 6/30/24: add partial data from NODEAP_17
# update 12/3/24: added some cue-related columns to both xlsx files 
###############################

choice_dat1 = read.xlsx('../ProcessedData/Choices.xlsx')
choice_dat2 = read.xlsx('../ProcessedData/Choices_NODEAP_17.xlsx')

choice_dat = rbind(choice_dat1,choice_dat2) %>%
  mutate(PrePost=mapvalues(PrePost,from=c(1,2),to=c('Pre','Post')),
         StimLoc=mapvalues(StimLoc,from=c('Posterior','Anterior'),to=c('pOFC','aOFC')),
         Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('cTBS on 1st','cTBS on 1st',
                                         'cTBS on 2nd','cTBS on 3rd',
                                         'cTBS on 2nd','cTBS on 3rd')),
         StimOrder_day2 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('cTBS on 2nd','cTBS on 3rd','cTBS on 1st',
                                         'cTBS on 1st','cTBS on 3rd','cTBS on 2nd'))) %>%
  mutate(PrePost=factor(PrePost,levels=c('Pre','Post')),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Cond=factor(Cond),
         StimOrder_day1=factor(StimOrder_day1),
         StimOrder_day2=factor(StimOrder_day2),
         Set=factor(Set),
         StimLoc=factor(StimLoc,levels=c('aOFC','pOFC')))

# save the full choice data set
save(choice_dat,file = '../ProcessedData/choice_dat.RData')

