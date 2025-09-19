
rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

###############################
#### analyze conditioning data ####
###############################

conditioning_dat1 = read.xlsx('../ProcessedData/Conditioning.xlsx')
# not include NODEAP_17 due to lack of cTBS-sham session

conditioning_dat = conditioning_dat1 %>%
  mutate(SelectedCue = if_else(SelectedCue == 0, NA_integer_, SelectedCue)) %>% # replace zero with NA
  mutate(Set=mapvalues(Set,from=c(1,2),to=c('A','B')),
         ChoiceType=mapvalues(ChoiceType,from=c(1,2),to=c('Sweet','Savory')),
         StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                        to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                             'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
         Cond=mapvalues(Cond,from=c('X','Y','Z'),
                        to=c('cTBS-sham','sham-cTBS','sham-sham'))) %>%
  mutate(Cond_day1 = ifelse(Cond=='cTBS-sham','cTBS','sham'),
         Cond_day2 = ifelse(Cond=='sham-cTBS','cTBS','sham')) %>%
  mutate(StimOrder_day1 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('C-S-S','C-S-S','S-C-S','S-S-C','S-C-S','S-S-C')),
         StimOrder_day2 = mapvalues(StimOrder,from=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                                     'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                                    to=c('S-C-S','S-S-C','C-S-S','C-S-S','S-S-C','S-C-S'))) %>%
  mutate(SubID=factor(SubID),
         StimOrder=factor(StimOrder),
         StimOrder_day1=factor(StimOrder_day1),
         StimOrder_day2=factor(StimOrder_day2),
         Cond=factor(Cond,levels=c('sham-sham','cTBS-sham','sham-cTBS')),
         Set=factor(Set),
         Cond_day1 = factor(Cond_day1),
         Cond_day2 = factor(Cond_day2),
         ChoiceType=factor(ChoiceType),
         StimLoc=factor(StimLoc))

levels(conditioning_dat$StimLoc)=c('aOFC','pOFC')

levels(conditioning_dat$StimOrder_day1)=
  c('cTBS on 1st','cTBS on 2nd','cTBS on 3rd')

save(conditioning_dat,file = '../ProcessedData/Conditioning.RData')

