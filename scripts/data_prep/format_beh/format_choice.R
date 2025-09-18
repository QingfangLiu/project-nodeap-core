
########################################################
# Prepare Day 2 Choice Data for Analysis
#
# This script combines choice data from the main dataset
# and NODEAP_17, recodes key variables, adds derived 
# columns (e.g., stimulation order per day), and saves
# the full dataset for downstream analysis.
#
# Final output:
# - A cleaned and organized choice data frame saved as RData
# - Sorted by subject ID for consistency
########################################################

rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

########################################################

# Load probe choice data
choice_dat <- read.xlsx(file.path(project_folder,'data_cleaned','Choices_all.xlsx'))

# Combine and recode variables
choice_dat <- choice_dat %>%
  mutate(
    PrePost = mapvalues(PrePost, from = c(1, 2), to = c('Pre', 'Post')),
    StimLoc = mapvalues(StimLoc, from = c('Posterior', 'Anterior'), to = c('pOFC', 'aOFC')),
    Set     = mapvalues(Set, from = c(1, 2), to = c('A', 'B')),
    StimOrder = mapvalues(StimOrder,
                          from = c('123','132','213','231','312','321'),
                          to   = c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                   'SC-SS-CS','SS-CS-SC','SS-SC-CS')),
    Cond = mapvalues(Cond, from = c('X','Y','Z'),
                     to = c('cTBS-sham','sham-cTBS','sham-sham'))
  ) %>%
  mutate(
    StimOrder_day1 = mapvalues(StimOrder,
                               from = c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                        'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                               to   = c('cTBS on 1st','cTBS on 1st',
                                        'cTBS on 2nd','cTBS on 3rd',
                                        'cTBS on 2nd','cTBS on 3rd')),
    StimOrder_day2 = mapvalues(StimOrder,
                               from = c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                        'SC-SS-CS','SS-CS-SC','SS-SC-CS'),
                               to   = c('cTBS on 2nd','cTBS on 3rd',
                                        'cTBS on 1st','cTBS on 1st',
                                        'cTBS on 3rd','cTBS on 2nd'))
  ) %>%
  mutate(
    PrePost         = factor(PrePost, levels = c('Pre', 'Post')),
    Cond            = factor(Cond, levels = c('sham-sham','cTBS-sham','sham-cTBS')),
    StimOrder_day1  = factor(StimOrder_day1),
    StimOrder_day2  = factor(StimOrder_day2),
    Set             = factor(Set),
    StimLoc         = factor(StimLoc, levels = c('aOFC', 'pOFC'))
  )

# Sort data by subject ID for consistency
choice_dat <- choice_dat %>%
  arrange(SubID)

# Save cleaned dataset
save(choice_dat, file = '../ProcessedData/choice_dat.RData')

