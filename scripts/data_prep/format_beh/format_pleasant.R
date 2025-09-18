
rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# This script loads and preprocesses pleasantness rating data,
# recodes relevant variables for analysis, and prepares the data
# for plotting and further statistical modeling.

# 📥 Load data
dat0 <- read.xlsx('../../data_input/Pleasant_all.xlsx')

# 🧹 Combine and clean
dat <- dat0 %>%
  # Recode variables
  mutate(
    PrePost = mapvalues(PrePost, from = c(1, 2), to = c('Pre-meal', 'Post-meal')),
    IfDevalue = mapvalues(IfDevalue, from = c('0', '1'), to = c('non-sated', 'sated')),
    StimLoc = mapvalues(StimLoc, from = c('Anterior', 'Posterior'), to = c('aOFC', 'pOFC')),
    Cond = mapvalues(Cond, from = c('X', 'Y', 'Z'), to = c('cTBS-sham', 'sham-cTBS', 'sham-sham'))
  ) %>%
  # Derive day-specific condition
  mutate(
    Cond_day1 = ifelse(Cond == 'cTBS-sham', 'cTBS', 'sham'),
    Cond_day2 = ifelse(Cond == 'sham-cTBS', 'cTBS', 'sham')
  ) %>%
  # Set factor levels
  mutate(
    PrePost = factor(PrePost, levels = c('Pre-meal', 'Post-meal')),
    Cond = factor(Cond, levels = c('sham-sham', 'cTBS-sham', 'sham-cTBS')),
    Cond_day1 = factor(Cond_day1),
    Cond_day2 = factor(Cond_day2),
    IfDevalue = factor(IfDevalue, levels = c('non-sated', 'sated'))
  ) %>%
  # Sort
  arrange(SubID, Cond, IfDevalue)

# 🔢 Add pair index (every 2 rows get the same pair ID)
dat <- dat %>%
  mutate(pair = rep(seq_len(nrow(dat) / 2), each = 2))

# save the whole rating data for later use
Odor_ratings_dat = dat
save(Odor_ratings_dat,file = '../../data_beh_processed/Odor_ratings_dat.RData')

