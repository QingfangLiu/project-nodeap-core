
############################################################
# NODEAP: Subject Conditions, Calories, and TMS Ratings Prep
#
# Summary (high level):
# - Load subject metadata and session logs from Excel.
# - Clean/harmonize IDs and labels; add Age/Sex from a demo sheet.
# - Convert Excel dates to Date and derive within-/between-session gaps.
# - Apply a manual correction: NODEAP_73 StimLoc → Posterior.
# - Save the cleaned subject condition table:
#     beh_data_found/SubConds.xlsx
# - Read calories for each Day-2 session (S1D2/S2D2/S3D2),
#   compute basic summaries, and align them with subjects.
# - Map stimulation order codes (e.g., "123") to labels
#   (e.g., "CS-SC-SS") and cross-tab by Sex/StimLoc/Order.
# - For all six sessions (S1D1–S3D2), load TMS ratings
#   ("Uncomfortable", "Strong"); derive per-session TMS type
#   (cTBS vs sham) from the mapped order.
# - Build a long data frame (subject × session) and run
#   mixed-effects models testing effects of TMS type and, within cTBS,
#   StimLoc (aOFC vs pOFC).
# - Z-score ratings within subject, relabel factors, and save:
#     beh_data_found/Survey_uncomf_strong.xlsx
############################################################


rm(list = ls())

# load_setup.R
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))
dat_folder = '/Volumes/X9Pro/NODEAP'

# 1) Load & clean subject metadata --------------------------------------
#   - read NODEAP_DataCollectionSheet_QL.xlsx (main sheet)
#   - select/rename, strip spaces, add age/sex from Data Organization_QL.xlsx

sub_dat = read.xlsx(file.path(dat_folder,'experiment_metadata','NODEAP_DataCollectionSheet_QL.xlsx')) %>%
  select(Subject,`Anterior/Posterior`,
         Order,Odor.Selection,`Sweet/Savory`)
colnames(sub_dat) = c('SubID','StimLoc','StimOrder','Odors','StartOdor')

# clean those spaces
sub_dat$SubID = str_replace_all(sub_dat$SubID," ","")
sub_dat$StimLoc = str_replace_all(sub_dat$StimLoc," ","")
sub_dat$StartOdor = str_replace_all(sub_dat$StartOdor," ","")

sub_dat_demo = readxl::read_excel(file.path(dat_folder,'experiment_metadata','Data Organization_QL.xlsx'))
sub_dat$Age = round(sub_dat_demo$Age)
sub_dat$Sex = sub_dat_demo$Sex

# 2) Dates → intervals ---------------------------------------------------
#   - convert SxDy to Date
#   - compute btS1, btS2, btS3, btS1S2, btS2S3
#   - quick sanity stats/plots

S1D1 = as.Date(sub_dat_demo$S1D1, origin = "1899-12-30")
S1D2 = openxlsx::convertToDate(sub_dat_demo$S1D2)
S1D2[46] = as_date("2023-04-27")
S2D1 = as.Date(sub_dat_demo$S2D1, origin = "1899-12-30")
S2D2 = as.Date(sub_dat_demo$S2D2, origin = "1899-12-30")
S3D1 = as.Date(sub_dat_demo$S3D1, origin = "1899-12-30")
S3D2 = as.Date(sub_dat_demo$S3D2, origin = "1899-12-30")

sub_dat$btS1 = difftime(S1D2,S1D1,units = 'days')
sub_dat$btS2 = difftime(S2D2,S2D1,units = 'days')
sub_dat$btS3 = difftime(S3D2,S3D1,units = 'days')

sub_dat$btS1S2 = difftime(S2D1,S1D1,units = 'days') # b/t Sess 1 and 2
sub_dat$btS2S3 = difftime(S3D1,S2D1,units = 'days') # b/t Sess 2 and 3

gap = as.numeric(c(sub_dat$btS1S2,
                   sub_dat$btS2S3))
summary(gap) 
hist(gap)
sd(gap)

table(sub_dat$Sex)
summary(sub_dat$Age)
sd(sub_dat$Age)

# 3) Manual fixes & save SubConds ---------------------------------------
#   - NODEAP_73 StimLoc fix
#   - write beh_data_found/SubConds.xlsx

# --- Fix StimLoc for NODEAP_73 ---
# This participant was originally assigned to the Anterior group,
# but TMS was actually applied to Posterior-associated coordinates.
# Update the StimLoc value to 'Posterior' for accurate downstream analysis.

sub_dat$StimLoc[sub_dat$SubID == 'NODEAP_73'] <- 'Posterior'

# Save the corrected dataset
outfile <- file.path(project_folder, 'beh_data_found', 'SubConds.xlsx')
write.xlsx(sub_dat, file = outfile)

# 4) Calories (D2 only) --------------------------------------------------
#   - read S1D2/S2D2/S3D2 sheets → Calories matrix (n_sub × 3)
#   - clean known NA; summary stats

## for calculating calories
# each subject, each D2 session
sess_names = c('S1D2','S2D2','S3D2')
sub_dat = NULL
sess_ctr = 0
Calories = matrix(NA,48,length(sess_names))
for(sess in sess_names){
  sess_ctr = sess_ctr + 1
  tmp_dat = read.xlsx(file.path(dat_folder,'experiment_metadata','NODEAP_DataCollectionSheet_QL.xlsx'),sheet = sess)
  Calories[,sess_ctr] = as.numeric(tmp_dat$Calories)
  
}

# 5) StimOrder mapping & cross-tabs -------------------------------------
#   - read SubConds.xlsx and map 123→CS-SC-SS (etc.)
#   - tables by StimLoc/Sex/StimOrder

# put the calories with the sub cond together
sub_cond_dat = read.xlsx(file.path(project_folder, 'beh_data_found','SubConds.xlsx')) %>%
  mutate(StimOrder=mapvalues(StimOrder,from=c('123','132','213','231','312','321'),
                             to=c('CS-SC-SS','CS-SS-SC','SC-CS-SS',
                                  'SC-SS-CS','SS-CS-SC','SS-SC-CS')))

cbind(sub_cond_dat$SubID,Calories)

Calories[10,2] = NA # remove this session not included in the beh analysis

mean(Calories,na.rm=T)
sd(Calories,na.rm=T)/sqrt(48)

# overall of assignment
table(sub_cond_dat$StimLoc) # 23 aOFC, 25 pOFC
table(sub_cond_dat$Sex) # 32 females, 16 males
table(sub_cond_dat$Sex,sub_cond_dat$StimLoc) 

table(sub_cond_dat$StimOrder)
table(sub_cond_dat$Sex,sub_cond_dat$StimOrder)

# 6) TMS ratings (6 sessions) -------------------------------------------
#   - read Uncomfortable/Strong from all SxDy sheets → matrices
#   - derive TMS_types from mapped StimOrder (C/S sequence per session)
#   - build long DF (one row per subject×session)

# for each session
sess_names = c('S1D1','S1D2','S2D1','S2D2','S3D1','S3D2')
sub_dat = NULL
sess_ctr = 0
uncomfortable = matrix(NA,48,length(sess_names))
strong = matrix(NA,48,length(sess_names))
for(sess in sess_names){
  sess_ctr = sess_ctr + 1
  tmp_dat = read.xlsx(file.path(dat_folder,'experiment_metadata','NODEAP_DataCollectionSheet_QL.xlsx'),sheet = sess)
  uncomfortable[,sess_ctr] = as.numeric(tmp_dat$`Uncomfortable.TMS.(0-10)`)
  strong[,sess_ctr] = as.numeric(tmp_dat$`Strong.TMS.(0-10)`)
}

TMS_types <- vector("list", nrow(sub_cond_dat))
for (j in 1:nrow(sub_cond_dat)) {
  tmp_order <- sub_cond_dat$StimOrder[j]
  
  if (tmp_order == 'CS-SC-SS') TMS_types[[j]] <- c('C', 'S', 'S', 'C', 'S', 'S')
  if (tmp_order == 'CS-SS-SC') TMS_types[[j]] <- c('C', 'S', 'S', 'S', 'S', 'C')
  if (tmp_order == 'SC-CS-SS') TMS_types[[j]] <- c('S', 'C', 'C', 'S', 'S', 'S')
  if (tmp_order == 'SC-SS-CS') TMS_types[[j]] <- c('S', 'C', 'S', 'S', 'C', 'S')
  if (tmp_order == 'SS-CS-SC') TMS_types[[j]] <- c('S', 'S', 'C', 'S', 'S', 'C')
  if (tmp_order == 'SS-SC-CS') TMS_types[[j]] <- c('S', 'S', 'S', 'C', 'C', 'S')
}


TMS_rating_df = data.frame(SubID=rep(sub_cond_dat$SubID,each=6),
                           StimLoc = rep(sub_cond_dat$StimLoc,each=6),
                           sess_name = rep(sess_names,48),
                           TMSSess = rep(1:6,48), # TMS Sess from 1 to 6
                           TMS_type = unlist(TMS_types),
                           uncomfortable=as.numeric(t(uncomfortable)),
                           strong=as.numeric(t(strong)))
TMS_rating_df = TMS_rating_df %>%
  mutate(SubID=factor(SubID),
         StimLoc=factor(StimLoc),
         TMS_type=factor(TMS_type))

# 7) Descriptives & mixed-effects models --------------------------------
#   - group means by TMS_type (and ×StimLoc)
#   - lmer models for uncomfortable/strong; cTBS-only StimLoc test

TMS_rating_df %>%
  group_by(TMS_type) %>%
  dplyr::summarise(uncomfortable=mean(uncomfortable,na.rm = T),
                   strong=mean(strong,na.rm = T))

model_1 <- lmer(uncomfortable ~ TMS_type + (1|SubID), data = TMS_rating_df)
model_0 <- lmer(uncomfortable ~ (1|SubID), data = TMS_rating_df)
anova(model_0,model_1)
summary(model_1)

model_1 <- lmer(strong ~ TMS_type + (1|SubID), data = TMS_rating_df)
model_0 <- lmer(strong ~ (1|SubID), data = TMS_rating_df)
anova(model_0,model_1)
summary(model_1)

TMS_rating_df %>%
  group_by(TMS_type,StimLoc) %>%
  dplyr::summarise(uncomfortable=mean(uncomfortable,na.rm = T),
                   strong=mean(strong,na.rm = T))

# cTBS only, does aOFC or pOFC matter?
use_dat = subset(TMS_rating_df,TMS_type=='C')
model_1 <- lmer(uncomfortable ~ StimLoc + (1|SubID), data = use_dat)
model_0 <- lmer(uncomfortable ~ (1|SubID), data = use_dat)
anova(model_0,model_1)

model_1 <- lmer(strong ~ StimLoc + (1|SubID), data = use_dat)
model_0 <- lmer(strong ~ (1|SubID), data = use_dat)
anova(model_0,model_1)

# 8) Within-subject z-score & save --------------------------------------
#   - z-score per subject; relabel factors
#   - write beh_data_found/Survey_uncomf_strong.xlsx

# standardize those ratings (z-score)
TMS_rating_df <- TMS_rating_df %>%
  group_by(SubID) %>%
  mutate(uncomfortable = (uncomfortable - mean(uncomfortable, na.rm = TRUE)) / 
           sd(uncomfortable, na.rm = TRUE))  %>%
  mutate(strong = (strong - mean(strong, na.rm = TRUE)) / 
           sd(strong, na.rm = TRUE))
levels(TMS_rating_df$StimLoc)=c('aOFC','pOFC')
levels(TMS_rating_df$TMS_type)=c('cTBS','sham')

# people tend to rate cTBS as more uncomfortable and stronger
# than sham sessions, very clear pattern

# save TMS rating data (uncomfortable and strong)
outfile <- file.path(project_folder, 'beh_data_found', 'Survey_uncomf_strong.xlsx')
write.xlsx(TMS_rating_df, file = outfile)




