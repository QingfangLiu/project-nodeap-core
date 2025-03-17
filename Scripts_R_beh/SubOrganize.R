

rm(list = ls())
source('Setup.R')
HomeDir = '/Users/liuq13/NODEAP'
sub_dat = read.xlsx(file.path(HomeDir,'NODEAP_DataCollectionSheet_QL.xlsx'))

sub_dat = sub_dat %>%
  select(Subject,Gender,Birth.Year,`Anterior/Posterior`,
         Order,Odor.Selection,`Sweet/Savory`)
colnames(sub_dat) = c('SubID','Gender','Birth.Year','StimLoc','StimOrder','Odors','StartOdor')

# clean those spaces
sub_dat$SubID = str_replace_all(sub_dat$SubID," ","")
sub_dat$Gender = str_replace_all(sub_dat$Gender," ","")
sub_dat$StimLoc = str_replace_all(sub_dat$StimLoc," ","")
sub_dat$StartOdor = str_replace_all(sub_dat$StartOdor," ","")

sub_dat_demo = readxl::read_excel(file.path(HomeDir,'Data Organization_QL.xlsx'))
sub_dat$Age = sub_dat_demo$Age
sub_dat$Sex = sub_dat_demo$Sex
sub_dat$Include = 1 # whether or not including this subject

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

write.xlsx(sub_dat,file = '../ProcessedData/SubConds_raw.xlsx')

table(sub_dat$Sex)
summary(sub_dat$Age)
sd(sub_dat$Age)






