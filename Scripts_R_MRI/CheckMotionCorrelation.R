

# some quick code to correlate motion parameters from
# different phases of multi-echo preprocessing
# in an old approach

path = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_08/nii_for_TMS/day0_me'
path_func = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_08/nii_for_TMS/functional'
path_new = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_08/nii/D0_me'

path = '/Users/liuq13/NODEAP/MRI/NODEAP_06/nii/day0_me'
path_func = '/Users/liuq13/NODEAP/MRI/NODEAP_06/nii/functional'
#path_new = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_08/nii/D0_me'


txt_files <- list.files(path = path, pattern = "^rp_.*\\.txt$", full.names = TRUE, recursive = TRUE)
txt_func = list.files(path = path_func, pattern = "^rp_.*\\.txt$", full.names = TRUE, recursive = TRUE)
#txt_new = list.files(path = path_new, pattern = "^rp_.*\\.txt$", full.names = TRUE, recursive = TRUE)

rp1 = read.table(txt_files[1])
rp2 = read.table(txt_files[2])
rp3 = read.table(txt_files[3])
rp4 = read.table(txt_func)
#rp5 = read.table(txt_new)
  
df = data.frame(e1=rp1$V4,
                e2=rp2$V4,
                e3=rp3$V4)
pairs(df)
