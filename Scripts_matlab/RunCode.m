

home_dir = '/Users/liuq13/Documents/Projects/NODEAP_data_analysis/Scripts_matlab';
addpath(genpath(home_dir));

run(fullfile(home_dir,'MRI_data','Do_6_post_PAID.m'))

home_dir = '/Users/liuq13/Documents/Projects/NODEAP_data_analysis/Scripts_matlab';
run(fullfile(home_dir,'FuncConn_MNI','FuncConn_PAID.m'))