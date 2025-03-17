
clear; clc; close all

medir = '/Users/liuq13/Documents/Projects/NODEAP_data_analysis/ME_example';

curr_rest = 'D0';

% sliceViewer(Y) % to display images of multiple slices for multiple echoes

% Display multiple slices
sliceIndex = [4,8,16,24,32,40];
numSlices = length(sliceIndex); % Number of slices to display

figure;
rowLabels = {'Echo 1','Echo 2','Echo 3','Combined'};

tiledlayout(4, numSlices, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:4

if k==4
    file = fullfile(medir,'fvol_4d.nii');
else
    n = dir(fullfile(medir, sprintf('%s_rest*_e%d.nii',curr_rest,k)));
    file = fullfile(n.folder,n.name);
end

V = spm_vol(file);  
Y = spm_read_vols(V(1));

for i = 1:numSlices
    nexttile; % Move to the next tile in the layout
    %subplot(4, numSlices, i+numSlices*(k-1)); % Adjust the subplot grid as needed
    imagesc(Y(:, :, sliceIndex(i)));
    colormap gray;
    axis off; % Turn off all axes
    if k==1
        title(['Slice ', num2str(sliceIndex(i))], 'FontSize', 18);
    end
end

end

set(gcf, 'Position', [100, 100, 1200, 800]); 
saveas(gcf, fullfile(medir,'mri_slices.bmp'));

%%
% this part looks at tedana output

tiledlayout(4, numSlices, 'TileSpacing', 'compact', 'Padding', 'compact');

file = fullfile(medir, 'T2starmap.nii.gz');
V = spm_vol(file);  
Y = spm_read_vols(V);
for i = 1:numSlices
    nexttile; % Move to the next tile in the layout
    imagesc(Y(:, :, sliceIndex(i)));
    colormap gray;
    axis off; % Turn off all axes
    title(['Slice ', num2str(sliceIndex(i))], 'FontSize', 18);
end

file = fullfile(medir, 'S0map.nii.gz');
V = spm_vol(file);  
Y = spm_read_vols(V);
for i = 1:numSlices
    nexttile; % Move to the next tile in the layout
    imagesc(Y(:, :, sliceIndex(i)));
    colormap gray;
    axis off; % Turn off all axes
end
   
file = fullfile(medir, 'desc-optcom_bold.nii.gz');
V = spm_vol(file);  
Y = spm_read_vols(V);
for i = 1:numSlices
    nexttile; % Move to the next tile in the layout
    imagesc(Y(:, :, sliceIndex(i)));
    colormap gray;
    axis off; % Turn off all axes
end

file = fullfile(medir, 'desc-optcomDenoised_bold.nii.gz');
V = spm_vol(file);  
Y = spm_read_vols(V);
for i = 1:numSlices
    nexttile; % Move to the next tile in the layout
    imagesc(Y(:, :, sliceIndex(i)));
    colormap gray;
    axis off; % Turn off all axes
end

set(gcf, 'Position', [100, 100, 1200, 800]); 
saveas(gcf, fullfile(medir,'tedana_slices.bmp'));
    


