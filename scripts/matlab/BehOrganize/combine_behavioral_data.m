
% combine_behavioral_data.m
% This script combines behavioral data files: Choices, Conditioning, and Pleasant
% from main and NODEAP_17 sources.

% Define absolute path to data_input
input_dir = '/Users/liuq13/NODEAP_scripts/data_input';

% File pairs to combine
file_pairs = {
    'Choices.xlsx', 'Choices_NODEAP_17.xlsx', 'Choices_all.xlsx';
    'Conditioning.xlsx', 'Conditioning_NODEAP_17.xlsx', 'Conditioning_all.xlsx';
    'Pleasant.xlsx', 'Pleasant_NODEAP_17.xlsx', 'Pleasant_all.xlsx'
};

% Combine each file pair
for i = 1:size(file_pairs, 1)
    file1 = fullfile(input_dir, file_pairs{i, 1});
    file2 = fullfile(input_dir, file_pairs{i, 2});
    outfile = fullfile(input_dir, file_pairs{i, 3});

    % Read both files
    data1 = readtable(file1);
    data2 = readtable(file2);

    % Combine tables
    combined_data = [data1; data2];

    % Sort by SubID (ascending)
    if any(strcmp('SubID', combined_data.Properties.VariableNames))
        combined_data = sortrows(combined_data, 'SubID');
    else
        warning('SubID column not found in %s or %s', file1, file2);
    end

    % Write output
    writetable(combined_data, outfile);

    fprintf('✅ Combined and saved: %s\n', outfile);
end
