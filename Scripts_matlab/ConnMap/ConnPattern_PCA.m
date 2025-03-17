
%% reorganize conn maps by ROIs and do MDS

for subj = 13 %:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    
    % coding sham and cTBS as 1 and 2, respectively
    % null as 3
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [3,2,1,1,2,1,1];            
        case 132
            TMS_types = [3,2,1,1,1,1,2];       
        case 213
            TMS_types = [3,1,2,2,1,1,1];            
        case 231
            TMS_types = [3,1,2,1,1,2,1];            
        case 312
            TMS_types = [3,1,1,2,1,1,2];               
        case 321
            TMS_types = [3,1,1,1,2,2,1];
    end
    
    session_order = 1:7;
    
    all_avg_corr_map = nan(n_gm_vox,n_rest_names,4);
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
        load(mat_name)
        for k = 1:4 % 4 ROIs
            tmp_loc = all_idx(:,3)==k;
            tmp_avg_corr_map = mean(corr_map(tmp_loc,:),1);
            all_avg_corr_map(:,r,k) = tmp_avg_corr_map';
        end
    end
    
    % handle missing sessions
    % (be careful not to rerun this)
    n_use_sessions = sum(count_table(subj,:)==1);
    missing_loc = count_table(subj,:)<1;
    session_order(missing_loc) = [];
    TMS_types(missing_loc) = [];
    all_avg_corr_map(:,missing_loc,:) = [];
    
    % conduct MDS for each ROI
    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 700]); % plot without showing
    for k = 1:4
        
        % it is suboptimal to do pca with 7 observations and 58660 columns
        % better to do it without averaging voxels within ROIs to get more
        % components
        
        connectivity_matrix = all_avg_corr_map(:,:,k);
        connectivity_matrix = connectivity_matrix';   % transpose so sessions as rows (obs), voxels as columns (variables)
        
        distances = pdist(connectivity_matrix, 'euclidean');  % pair-wise distance
        distances_square = squareform(distances);  % Convert the distances into a square form
        opts = optimset('MaxIter', 1e5);
        [Y, stress] = mdscale(distances_square, 2, 'Options', opts); % 2D MDS
    
        subplot(2,2,k)
        hold on;
        colors = lines(3); 
        for i = 1:n_use_sessions
            tmp_size = session_order(i)*25;
            if TMS_types(i) == 1
                scatter(Y(i, 1), Y(i, 2), tmp_size, 'filled', 'MarkerFaceColor', colors(1, :));
            elseif TMS_types(i) == 2
                scatter(Y(i, 1), Y(i, 2), tmp_size, 'filled', 'MarkerFaceColor', colors(2, :));
            else
                scatter(Y(i, 1), Y(i, 2), tmp_size, 'filled', 'MarkerFaceColor', colors(3, :));
            end
            text(Y(i, 1), Y(i, 2), num2str(session_order(i)), 'FontSize', 18, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end

        xlabel('Dimension 1');
        ylabel('Dimension 2');
        title(sprintf('MDS of Conn. Patterns (%s)',labels{k}));
        subtitle(sprintf('Num. PC=%d, Stress=%f',num_components,stress))
        hold off;
    end
    
    img_name = fullfile(func_conn_dir,'MDS_conn_patterns_sessions.bmp');
    saveas(gcf, img_name);
    
end
