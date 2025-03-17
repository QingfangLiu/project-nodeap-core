% Computes (unnormalized) cross-covariance function out to +/- L lags 
% (TR shifts) between each column of Avg1 and Avg2

% Avg1 and Avg2 are time x region matrices/vectors

% See tdmx_template.m for appropriate normalization of lagged_cov output

function  r = lagged_cov(Avg1,Avg2,L)

	L1 = size(Avg1,2);
	L2 = size(Avg2,2);
	r = single(zeros(L1,L2,2*L+1));

    parfor k = 1:(2*L+1)
        i = k - (L + 1);  % Calculate the lag index (from -L to L)
		tau = abs(i);
        
        % when L is positive, Avg1 contains past values, Avg2 contains
        % future values
        % when L is negative, Avg1 contains future values, Avg2 contains
        % past values
         
        if i >=0
            Avg1_lagged = Avg1(1:end-tau,:);
            Avg2_lagged = Avg2(1+tau:end,:);
        else
            Avg1_lagged = Avg1(1+tau:end,:);
            Avg2_lagged = Avg2(1:end-tau,:);
        end    
        
        % Compute the cross-covariance for this lag
		r(:,:,k) = Avg1_lagged' * Avg2_lagged;
        
        % understanding of r:
        % with positive i, r(1,2) indicates past of voxel 1 correlating 
        % with future voxel 2;
        % with negative i, r(1,2) indicates past of voxel 2 correlating
        % with future voxel 1;
        
	end

end
