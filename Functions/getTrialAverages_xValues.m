function [across_mean,across_SE,x_values] = getTrialAverages_xValues(NumTrials,across_mean,across_SE,var)
%GET TRIAL AVERAGES: Get averages and SEs across catch trials

% Calculate group average on each catch trial, SE, and X values for plot
for TrialNum = 1:NumTrials
    % Mean
    across_mean(1,TrialNum) = nanmean(var(:,TrialNum));
    
    % Standard error
    across_SE(1,TrialNum) = nanstd(var(:,TrialNum))/sqrt(sum(~isnan(var(:,TrialNum))));
end

% Get x-values
x_values = find(~isnan(across_mean));

% Remove nan values
across_mean(isnan(across_mean)) = [];
across_SE(isnan(across_SE)) = [];

end
