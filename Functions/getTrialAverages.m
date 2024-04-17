function [across_mean,across_SE] = getTrialAverages(NumTrials,across_mean,across_SE,var)
%GET TRIAL AVERAGES: Get averages and SEs across catch trials

% Calculate group averages on each catch trial and SE
for TrialNum = 1:NumTrials
    % Mean
    across_mean(1,TrialNum) = nanmean(var(:,TrialNum));
    
    % Standard error
    across_SE(1,TrialNum) = nanstd(var(:,TrialNum))/sqrt(sum(~isnan(var(:,TrialNum))));
end

% Remove nan values
across_mean(isnan(across_mean)) = [];

across_SE(isnan(across_SE)) = [];
end

