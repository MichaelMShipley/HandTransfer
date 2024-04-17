% Analyze hand transfer experiment

%% Hand transfer (3_9) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Add function and data path
addpath([cd,'\Functions'])
addpath([cd,'\Data'])

NumSubjects = 14;
NumCatchTrials = 10;
NumTrials = 60;

subjectIDs = {'eeo18e','zcuwa9','mfcpac','h4pd5a','ndvqv3','fe3ovs','swg4bs','f3j3y7','yremkm','gufrfy','k0r8cx','sipxpy','i91346','g8964i'};
baseline_weight = 3;
catch_weight = 9;
half_baseline_weight = 1.5;

% Import scored data table
fullTable = readtable('Hand_Transfer_Scored_final.txt');


%% Calculate measures from scored values and raw data: 3_9

% Initialize variables
VFR_at_half_3 = nan(NumSubjects,NumTrials);         % initial VFR at half baseline VF
VFR_at_half_9 = nan(NumSubjects,NumTrials);

liftHeight_3 = nan(NumSubjects,NumTrials);          % lift height
liftHeight_9 = nan(NumSubjects,NumTrials);

corr_VF_slope_9 = nan(NumSubjects,NumTrials);         % slope from corr_onset to lift off

correctDelay_9 = nan(NumSubjects,NumTrials);          % time from expected lift off to corrective response onset

trial_hand = nan(NumSubjects,NumTrials);

% Initialize mean arrays
mean_VFR_at_half_3 = nan(NumSubjects,4);    % columns = 1st hand, 2nd hand, Dom hand, Non-dom hand
mean_VFR_at_half_9 = nan(NumSubjects,4);    % columns = 1st hand, 2nd hand, Dom hand, Non-dom hand

mean_liftHeight_3 = nan(NumSubjects,4);     % columns = 1st hand, 2nd hand, Dom hand, Non-dom hand
mean_liftHeight_9 = nan(NumSubjects,4);     % columns = 1st hand, 2nd hand, Dom hand, Non-dom hand

mean_corr_VF_slope_9 = nan(NumSubjects,1);

mean_correctDelay_9 = nan(NumSubjects,4);   % columns = 1st hand, 2nd hand, Dom hand, Non-dom hand

for SubjectNum = 1:NumSubjects
    
    close all
    
    thisSubject = subjectIDs{SubjectNum};
        
    % Pull scored data times and values
    subj = fullTable.subj(:);
    subjectTable = fullTable(ismember(subj,thisSubject),:);         % get subject table
       
    trial = subjectTable.trial(:);                                % get trial column
    outlier = subjectTable.outlier(:);
    catchTrial = subjectTable.Catch(:);                           % get catch trial column
    catchTrialNum = subjectTable.catchTrial1_20(:);
    hand = subjectTable.hand_dom_1_Non_dom_2_(:);
    lf_on_s = subjectTable.lf_on_t(:);                            % get lf_on_s column
    peak_VFR_s = subjectTable.peak_VFR_t(:);                      % get peak_VFR_s column
    lift_on_s = subjectTable.lift_on_t(:);                        % get lift_on_s column
    corr_on_s = subjectTable.corr_resp_on_t(:);                   % get corr_resp_on_t column
    
    % Import raw data mat file
    load(['S',num2str(SubjectNum),'.mat']);
    
    % Pull raw data variables
    end_time = D.trial_end_time;
    weight = D.Weight;
    time = D.Time;
    fx1 = D.FX1;
    fx2 = D.FX2;
    fz1 = D.FZ1;
    fz2 = D.FZ2;
    posZ = D.PosZ;
    
    % Loop through trials to filter and calculate vertical force from raw
    % data
    for TrialNumber = 1:length(weight)
        
        end_time_str = cell2mat(end_time(TrialNumber,:));
        end_time_num = str2double(end_time_str(1:2))*3600 + str2double(end_time_str(4:5))*60 + str2double(end_time_str(7:8));     % convert trial end time into seconds
        
        % Raw variables
        trialTime = nonzeros(time(TrialNumber,:));
        LastSample = length(trialTime);
        
        trialFX1 = fx1(TrialNumber,1:LastSample);
        trialFX2 = fx2(TrialNumber,1:LastSample);
        
        trialFZ1 = fz1(TrialNumber,1:LastSample);
        trialFZ2 = fz2(TrialNumber,1:LastSample);
        
        trialPosZ = 0.1*(posZ(TrialNumber,1:LastSample)); % divide by 10 (convert from mm to cm)
        
        % Construct fourth-order butterworth filter with cutoff frequency
        % 14 Hz (force data was sampled at 1000 Hz)
        [b,a] = butter(4,14/(1000/2));
        pad_length = 1000;
        
        trialFX1 = [trialFX1(1)*ones(pad_length,1);trialFX1(:)];                   % Pad beginning with 1000 samples of initial value
        trialFX1 = filtfilt(b,a,trialFX1);    % filter trialFX1
        trialFX1 = trialFX1(pad_length+1:end);
        
        trialFX2 = [trialFX2(1)*ones(pad_length,1);trialFX2(:)];                   % Pad beginning with 1000 samples of initial value
        trialFX2 = filtfilt(b,a,trialFX2);    % filter trialFX2
        trialFX2 = trialFX2(pad_length+1:end);
        
        trialFZ1 = [trialFZ1(1)*ones(pad_length,1);trialFZ1(:)];                   % Pad beginning with 1000 samples of initial value
        trialFZ1 = filtfilt(b,a,trialFZ1);    % filter trialFX1
        trialFZ1 = trialFZ1(pad_length+1:end);
        
        trialFZ2 = [trialFZ2(1)*ones(pad_length,1);trialFZ2(:)];                   % Pad beginning with 1000 samples of initial value
        trialFZ2 = filtfilt(b,a,trialFZ2);    % filter trialFX2
        trialFZ2 = trialFZ2(pad_length+1:end);
        
        trialPosZ = [trialPosZ(1)*ones(pad_length,1);trialPosZ(:)];                   % Pad beginning with 1000 samples of initial value
        trialPosZ = filtfilt(b,a,trialPosZ);    % filter PosZ
        trialPosZ = trialPosZ(pad_length+1:end);                                 % Remove padding
        
        verticalForce = trialFX1 + trialFX2;                                       % Calulate vertical force
        VFRate = 1000 * diff(verticalForce);                                       % Calulate vertical force rate
        gripForce = ((-1)*trialFZ1 + (-1)*trialFZ2)/2;                             % Calculate grip force
        GFRate = 1000 * diff(gripForce);                                           % Calulate grip force rate
        
        
        % Pull out variables from raw data using scored times
        if find(ismember(trial,TrialNumber)) && outlier(ismember(trial,TrialNumber))~=1     % in case trial was skipped during scoring or its an outlier
            
            
            % Get indices from scored times
            lfOnset_i = round(1+1000*lf_on_s(ismember(trial,TrialNumber)));                        % get lfOnset index (multiply by 1000 to get index (1000 Hz),+1 since time starts at 0)
            peakVFR_i = round(1+1000*peak_VFR_s(ismember(trial,TrialNumber)));                     % get peakVFR index
            corr_on_i = round(1+1000*corr_on_s(ismember(trial,TrialNumber)));                      % get corrective onset index
            scored_liftOnset_i = round(1+1000*lift_on_s(ismember(trial,TrialNumber)));                    % get liftOnset index
            
            % Get hand on this trial and store
            trial_hand(SubjectNum,TrialNumber) = hand(ismember(trial,TrialNumber));
            
            % Store duration and end of trial
            trial_dur(SubjectNum,TrialNumber) = LastSample/1000;           %#ok<SAGROW>
            trial_end_time(SubjectNum,TrialNumber) = end_time_num;         %#ok<SAGROW>
            
            
            % %%%%%%%%%%%%%%%%% VFR at half baseline weight %%%%%%%%%%%%%%%%%%%%%%%%%%            
            m = 1;
            while verticalForce(m) < half_baseline_weight
                m = m + 1;
            end
            
            VFR_at_half = VFRate(m);
            
            % %%%%%%%%%%%%%%%%% Lift Height %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set baseline posZ at 0
            trialPosZ = trialPosZ - median(trialPosZ(1:lfOnset_i));
            
            % Find local peak of posZ after liftOnset
            t = scored_liftOnset_i;                           % begining of first window at lift onset index
            while trialPosZ(t+50)> trialPosZ(t)      % find posZ reversal using windows of 50ms
                t = t+1;
            end
            
            [liftHeight,~] = max(trialPosZ(t:t+50));    % get maximum peak of tralPosZ in window
            
            if catchTrial(ismember(trial,TrialNumber)) == 1     % calculate catch trial specific measures
                
                
                % %%%%%%%%%%%%%%%%% Corrective VF Slope %%%%%%%%%%%%%%%%%%%%%%%%%%
                m = 1;
                while verticalForce(m) < catch_weight
                    m = m + 1;
                end
                lift_onset_i = m;
                
                if lift_onset_i>corr_on_i       % Only calculate slope if corrective response starts before they reach catch weight
                    if (lift_onset_i - corr_on_i) < 200             % if time between corrective onset and lift off is less than 200 ms
                        corr_slope = 1000*(verticalForce(lift_onset_i)-verticalForce(corr_on_i))./(lift_onset_i-corr_on_i);   % *1000 since denominator is ms
                    elseif (lift_onset_i - corr_on_i) > 200         % if time between corrective onset and lift off is more than 200 ms
                        corr_slope = 1000*(verticalForce(corr_on_i+200)-verticalForce(corr_on_i))./200;                       % *1000 since denominator is ms
                    end
                else
                    corr_slope = nan;
                end
                
                % % %%%%%%%%%%%%%%%%% Corrective Response Delay %%%%%%%%%%%%%%%%%%%%%%%%%%
                m = 1;
                while verticalForce(m) < baseline_weight
                    m = m + 1;
                end
                
                correctDelay = corr_on_i - m;      % time window between force reaching baseline weight and corr onset
            end
            
            % Sort and store variables
            if weight(TrialNumber) == 3
                VFR_at_half_3(SubjectNum,TrialNumber) = VFR_at_half;
                liftHeight_3(SubjectNum,TrialNumber) = liftHeight;
            elseif weight(TrialNumber) == 9
                VFR_at_half_9(SubjectNum,TrialNumber) = VFR_at_half;
                liftHeight_9(SubjectNum,TrialNumber) = liftHeight;
                corr_VF_slope_9(SubjectNum,TrialNumber) = corr_slope;
                correctDelay_9(SubjectNum,TrialNumber) = correctDelay;
            end
            
        end
    end
    
    
    % Subject means
    mean_VFR_at_half_3(SubjectNum,1) = nanmean(VFR_at_half_3(SubjectNum,1:30));         % 1st hand
    mean_VFR_at_half_3(SubjectNum,2) = nanmean(VFR_at_half_3(SubjectNum,31:end));       % 2nd hand
    mean_VFR_at_half_9(SubjectNum,1) = nanmean(VFR_at_half_9(SubjectNum,1:30));         % 1st hand
    mean_VFR_at_half_9(SubjectNum,2) = nanmean(VFR_at_half_9(SubjectNum,31:end));       % 2nd hand
    
    mean_liftHeight_3(SubjectNum,1) = nanmean(liftHeight_3(SubjectNum,1:30));       % 1st hand
    mean_liftHeight_3(SubjectNum,2) = nanmean(liftHeight_3(SubjectNum,31:end));     % 2nd hand
    mean_liftHeight_9(SubjectNum,1) = nanmean(liftHeight_9(SubjectNum,1:30));       % 1st hand
    mean_liftHeight_9(SubjectNum,2) = nanmean(liftHeight_9(SubjectNum,31:end));     % 2nd hand
    
    mean_corr_VF_slope_9(SubjectNum,1) = nanmean(corr_VF_slope_9(SubjectNum,1:end));        % All trials
            
    mean_correctDelay_9(SubjectNum,1) = nanmean(correctDelay_9(SubjectNum,1:30));       % 1st hand
    mean_correctDelay_9(SubjectNum,2) = nanmean(correctDelay_9(SubjectNum,31:end));     % 2nd hand
    
    % Means separated by Dom/Non-dom hand for stats
    mean_VFR_at_half_3(SubjectNum,3) = nanmean(VFR_at_half_3(SubjectNum,trial_hand(SubjectNum,:)==1));  % Dom hand
    mean_VFR_at_half_3(SubjectNum,4) = nanmean(VFR_at_half_3(SubjectNum,trial_hand(SubjectNum,:)==2));  % Non-dom hand
    mean_VFR_at_half_9(SubjectNum,3) = nanmean(VFR_at_half_9(SubjectNum,trial_hand(SubjectNum,:)==1));  % Dom hand
    mean_VFR_at_half_9(SubjectNum,4) = nanmean(VFR_at_half_9(SubjectNum,trial_hand(SubjectNum,:)==2));  % Non-dom hand
    
    mean_liftHeight_3(SubjectNum,3) = nanmean(liftHeight_3(SubjectNum,trial_hand(SubjectNum,:)==1));  % Dom hand
    mean_liftHeight_3(SubjectNum,4) = nanmean(liftHeight_3(SubjectNum,trial_hand(SubjectNum,:)==2));  % Non-dom hand
    mean_liftHeight_9(SubjectNum,3) = nanmean(liftHeight_9(SubjectNum,trial_hand(SubjectNum,:)==1));  % Dom hand
    mean_liftHeight_9(SubjectNum,4) = nanmean(liftHeight_9(SubjectNum,trial_hand(SubjectNum,:)==2));  % Non-dom hand
    
    mean_correctDelay_9(SubjectNum,3) = nanmean(correctDelay_9(SubjectNum,trial_hand(SubjectNum,:)==1));  % Dom hand
    mean_correctDelay_9(SubjectNum,4) = nanmean(correctDelay_9(SubjectNum,trial_hand(SubjectNum,:)==2));  % Non-dom hand
    
    % Mean effect of weight: pred VFR, and lift height
    ME_weight_VFR_at_half(SubjectNum,:) = [nanmean(VFR_at_half_3(SubjectNum,:)) nanmean(VFR_at_half_9(SubjectNum,:))];      %#ok<SAGROW>
    ME_weight_liftHeight(SubjectNum,:) = [nanmean(liftHeight_3(SubjectNum,:)) nanmean(liftHeight_9(SubjectNum,:))];         %#ok<SAGROW>
    
    % Grand mean latency
    GM_latency(SubjectNum,1) = nanmean(correctDelay_9(SubjectNum,:));   %#ok<SAGROW>
    
    % Trial and inter-trial duration
    mean_trial_dur(SubjectNum,1) = nanmean(trial_dur(SubjectNum,:));   %#ok<SAGROW>
    
    inter_dur(SubjectNum,:) = diff(trial_end_time(SubjectNum,:));      %#ok<SAGROW>
    mean_inter_dur(SubjectNum,1) = nanmean(inter_dur(SubjectNum,:));   %#ok<SAGROW>
end

% Subtract trial duration from end of trial intervals
mean_inter_dur = mean_inter_dur - mean_trial_dur;

cts_corr_VF_slope_9 = corr_VF_slope_9(:,~all(isnan(corr_VF_slope_9)));  % catch trial only slope variable


%% 3_9 across trials figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Corrective slope variables
Across_corr_VF_slope_9 = nan(1,NumTrials);

SE_Across_corr_VF_slope_9 = nan(1,NumTrials);

[Across_corr_VF_slope_9,SE_Across_corr_VF_slope_9] = getTrialAverages(NumTrials,Across_corr_VF_slope_9,SE_Across_corr_VF_slope_9,corr_VF_slope_9);

% For x-axis as the catch trial numbers: 10
end_1st_half_c = 5;
start_2nd_half_c = 6;                       
                         
% Across trials corrective slope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Group Catch Trial # vs corrective Slope
% Plot first half
figure(1)
plot(1:end_1st_half_c,Across_corr_VF_slope_9(1:end_1st_half_c),'Marker','.','MarkerSize',15,'Color',[0,0.4470,0.7410])
hold on
fill([1:end_1st_half_c,flip(1:end_1st_half_c)],[Across_corr_VF_slope_9(1:end_1st_half_c) - SE_Across_corr_VF_slope_9(1:end_1st_half_c),flip(Across_corr_VF_slope_9(1:end_1st_half_c) + SE_Across_corr_VF_slope_9(1:end_1st_half_c))],[0,0.4470,0.7410],'linestyle','none','FaceAlpha',.3,'EdgeAlpha',.3,'HandleVisibility','off');   % Error patch
% Plot second half
plot(start_2nd_half_c:NumCatchTrials,Across_corr_VF_slope_9(start_2nd_half_c:NumCatchTrials),'Marker','.','MarkerSize',15,'Color',[0.8500 0.3250 0.0980])
hold on
fill([start_2nd_half_c:NumCatchTrials,flip(start_2nd_half_c:NumCatchTrials)],[Across_corr_VF_slope_9(start_2nd_half_c:NumCatchTrials) - SE_Across_corr_VF_slope_9(start_2nd_half_c:NumCatchTrials),flip(Across_corr_VF_slope_9(start_2nd_half_c:NumCatchTrials) + SE_Across_corr_VF_slope_9(start_2nd_half_c:NumCatchTrials))],[0.8500 0.3250 0.0980],'linestyle','none','FaceAlpha',.3,'EdgeAlpha',.3,'HandleVisibility','off');   % Error patch
% Labels
xlabel('Catch Trial Number')
ylabel('Corrective Response VF Slope (N/s)')
% Resize figure
set(gcf,'PaperPositionMode','auto')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position',[50 600 550 300]);


%% Mean group summary figure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, 5th, and 6th catch trials for VF Slope. First hand/second hand for
% correctDelay, VFR at half baseline, and lift height

% Calculate group average correct delay for first/second hand
Group_mean_correctDelay_9 = nanmean(mean_correctDelay_9(:,1:2));
Group_mean_VFR_at_half_3 = nanmean(mean_VFR_at_half_3(:,1:2));
Group_mean_liftHeight_3 = nanmean(mean_liftHeight_3(:,1:2));

% Calculate SEs for first/second hand
SE_mean_correctDelay_9 = nanstd(mean_correctDelay_9(:,1:2))./sqrt(sum(~isnan(mean_correctDelay_9(:,1:2))));
SE_mean_VFR_at_half_3 = nanstd(mean_VFR_at_half_3(:,1:2))./sqrt(sum(~isnan(mean_VFR_at_half_3(:,1:2))));
SE_mean_liftHeight_3 = nanstd(mean_liftHeight_3(:,1:2))./sqrt(sum(~isnan(mean_liftHeight_3(:,1:2))));

%%%%%%%%% Plot group average bars
figure(2);
subplot(1,4,1)
tmpData = [Across_corr_VF_slope_9(1),Across_corr_VF_slope_9(5),Across_corr_VF_slope_9(6)];
tmpSE = [SE_Across_corr_VF_slope_9(1),SE_Across_corr_VF_slope_9(5),SE_Across_corr_VF_slope_9(6)];
a = bar(tmpData);
a.FaceColor = 'flat';
a.EdgeColor = 'flat';
a.CData(1:2,:) = [126, 176, 217; 126, 176, 217]./255;
a.CData(3,:) = [235, 152, 127]./255;
ylabel('Corrective VF Slope (N/s)')
xlabel('Catch Trial')
xticklabels({'1st','5th','6th'})
hold on
er = errorbar(tmpData(1:2),tmpSE(1:2));
er.Color = [0 0.4470 0.7410];
er.LineStyle = 'none';
er = errorbar(3,tmpData(3),tmpSE(3));
er.Color = [0.8500 0.3250 0.0980];
er.LineStyle = 'none';

subplot(1,4,2)
tmpData = Group_mean_correctDelay_9;
tmpSE = SE_mean_correctDelay_9;
a = bar(tmpData);
a.FaceColor = 'flat';
a.EdgeColor = 'flat';
a.CData(1,:) = [126, 176, 217]./255;
a.CData(2,:) = [235, 152, 127]./255;
ylabel('Reponse Latency (ms)')
xlabel('Object Weight')
xticklabels({'9N','9N'})
hold on
er = errorbar(tmpData(1),tmpSE(1));
er.Color = [0 0.4470 0.7410];
er.LineStyle = 'none';
er = errorbar(2,tmpData(2),tmpSE(2));
er.Color = [0.8500 0.3250 0.0980];
er.LineStyle = 'none';

subplot(1,4,3)
tmpData = [Group_mean_VFR_at_half_3(1),Group_mean_VFR_at_half_3(2)];
tmpSE = [SE_mean_VFR_at_half_3(1),SE_mean_VFR_at_half_3(2)];
a = bar(tmpData);
a.FaceColor = 'flat';
a.EdgeColor = 'flat';
a.CData(1,:) = [126, 176, 217]./255;
a.CData(2,:) = [235, 152, 127]./255;
xlabel('Object Weight')
ylabel('Predictive VF Rate (N/s)')
ylim([0 25])
xticklabels({'3N','3N'})
hold on
er = errorbar(tmpData(1),tmpSE(1));
er.Color = [0 0.4470 0.7410];
er.LineStyle = 'none';
er = errorbar(2,tmpData(2),tmpSE(2));
er.Color = [0.8500 0.3250 0.0980];
er.LineStyle = 'none';

subplot(1,4,4)
tmpData = [Group_mean_liftHeight_3(1),Group_mean_liftHeight_3(2)];
tmpSE = [SE_mean_liftHeight_3(1),SE_mean_liftHeight_3(2)];
a = bar(tmpData,'FaceColor',[0.8,0.8,0.8]);
a.FaceColor = 'flat';
a.EdgeColor = 'flat';
a.CData(1,:) = [126, 176, 217]./255;
a.CData(2,:) = [235, 152, 127]./255;
xlabel('Object Weight')
ylabel('Lift Height (cm)')
ylim([0 4])
xticklabels({'3N','3N'})
hold on
er = errorbar(tmpData(1),tmpSE(1));
er.Color = [0 0.4470 0.7410];
er.LineStyle = 'none';
er = errorbar(2,tmpData(2),tmpSE(2));
er.Color = [0.8500 0.3250 0.0980];
er.LineStyle = 'none';


%%%%%%%%% Individual participant lines and points
for SubjectNum = 1:NumSubjects
    subplot(1,4,1)
    plot([cts_corr_VF_slope_9(SubjectNum,1),cts_corr_VF_slope_9(SubjectNum,5),cts_corr_VF_slope_9(SubjectNum,6)],'Marker','none','Color',[0.4,0.4,0.4]);
    
    subplot(1,4,2)
    plot(mean_correctDelay_9(SubjectNum,1:2),'Marker','none','Color',[0.4,0.4,0.4]);
    
    subplot(1,4,3)
    plot([mean_VFR_at_half_3(SubjectNum,1),mean_VFR_at_half_3(SubjectNum,2)],'Marker','none','Color',[0.4,0.4,0.4]);
    
    subplot(1,4,4)
    plot([mean_liftHeight_3(SubjectNum,1),mean_liftHeight_3(SubjectNum,2)],'Marker','none','Color',[0.4,0.4,0.4]);
end

set(gcf,'PaperPositionMode','auto')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position',[50 100 1200 400]);



