% Plot individual hand transfer traces
clear
close all

% Add data path
addpath([cd,'\Data'])

mode = 2;   % 1: plot all individual traces, 2: plot example trials
% examples traces: participant 3 (mfcpac) trials 41(catch) and 43(base)

subjectIDs = {'eeo18e','zcuwa9','mfcpac','h4pd5a','ndvqv3','fe3ovs','swg4bs','f3j3y7','yremkm','gufrfy','k0r8cx','sipxpy','i91346','g8964i'};

% Import scored data table
fullTable = readtable('Hand_Transfer_Scored_final.txt');

NumSubjects = 14;
FirstTrial = 1;
LastTrial = 60;
plot_length = 6000;

green = [81, 191, 59]/255;                    % Plot colors
blue = [15, 105, 189]/255;
red = [212, 22, 19]/255;

% Plot individual traces
for participant = 1:NumSubjects
    this_ID = subjectIDs{participant};
    
    subj = fullTable.subj(:);
    subjectTable = fullTable(ismember(subj,this_ID),:);         % get subject table
    
    trial = subjectTable.trial(:);                              % get trial column
    catchTrial = subjectTable.Catch(:);                         % get catch trial column
    catchTrialNum = subjectTable.catchTrial1_20(:);             % get catch trial number
    corr_on_t = subjectTable.corr_resp_on_t(:);                 % get corr_resp_on_t column
    outlier = subjectTable.outlier(:);                          % outlier column
    
    load(['S',num2str(participant),'.mat']);      % Import raw data mat file
    
    for TrialNumber = FirstTrial:LastTrial
        
        % Construct fourth-order low pass butterworth filter with cutoff frequency
        % 14 Hz (force data was sampled at 1000 Hz)
        [b,a] = butter(4,14/(1000/2));
        pad_length = 1000;
        
        Tmp = D.Time(TrialNumber,:);
        Time = Tmp(:);
        Time = nonzeros(Time);      % trim zeros
        LastSample = length(Time);
        
        Time = 0:length(Time);  % create standardized time array in ms
        
        Weight = D.Weight(TrialNumber);
        
        Tmp = D.FX1(TrialNumber,1:LastSample);
        FX1 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FX1 = filtfilt(b,a,FX1);    % filter FX1
        FX1 = FX1(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.FY1(TrialNumber,1:LastSample);
        FY1 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FY1 = filtfilt(b,a,FY1);    % filter FY1
        FY1 = FY1(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.FZ1(TrialNumber,1:LastSample);
        FZ1 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FZ1 = filtfilt(b,a,FZ1);    % filter FZ1
        FZ1 = FZ1(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.FX2(TrialNumber,1:LastSample);
        FX2 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FX2 = filtfilt(b,a,FX2);    % filter FX2
        FX2 = FX2(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.FY2(TrialNumber,1:LastSample);
        FY2 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FY2 = filtfilt(b,a,FY2);    % filter FY2
        FY2 = FY2(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.FZ2(TrialNumber,1:LastSample);
        FZ2 = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        FZ2 = filtfilt(b,a,FZ2);    % filter FZ2
        FZ2 = FZ2(pad_length+1:end);                                  % Remove padding
        
        Tmp = D.PosZ(TrialNumber,1:LastSample);
        PosZ = [Tmp(1)*ones(pad_length,1);Tmp(:)];                   % Pad beginning with 1000 samples of initial value
        PosZ = filtfilt(b,a,PosZ);    % filter PosZ
        PosZ = PosZ(pad_length+1:end);                                  % Remove padding
        PosZ = (PosZ - PosZ(10))/10;                     % Transform PosZ and convert from mm to cm + subtract baseline height
        
        VF = FX1 + FX2;                                              % Calulate vertical force
        VFRate = diff(VF);                                           % Derivative of vertical force
        VFRate = VFRate*1000;                                        % Convert from N/ms to N/s (*1000)
        
        VFRate_rate = diff(VFRate)*1000;
        
        GF = ((-1)*FZ1 + (-1)*FZ2)/2;                                % Calculate grip force
        GFRate = diff(GF);                                           % Derivative of grip force
        GFRate = GFRate*1000;
        
        VelZ = diff(PosZ)*1000;                                      % Calculate velocity
        
        
        % Line up start of VF/VFRate traces by the scored corrective
        % response onset
        if mode == 1
            if find(ismember(trial,TrialNumber)) && outlier(ismember(trial,TrialNumber))~=1 && catchTrial(ismember(trial,TrialNumber)) == 1   % if trial is not missing and is not an outlier and is a catch trial
                onset_time = corr_on_t(ismember(trial,TrialNumber));             % get onset_time
                onset_time = 1000*onset_time;
                
                [startMin,Start1] = min(abs(Time - onset_time));
                
                timeArray = (Start1-1)*-1:plot_length;
                Start = 1;
                
                % Plot VF traces on 1st, 5th and 6th catch trials colour coded by pre/post transfer              
                if Weight == 9
                    if catchTrialNum(ismember(trial,TrialNumber)) == 1
                        % VF
                        subplot(1,3,1)
                        hold on
                        plot(timeArray(1:length(VF(Start:end))),VF(Start:end),'color',blue,'LineWidth',0.25,'HandleVisibility','off')
                        
                    elseif catchTrialNum(ismember(trial,TrialNumber)) == 5
                        % VF
                        subplot(1,3,2)
                        hold on
                        plot(timeArray(1:length(VF(Start:end))),VF(Start:end),'color',blue,'LineWidth',0.25,'HandleVisibility','off')
                      
                    elseif catchTrialNum(ismember(trial,TrialNumber)) == 6
                        % VF
                        subplot(1,3,3)
                        hold on
                        plot(timeArray(1:length(VF(Start:end))),VF(Start:end),'color',red,'LineWidth',0.25,'HandleVisibility','off')
                      
                    end                    
                end
            end
        elseif mode == 2 && participant == 3 && (TrialNumber == 41 || TrialNumber == 43)
            if TrialNumber == 41
                color_tmp = red;
                
                onset_time = corr_on_t(ismember(trial,TrialNumber));             % get onset_time
                onset_time = 1000*onset_time;
                
                [startMin,Start2] = min(abs(Time - onset_time));
                
                timeArray = (Start2-1)*-1:plot_length;
                Start = 1;
            elseif TrialNumber == 43
                color_tmp = green;
                
                timeArray = (Start2+110)*-1:plot_length;
                Start = 1;
            end
            
            % VF
            subplot(5,1,1)
            hold on
            plot(timeArray(1:length(VF(Start:end))),VF(Start:end),'color',color_tmp,'LineWidth',0.25,'HandleVisibility','off')
            
            % VFRate
            subplot(5,1,2)
            hold on
            plot(timeArray(1:length(VFRate(Start:end))),VFRate(Start:end),'color',color_tmp,'LineWidth',0.25,'HandleVisibility','off')
            
            % VFRate_rate
            subplot(5,1,3)
            hold on
            plot(timeArray(1:length(VFRate_rate(Start:end))),VFRate_rate(Start:end),'color',color_tmp,'LineWidth',0.25,'HandleVisibility','off')
            
            % Position
            subplot(5,1,4)
            hold on
            plot(timeArray(1:length(PosZ(Start:end))),PosZ(Start:end),'color',color_tmp,'LineWidth',0.25,'HandleVisibility','off')
            
            % Velocity
            subplot(5,1,5)
            hold on
            plot(timeArray(1:length(VelZ(Start:end))),VelZ(Start:end),'color',color_tmp,'LineWidth',0.25,'HandleVisibility','off')
        end
    end % for TrialNumber = FirstTrial:LastTrial
end


%% Figure sizes, titles, axis labels etc.
if mode == 1
    figure(1)
    subplot(1,3,1)
    plot([-400,1000],[0,0],'color','k','LineStyle','--')        % x-line  
    plot([0,0],[-1,14],'color','k','LineStyle','--')            % y-line
    patch([0,200,200,0],[-1,-1,14,14],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    title('Catch 1')
    xlim([-500,600])
    xticks([-400 -200 0 200 400 600 800])
    xticklabels({'-400','-200','0','200','400','600','800'})
    ylim([-1,14])
    yticks([0 4 8 12])
    yticklabels({'0','4','8','12'})
    ylabel('Vertical Force (N)')
    
    subplot(1,3,2)
    plot([-400,1000],[0,0],'color','k','LineStyle','--')        % x-line  
    plot([0,0],[-1,14],'color','k','LineStyle','--')            % y-line
    patch([0,200,200,0],[-1,-1,14,14],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    title('Catch 5')
    xlim([-500,600])
    xticks([-400 -200 0 200 400 600 800])
    xticklabels({'-400','-200','0','200','400','600','800'})
    xlabel('Time re: Start of Corrective Response (ms)')
    ylim([-1,14])
    yticks([0 4 8 12])
    yticklabels({'0','4','8','12'})
    
    subplot(1,3,3)
    plot([-400,1000],[0,0],'color','k','LineStyle','--')        % x-line  
    plot([0,0],[-1,14],'color','k','LineStyle','--')            % y-line
    patch([0,200,200,0],[-1,-1,14,14],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    title('Catch 6')
    xlim([-500,600])
    xticks([-400 -200 0 200 400 600 800])
    xticklabels({'-400','-200','0','200','400','600','800'})
    ylim([-1,14])
    yticks([0 4 8 12])
    yticklabels({'0','4','8','12'})
    
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[50 75 1200 500]);
    
elseif mode == 2
    figure(1)
    subplot(5,1,1)
    plot([0,0],[-1,14],'color','k','LineStyle','--')
    plot([200,200],[-1,14],'color','k','LineStyle','--')
    patch([0,200,200,0],[-1,-1,14,14],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    xlim([-400,800])
    xticks([-400 100])
    xticklabels({'0','500'})
    yticks([0 4])
    yticklabels({'0','4'})
    ylabel('Vertical Force (N)')
    
    subplot(5,1,2)
    plot([0,0],[-20,30],'color','k','LineStyle','--')
    plot([200,200],[-20,30],'color','k','LineStyle','--')
    patch([0,200,200,0],[-20,-20,30,30],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    xlim([-400,800])
    xticks([-400 100])
    xticklabels({'0','500'})
    yticks([0 20])
    yticklabels({'0','20'})
    ylabel('Vertical Force Rate (N/s)')
    
    subplot(5,1,3)
    plot([0,0],[-400,800],'color','k','LineStyle','--')
    plot([200,200],[-400,800],'color','k','LineStyle','--')
    patch([0,200,200,0],[-400,-400,800,800],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    xlim([-400,800])
    xticks([-400 100])
    xticklabels({'0','500'})
    yticks([0 400])
    yticklabels({'0','400'})
    ylabel('Vertical Force Rate Rate (N/s^2)')
    
    subplot(5,1,4)
    plot([0,0],[-1,3],'color','k','LineStyle','--')
    plot([200,200],[-1,3],'color','k','LineStyle','--')
    patch([0,200,200,0],[-1,-1,3,3],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    xlim([-400,800])
    xticks([-400 100])
    xticklabels({'0','500'})
    yticks([0 2])
    yticklabels({'0','2'})
    ylabel('Vertical Position (cm)')
    
    subplot(5,1,5)
    plot([0,0],[-5,20],'color','k','LineStyle','--')
    plot([200,200],[-5,20],'color','k','LineStyle','--')
    patch([0,200,200,0],[-5,-5,20,20],[.83,.83,.83],'FaceColor',[.83,.83,.83],'EdgeColor',[.83,.83,.83],'FaceAlpha',.3,'EdgeAlpha',.3)      % shaded region
    xlim([-400,800])
    xticks([-400 100])
    xticklabels({'0','500'})
    xlabel('Time (ms)')
    yticks([0 10])
    yticklabels({'0','10'})
    ylabel('Vertical Velocity (cm/s)')
    
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[50 75 300 800]);
end

