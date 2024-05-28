function [angle] = calculateHDfromDLC(varargin)

% This function calculates the head direction of the animals in the form of 
% an angle w.r.t. the x-direction using the tracked coordinates from 
% DeepLabCut. It also determines the direction of turn of the animals after
% the lick at the home port to start a new trial. 

addpath('I:\Videos');

p = inputParser;
addParameter(p,'root',pwd,@isstr);
addParameter(p,'dlcPath','F:\DLC_ACgN_HeadDirection-athina-2023-03-03\videos_new',@isstr);

parse(p,varargin{:});
root = p.Results.root;
dlcPath = p.Results.dlcPath;



sessions = {'IZ39\Final\IZ39_220622_sess8', 'IZ39\Final\IZ39_220629_sess12', ...
'IZ39\Final\IZ39_220705_sess16', 'IZ39\Final\IZ39_220707_sess17', ...
'IZ40\Final\IZ40_220705_sess15', 'IZ40\Final\IZ40_220707_sess16', ...
'IZ43\Final\IZ43_220826_sess2', 'IZ43\Final\IZ43_220830_sess6', ...
'IZ43\Final\IZ43_220901_sess8', 'IZ43\Final\IZ43_220911_sess9', ...
'IZ43\Final\IZ43_220915_sess13', 'IZ43\Final\IZ43_220920_sess15',...    
'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5', ...
'IZ44\Final\IZ44_220829_sess6', 'IZ44\Final\IZ44_220912_sess10', ...
'IZ44\Final\IZ44_220915_sess13', 'IZ44\Final\IZ44_220920_sess15'}; 

for s = 1:length(sessions)
    basepath = fullfile(root, sessions{s});

    destination = strcat('I:\Videos\', sessions{s});
    if ~exist(destination)
        mkdir(destination)
    end

    csvFile = dir([dlcPath filesep extractAfter(sessions{s},'Final\') '*.csv']);
    copyfile(fullfile(dlcPath, csvFile(1).name), destination);
    csvFileNew = fullfile(destination, csvFile(1).name);
    
    cd(sessions{s})

    %% Load data
    if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
        disp('Tracking already detected! Loading file.');
        file = dir([basepath filesep '*.Tracking.Behavior.mat']);
        copyfile(file(1).name, destination);
        load(fullfile(destination, file(1).name));
    end
    
    if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
        disp('Behavior already detected! Loading file.');
        file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
        copyfile(file(1).name, destination);
        load(fullfile(destination, file(1).name));
    end
    
    if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
        disp('Spikes already detected! Loading file.');
        file = dir([basepath filesep '*.spikes.cellinfo.mat']);
        copyfile(file(1).name, destination);
        load(fullfile(destination, file.name));
    end
        
    %% Import coordinates, session and trial timestamps
    % if strcmp(csvFile, "IZ40_220714_sess18")
    %     csvFile = 'F:\DLC_ACgN_HeadDirection-athina-2023-03-03\videos\iter-1\IZ40_220714_sess18DLC_resnet50_DLC_ACgN_HeadDirectionMar3shuffle1_500000_filtered.csv';
    %     disp(csvFile);
    % end
    opts = detectImportOptions(csvFileNew, ...
        'NumHeaderLines', 3, 'VariableNamingRule', 'preserve');
    DLC = readtable(csvFileNew, opts);
    
    numTrials = size(behavTrials.timestamps, 1);   
    numFrames = size(DLC, 1);
    
    angle = zeros(numFrames, 1);
    centroid_y = zeros(numFrames, 1);
    rotations = zeros(numTrials, 1);
       
    % Define the line line_track = two ends of track 
    % 14 is the midpoint of x coordinate in the video in pixels and 460 is 
    % the length of y in the video in pixels. 460 is the beginning of the
    % track. 
    line_track = [14, 460, 0] - [14, 0, 0];
    
    for i = 1:numFrames
    
        x1 = DLC.Var2(i); % leftear
        x2 = DLC.Var5(i); % rightear
        x3 = DLC.Var8(i); % snout
        x4 = DLC.Var11(i); % tailbase 
    
        y1 = DLC.Var3(i);
        y2 = DLC.Var6(i);
        y3 = DLC.Var9(i);
        y4 = DLC.Var12(i); 
    
        p1 = DLC.Var4(i);
        p2 = DLC.Var7(i);
        p3 = DLC.Var10(i);
        p4 = DLC.Var13(i); 
    
    
        %% Fix incorrectly tracked values
        % Set NaN all values where probability in at least 2 coordinates is < 0.9.
        % This is done for longer "chunks" of the video (across 10 frames).
    %         if i > 5 && i < numFrames-5
    %             if ((p1 < 0.9 && p2 < 0.9)||(p1 < 0.9 && p3 < 0.9)||(p1 < 0.9 && p4 < 0.9)...
    %                     ||(p2 < 0.9 && p3 < 0.9)||(p2 < 0.9 && p4 < 0.9)||(p3 < 0.9 && p4 < 0.9)...
    %                     ||(p1 < 0.9 && p2 < 0.9 && p3 < 0.9)||(p1 < 0.9 && p2 < 0.9 && p4 < 0.9)...
    %                     ||(p1 < 0.9 && p3 < 0.9 && p4 < 0.9)||(p2 < 0.9 && p3 < 0.9 && p4 < 0.9)...
    %                     ||(p1 < 0.9 && p2 < 0.9 && p3 < 0.9 && p4 < 0.9)) && ...
    %                     ((median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var7((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var10((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9))||...
    %                     (median(DLC.Var4((i-5):(i+5)) < 0.9) && median(DLC.Var7((i-5):(i+5)) < 0.9) && median(DLC.Var10((i-5):(i+5)) < 0.9) && median(DLC.Var13((i-5):(i+5)) < 0.9)))
    %         
    %                 
    %                 x1 = NaN;
    %                 x2 = NaN;
    %                 x3 = NaN;
    %                 x4 = NaN; 
    %                 y1 = NaN; 
    %                 y2 = NaN; 
    %                 y3 = NaN; 
    %                 y4 = NaN;                 
    %             end
    %         end
    
        % Smooth out frames with wrong values. This is done for a few frames (a
        % couple of consecutive frames) for each coordinate. 
        if i > 2 && i < numFrames-2
            if ((p1 < 0.9) && (DLC.Var4(i-1) > 0.9))             
                x1 = median(DLC.Var2((i-2):(i+2)));
                y1 = median(DLC.Var3((i-2):(i+2)));
        
            elseif ((p2 < 0.9) && (DLC.Var7(i-1) > 0.9))
                x2 = median(DLC.Var5((i-2):(i+2)));
                y2 = median(DLC.Var6((i-2):(i+2)));
            
            elseif ((p3 < 0.9) && (DLC.Var10(i-1) > 0.9)) 
                x3 = median(DLC.Var8((i-2):(i+2)));
                y3 = median(DLC.Var9((i-2):(i+2)));
        
            elseif ((p4 < 0.9) && (DLC.Var13(i-1) > 0.9))
                x4 = median(DLC.Var11((i-2):(i+2)));
                y4 = median(DLC.Var12((i-2):(i+2)));
            end
        end
        
        %% Calculate the angle of the head compared to the y-axis i.e., the 
        % line connecting the beginning and end of the linear track. 
    
        % Find midpoint between the two ears
        x_ear = mean([x1, x2]);
        y_ear = mean([y1, y2]);
        
        % Define the line line_head = nose - earMidline 
        line_head = [x3, y3, 0] - [x_ear, y_ear, 0];
                
        % Calculate the angle
        cross_product = cross(line_head, line_track);
    
        angle(i) = sign(cross_product(3)) * atan2d(norm(cross(line_head, ...
            line_track)), dot(line_head, line_track));
        
        centroid_y(i) = mean([y1,y2,y3,y4]);
        
    end
    
    %% Calculate direction of rotations -1 = counterclockwise, 1 = clockwise
    % The direction of rotations is determined based on when the head
    % direction passes a threshold angle, which means that the mice are in
    % the process of rotating. We then find the sign of the angle when the
    % threshold is crossed. 
    
    angle_thres = 100;
    angle_thres_idx = zeros(numTrials, 1);
    
    dtime = mean(diff(tracking.timestamps));
    
    win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
    spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);
    
    spkMat = spkData.data';
    timestamps = spkData.timestamps';

    for i = 1:numTrials
        [~,startIdx(i)] = min(abs(timestamps - behavTrials.timestamps(i,1))); 
        [~,endIdx(i)] = min(abs(timestamps - behavTrials.timestamps(i,2))); 

%         angle_thres_idx(i) = find(abs(angle(startIdx(i):endIdx(i))) ...
%             < angle_thres, 1, 'first');
%         angle_thres_idx(i) = angle_thres_idx(i) + startIdx(i);
%               
%         rotations(i) = sign(angle(angle_thres_idx(i)));     
    end
    
    %% Store the angle and rotations in the trackingFile
    tracking.position.angle = angle(:);
    tracking.position.rotations = rotations(:);
    save(trackingFile, 'tracking');
    
%     trackRotation.direction = rotations(:);
%     trackRotation.timestamps = tracking.timestamps(angle_thres_idx(:));
%     trackRotation.linTrial = behavTrials.linTrial;
%     trackRotation.lickLoc = behavTrials.lickLoc + 1;
%     save(fullfile(destination, strcat(extractAfter(sessions{s},'Final\'), '.rotation.mat')), 'trackRotation');
    
    
    %% Plot direction of rotations as a function of lick location
    lickLoc = behavTrials.lickLoc + 1;
    ports = 1:6;
    
    % Tone trials
    [freqsC_tone,~] = histcounts(lickLoc(rotations == 1 & ...
        behavTrials.linTrial == 0), 'NumBins', length(ports));
    [freqsCC_tone,~] = histcounts(lickLoc(rotations == -1 & ...
        behavTrials.linTrial == 0), 'NumBins', length(ports));
    
    % No-tone trials
    [freqsC_notone,~] = histcounts(lickLoc(rotations == 1 & ...
        behavTrials.linTrial == 1), 'NumBins', length(ports));
    [freqsCC_notone,~] = histcounts(lickLoc(rotations == -1 & ...
        behavTrials.linTrial == 1), 'NumBins', length(ports));
    
    figure;
    subplot(1,2,1);
    bar(ports, [freqsC_tone' freqsCC_tone'], 'stacked')
    legend('Clockwise', 'Countercloskwise')
    xlabel('Ports')
    ylabel('Number of rotations - Tone')
    
    subplot(1,2,2);
    bar(ports, [freqsC_notone' freqsCC_notone'], 'stacked')
    legend('Clockwise', 'Countercloskwise')
    xlabel('Ports')
    ylabel('Number of rotations - No-tone')
    
    saveas(gcf, fullfile(behaviorPath, 'rotations'), 'fig');
    saveas(gcf, fullfile(behaviorPath, 'rotations'), 'png');
    
    %% Plot the head direction in angles
    %     figure;
    %     x = [1:numFrames]; % 60 sec/min * 30 frames/sec = 1800 frames/min
    % %     idx1 = find(abs(angle)>110 | abs(angle)<70);
    % %     idx2 = find(abs(angle)<=110 & abs(angle)>=70);
    %     
    %     hold on
    %     yyaxis left
    % %     scatter(x, angle, 'black', 'filled')
    % %     scatter(x(endFrames), angle(endFrames), 'red', 'filled');
    %     scatter(numTrials, rotations, 'blue', 'filled');
    % %     plot(x, speed, '-')
    % %     scatter(x(idx1), angle(idx1), 'black', 'filled')
    % %     scatter(x(idx2), angle(idx2), 'red', 'filled')
    % 
    %     yyaxis right       
    %     plot(x, centroid_y)
    % 
    %     yyaxis left
    %     title('Head direction relative to the linear track axis', FontSize=20);
    %     xlabel('Time (min)', FontSize=18);
    %     ylabel('Head direction (degrees)', FontSize=18);
    %     ylim([-180,180]);
    % 
    %     yyaxis right 
    %     ylabel('Y position (pixels)', FontSize=18);
    %     ylim([10,460]);
    % 
    %     hold off
    % end
    cd(root)
        
end

end