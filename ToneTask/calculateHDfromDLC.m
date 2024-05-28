function [angle] = calculateHDfromDLC(varargin)

% This function calculates the head direction of the animals in the form of 
% an angle w.r.t. the x-direction using the tracked coordinates from 
% DeepLabCut. It also determines the direction of turn of the animals after
% the lick at the home port to start a new trial. 

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'plotfig',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
plotfig = p.Results.plotfig;

cd(basepath)

if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Tracking already detected! Loading file.');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

% Extract DLC path
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));

for ii = 1:size(MergePoints.foldernames,2)
     if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'test*']))   
        filename = dir([basepath filesep MergePoints.foldernames{ii} filesep 'test*' '_DLC.csv']);
        videoname = dir([basepath filesep MergePoints.foldernames{ii} filesep 'test*' '.avi']);
        videoFile = strcat(videoname.folder,'\',videoname.name);
        videoObj = VideoReader(videoFile);

        csvFile = strcat(filename.folder,'\',filename.name);

        opts = detectImportOptions(csvFile, ...
            'NumHeaderLines', 3, 'VariableNamingRule', 'preserve');
        DLC = readtable(csvFile, opts);

        numTrials = size(behavTrials.timestamps, 1);   
        numFrames = size(DLC, 1);
        
        angle = zeros(numFrames, 1);
        centroid_y = zeros(numFrames, 1);
        rotations = zeros(numTrials, 1);
       
        % Define the line line_track = two ends of track 
        % 14 is the midpoint of x coordinate in the video in pixels and 460 is 
        % the length of y in the video in pixels. 460 is the beginning of the
        % track. 
        line_track = atan2((460-0),(14-14));

        if plotfig
            figure
        end
        for i = 1:numFrames
    
            x1 = DLC.Var11(i); % rightear
            x2 = DLC.Var14(i); % leftear
            x3 = DLC.Var2(i); % snout
            x4 = DLC.Var5(i); % tailbase 

            y1 = DLC.Var12(i);
            y2 = DLC.Var15(i);
            y3 = DLC.Var3(i);
            y4 = DLC.Var6(i); 
        
            p1 = DLC.Var13(i);
            p2 = DLC.Var16(i);
            p3 = DLC.Var4(i);
            p4 = DLC.Var7(i); 

            % Smooth out frames with wrong values. This is done for a few frames (a
            % couple of consecutive frames) for each coordinate. 
            if i > 2 && i < numFrames-2
                if ((p1 < 0.9) && (DLC.Var13(i-1) > 0.9))             
                    x1 = median(DLC.Var11((i-2):(i+2)));
                    y1 = median(DLC.Var12((i-2):(i+2)));
            
                elseif ((p2 < 0.9) && (DLC.Var16(i-1) > 0.9))
                    x2 = median(DLC.Var14((i-2):(i+2)));
                    y2 = median(DLC.Var15((i-2):(i+2)));
                
                elseif ((p3 < 0.9) && (DLC.Var7(i-1) > 0.9)) 
                    x3 = median(DLC.Var5((i-2):(i+2)));
                    y3 = median(DLC.Var6((i-2):(i+2)));
            
                elseif ((p4 < 0.9) && (DLC.Var13(i-1) > 0.9))
                    x4 = median(DLC.Var11((i-2):(i+2)));
                    y4 = median(DLC.Var12((i-2):(i+2)));
                end
            end
            
            if plotfig
                temp_frames = read(videoObj,i);                 
                imagesc(temp_frames);
                hold on
                scatter(x1,y1,[],'b')
                scatter(x2,y2,[],'g')
                scatter(x3,y3,[],'r')
                scatter(x4,y4,[],'k')
            end            
            %% Calculate the angle of the head compared to the y-axis i.e., the 
            % line connecting the beginning and end of the linear track. 
        
            % Find midpoint between the two ears
            x_ear = mean([x1, x2]);
            y_ear = mean([y1, y2]);
            
            % Define the line line_head = nose - earMidline 
            %line_head = [x3, y3, 0] - [x_ear, y_ear, 0];
                  
            line_head = atan2((y3-y_ear),(x3-x_ear));

            % Calculate the angle
            % cross_product = cross(line_head, line_track);
        
            % angle(i) = sign(cross_product(3)) * atan2d(norm(cross(line_head, ...
            %     line_track)), dot(line_head, line_track));


            angle(i) = rad2deg(line_head);%-rad2deg(line_track);
            centroid_y(i) = mean([y1,y2,y3,y4]);
                  
        end
    
     end     
end

%% Calculate direction of rotations -1 = counterclockwise, 1 = clockwise
% The direction of rotations is determined based on when the head
% direction passes a threshold angle, which means that the mice are in
% the process of rotating. We then find the sign of the angle when the
% threshold is crossed. 

angle_thres = 100;
angle_thres_idx = zeros(numTrials, 1);      

for i = 1:numTrials
    [~,startIdx(i)] = min(abs(tracking.timestamps - behavTrials.timestamps(i,1))); 
    [~,endIdx(i)] = min(abs(tracking.timestamps - behavTrials.timestamps(i,2))); 

    angle_thres_idx(i) = find(abs(angle(startIdx(i):endIdx(i))) ...
        < angle_thres, 1, 'first');
    angle_thres_idx(i) = angle_thres_idx(i) + startIdx(i);

    rotations(i) = sign(angle(angle_thres_idx(i)));     
end

%     %% Store the angle and rotations in the trackingFile
%     tracking.position.angle = angle(:);
%     tracking.position.rotations = rotations(:);
%     save(trackingFile, 'tracking');
% 
% %     trackRotation.direction = rotations(:);
% %     trackRotation.timestamps = tracking.timestamps(angle_thres_idx(:));
% %     trackRotation.linTrial = behavTrials.linTrial;
% %     trackRotation.lickLoc = behavTrials.lickLoc + 1;
% %     save(fullfile(destination, strcat(extractAfter(sessions{s},'Final\'), '.rotation.mat')), 'trackRotation');
% 
% 
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
% 
%     saveas(gcf, fullfile(behaviorPath, 'rotations'), 'fig');
%     saveas(gcf, fullfile(behaviorPath, 'rotations'), 'png');
    
%    Plot the head direction in angles
        figure;
        x = [1:numFrames]; % 60 sec/min * 30 frames/sec = 1800 frames/min
        idx1 = find(abs(angle)>110 | abs(angle)<70);
        idx2 = find(abs(angle)<=110 & abs(angle)>=70);

        hold on
        yyaxis left
        scatter(x, angle, 'black', 'filled')
        %scatter(x(endFrames), angle(endFrames), 'red', 'filled');
        scatter(numTrials, rotations, 'blue', 'filled');
        plot(x, speed, '-')
        scatter(x(idx1), angle(idx1), 'black', 'filled')
        scatter(x(idx2), angle(idx2), 'red', 'filled')

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
%    cd(root)
        
end
