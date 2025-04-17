
function [tracking] = Pos2Tracking_Cue(varargin)
% Get position tracking
%
% USAGE
%
%   [behavior] = Pos2Tracking_Cue(varagin)
%
% INPUTS
%   excelfile      Look for it in the basepath
%   It requires an excel file (deeplabcut Output) and a digitalin.dat file in the
%   basepath folder.
% 
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in
%                   buzcode format.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveFrames     - Creates mat file containin all frames (default false).
%   forceReload    - default false.
%   bazlerTTL      - Rx1 with timestamps from bazler ttl pulses. If not
%                   provided try to extract ttl pulses from digitalIn
%                   channel 1. If it fails, gives video time.
%   saveMat        - default true
%   artifactThreshold - max allow movements per frame (in cm, default 3).
%                   Disabled if not convFact is provided.
% 
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%   samplingRate    - in Hz
%   averageFrame    - 
%   
%
%   Manu Valero 2019
% TO DO: Generalize for non-T maze behaviour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'convFact',0.133,@isnumeric); % 0.1149
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'thresh',40,@isnumeric)
addParameter(p,'artifactThreshold',10,@isnumeric)
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'useDLC',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
convFact = p.Results.convFact;
forceReload = p.Results.forceReload;
artifactThreshold = p.Results.artifactThreshold;
thresh = p.Results.thresh;
saveMat = p.Results.saveMat;
useDLC = p.Results.useDLC;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
else

    if ~isempty(dir([basepath filesep 'Top*avi']))
        aviFile = dir([basepath filesep 'Top*avi']); 
        aviFile = erase(aviFile.name,'.avi');

        aviFile2 = dir([basepath filesep 'Front*avi']); 
        aviFile2 = erase(aviFile2.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end    

    if useDLC
        csvFile = dir([basepath filesep 'DLC_Top*csv']);
        csvFile = csvFile.name;
    end
end

%% Save the top camera frames
if ~exist([basepath filesep aviFile '.mat'],'file') || forceReload
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
end


%% DETECT MOUSE POSITION
disp('Detect mouse position...');

if useDLC
    opts = detectImportOptions(csvFile, ...
        'NumHeaderLines', 3, 'VariableNamingRule', 'preserve');
    DLC = readtable(csvFile, opts);
    x1 = DLC.Var2; % leftear
    y1 = DLC.Var3;
    p = DLC.Var4;
    x1(p<0.9995) = nan;
    y1(p<0.9995) = nan;
    pos = [x1 y1];
else
    tic
    f = waitbar(0,'Detecting mouse position...');
    for ii = 1:numFrames
        waitbar(ii/numFrames,f)
        temp_frames = read(videoObj,ii);
        fr = rgb2gray(temp_frames(100:950,550:850,:));
    
        % Apply a mask around regions that are not off interest
        mask_fr = fr;
        mask_fr(1:735,1:135) = 255;
        mask_fr(1:735,195:end) = 255;
        mask_fr(800:end,:) = 255;
    
        bin_fr = ~(imbinarize(double(mask_fr),thresh)); %
        bin_fr = bwareafilt(bin_fr,1);
        stats_fr = regionprops(bin_fr);
        maxBlob = find([stats_fr.Area]== max([stats_fr.Area]),1);
        if ~isempty(maxBlob)
            sz_fr(ii) = stats_fr(maxBlob).Area;
            Rr_x(ii) = stats_fr(maxBlob).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlob).Centroid(2);
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
    end
    close(f)
    toc
    
    pos = [Rr_x; Rr_y]';
end

%% Load digitalInputs
items = dir;

% Step 2: Filter for directories whose names start with 'T'
folders = items([items.isdir] & startsWith({items.name}, 'T'));
folderNames = {folders.name};
load(strcat(folderNames{1},'\digitalIn.events.mat'))
bazlerTtl = digitalIn.timestampsOn{12};


%% postprocessing of LED position 
fs = 1./mode(diff(bazlerTtl));
pos = pos * convFact;                                   % cm or normalized
art = find(sum(abs(diff(pos))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos(art,:) = NaN;

xt = linspace(0,size(pos,1)/fs,size(pos,1));            % kalman filter
[t,x,y,vx,vy,ax,ay] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',xt,0);
art = find(sum(abs(diff([x y]))>artifactThreshold,2))+1;
art = [art - 2 art - 1 art art + 1 art + 2];
x(art(:)) = NaN; y(art(:)) = NaN;
F = fillmissing([x y],'linear');
x = F(:,1); y = F(:,2);

h2 = figure;
hold on
scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
caxis([t(1) t(end)])
xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');

tracking.framesDropped = length(bazlerTtl) - length(x);

% match basler frames con ttl pulses
if length(bazlerTtl) == length(x)
    disp('Number of frames match!!');
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) > length(x) + 5 && ~isempty(dir([basepath filesep 'Top*csv'])) %% if a csv file exists
    
    txtFile = dir([basepath filesep 'Top*csv']);
    fid = fopen(txtFile.name, 'r');
    behavInfo = textscan(fid,'%s','delimiter','\n');
    TS = cellfun(@(x)strsplit(x,'_'),behavInfo{1},'UniformOutput',false);       
    TS = vertcat(TS{:});
    TS = str2double(TS);
    TS = round(TS-TS(1));
    if length(TS) ~= length(x)
        disp('Yikes!! CSV file does not match');
    end
    %Determine how the dropped frames align
    tsBas = round((bazlerTtl-bazlerTtl(1))*1000);    
    Xcorr(1:length(tsBas)) = nan;
    Ycorr(1:length(tsBas)) = nan;
    
    for tt = 1:length(TS)
        [~, idx] = min(abs(tsBas-TS(tt)));
        Xcorr(idx) = x(tt);
        Ycorr(idx) = y(tt);
    end
    
    fprintf('%3.i frames were dropped, fixed using csv file \n',...
        length(bazlerTtl) - length(x));
    x = Xcorr;
    y = Ycorr;
    
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) <= length(x) + 30 * 1 
    fprintf('%3.i frames were dropped, probably at the end of the recording. Skipping... \n',...
        length(bazlerTtl) - length(x));
    bazlerTtl = bazlerTtl(1:length(x));
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) > length(x) + 30 * 1 
    fprintf('%3.i frames were dropped!! Examine the recording carefully!... \n',...
        length(bazlerTtl) - length(x));
elseif length(bazlerTtl) < length(x) 
    fprintf('%3.i video frames without TTL... Was the recording switched off before the camera?. Skipping... \n',...
        length(x) - length(bazlerTtl));
    x = x(1:length(bazlerTtl));
    y = y(1:length(bazlerTtl));
    vx = vx(1:length(bazlerTtl));
    vy = vy(1:length(bazlerTtl));
    ax = ax(1:length(bazlerTtl));
    ay = ay(1:length(bazlerTtl)); 
elseif isempty(bazlerTtl)
    bazlerTtl = xt;
end

[~,fbasename,~]=fileparts(pwd);

% Get instantaneous speed
[~,~,~,vx,vy,~,~] = KalmanVel(x,y,bazlerTtl,2);
v = sqrt(vx.^2+vy.^2);


tracking.position.x =x;
tracking.position.y = y;
tracking.position.vx =vx;
tracking.position.vy = vy;
tracking.position.v =v;

tracking.description = '';
tracking.timestamps = bazlerTtl;
tracking.originalTimestamps = [];
tracking.folder = fbasename;
tracking.samplingRate = fs;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end