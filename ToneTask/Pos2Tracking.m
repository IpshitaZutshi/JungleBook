
function [tracking] = Pos2Tracking(aviFile,varargin)
% Get position tracking
%
% USAGE
%
%   [behavior] = Pos2Tracking(varagin)
%
% INPUTS
%   aviFile     Avi format video. If not provided, look for it in the
%   basepath f
%   It requires an avi format video and a digitalin.dat file in the
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
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'artifactThreshold',10,@isnumeric);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'verbose',false,@islogical);
addParameter(p,'thresh',25,@isnumeric)
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
convFact = p.Results.convFact;
forceReload = p.Results.forceReload;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
bazlerTtl = p.Results.bazlerTTL;
saveMat = p.Results.saveMat;
% RGBChannel = p.Results.RGBChannel;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

if ~exist('aviFile') || isempty(aviFile)
    if ~isempty(dir([basepath filesep 'test*avi']))
        aviFile = dir([basepath filesep 'test*avi']); 
        aviFile = erase(aviFile.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end
end

if ~exist([basepath filesep aviFile '.mat'],'file') || forceReload
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
    frames = [];
    tic
    f = waitbar(0,'Getting frames...');
    for ii = 1:numFrames
        waitbar(ii/numFrames,f)
        temp_frames = read(videoObj,ii);        % get all frames
        frames(:,:,ii) = rgb2gray(temp_frames);          % convert to grayscale
    end
    close(f)
    toc
    
    if saveFrames
        disp('Saving frames...');
        save([basepath filesep aviFile '.mat'],'frames','-v7.3');
    end
else
    disp('Loading frames from mat file...');
    load([basepath filesep aviFile '.mat'],'frames');
end

if isempty(convFact)                             % if convFact not provided, normalize to 1 along the longest axis
    convFact = 1/max([size(frames,1) size(frames,2)]);
    artifactThreshold = Inf;
end
xMaze = [0 size(frames,2) * convFact];
yMaze = [0 size(frames,1) * convFact];


%% DETECT MOUSE POSITION
disp('Detect mouse position...');
if ~verbose
    tic
    f = waitbar(0,'Detecting mouse position...');
    for ii = 1:size(frames,3)
        waitbar(ii/size(frames,3),f)
        fr = frames(:,:,ii);
        bin_fr = ~(imbinarize(double(fr),thresh)); %
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
else
    h1 = figure;
    hold on
    tic
    for ii = 1:size(frames,3)
        fr = frames(:,:,ii);
        bin_fr = ~(imbinarize(double(fr),thresh)); %
        bin_fr = bwareafilt(bin_fr,1);
        stats_fr = regionprops(bin_fr);
        maxBlob = find([stats_fr.Area]== max([stats_fr.Area]),1);
        if ~isempty(maxBlob)
            sz_fr(ii) = stats_fr(maxBlob).Area;
            Rr_x(ii) = stats_fr(maxBlob).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlob).Centroid(2);
            cla
            imagesc(fr)
            plot(Rr_x(ii),Rr_y(ii),'or')
            drawnow;
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
    end
    toc
    close(h1);
end

pos = [Rr_x; Rr_y]';

%% postprocessing of LED position 
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
xlim(xMaze); ylim(yMaze);
mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');
if ~verbose
    close(h2);
end


%% Get basler TTL
% digital_in legend: 1. Basler, 2. maze LEd, 3. Left Alternation, 4.Righ
% Alternation, 5. Home Delay, 6. Is alternation forzed?
if isempty(bazlerTtl)
    digitalIn = bz_getDigitalIn;
    bazlerTtl = digitalIn.timestampsOn{1};
end

% match basler frames con ttl pulses
if length(bazlerTtl) == length(x)
    disp('Number of frames match!!');
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) <= length(x) + 15 * 1 
    fprintf('%3.i frames were dropped, probably at the end of the recording. Skipping... \n',...
        length(bazlerTtl) - length(x));
    bazlerTtl = bazlerTtl(1:length(x));
elseif length(bazlerTtl) < length(x) && (length(x)-length(bazlerTtl)) < 30 * 4
    fprintf('%3.i video frames without TTL... Was the recording switched off before the camera?. Skipping... \n',...
        length(x) - length(bazlerTtl));
    x = x(1:length(bazlerTtl));
    y = y(1:length(bazlerTtl));
    vx = vx(1:length(bazlerTtl));
    vy = vy(1:length(bazlerTtl));
    v = v(1:length(bazlerTtl));
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
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

% maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 
