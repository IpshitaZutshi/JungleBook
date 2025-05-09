function combineFrontandTopvideos(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;

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

 if ~isempty(dir([basepath filesep 'Top*csv']))
    csvFile = dir([basepath filesep 'Top*csv']); 
    csvFile2 = dir([basepath filesep 'Front*csv']);       
else
    warning('No csv file!!');
    tracking = [];
    return
 end

 %% Load videos and csv files
videoObj = VideoReader([aviFile '.avi']);   % get video
numFrames = get(videoObj, 'NumFrames');

videoObj2 = VideoReader([aviFile2 '.avi']);   % get video
numFrames2 = get(videoObj2, 'NumFrames');

fid = fopen(csvFile.name, 'r');
behavInfo = textscan(fid,'%s','delimiter','\n');
TS = cellfun(@(x)strsplit(x,'_'),behavInfo{1},'UniformOutput',false);       
TS = vertcat(TS{:});
TS = str2double(TS);

fid2 = fopen(csvFile2.name, 'r');
behavInfo2 = textscan(fid2,'%s','delimiter','\n');
TS2 = cellfun(@(x)strsplit(x,'_'),behavInfo2{1},'UniformOutput',false);       
TS2 = vertcat(TS2{:});
TS2 = str2double(TS2);

fs = round(1000./mode(diff(TS)));

%% open the video writer
fig = figure;
set(fig,'Position',[990 333 807 564])
set(fig,'Color','k')

subplot(3,3,[1 4 7]);
set(gca,'visible','off')
axis off

subplot(3,3,[2 3 5 6 8 9]);
set(gca,'visible','off')
axis off

writerObj = VideoWriter('CombinedVideo','MPEG-4');
writerObj.FrameRate = fs;
open(writerObj);

%% Align frames
for ii = 18000:length(TS2) % Generally TS2 has fewer frames that TS, so align TS to TS2
    %Find TS index closest to current TS2 frame
    [~,idx] = min(abs(TS-TS2(ii)));
    frame2 = read(videoObj2,ii);
    frame1 = read(videoObj,idx);

    subplot(3,3,[1 4 7]);
    imagesc(frame1(120:960,580:840,:))  
    axis off

    subplot(3,3,[2 3 5 6 8 9]);
    imagesc(frame2(200:770,1:900,:)) 
    camroll(90)
    axis off 

    M = getframe(fig);   
    writeVideo(writerObj, M);
end

close(writerObj)
end
