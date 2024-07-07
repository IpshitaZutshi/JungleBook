function Fig4_plotRastermapUMAPDecoding

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrow = 6;
numcol = 5;

RdPu=cbrewer('seq', 'RdPu', 11);
Purples=cbrewer('seq', 'PuBu', 11);

%% Go to the directory
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)

%% First plot UMAP
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230714_sess28\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
umap_name = 'behavior_speed_1_smooth_5';
A = 0;
E = 90;
% First, entire manifold
TRIAL_TYPE = [0 1 2 3 4 5 6 7 8];
manifoldPlot('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

% Next, tone trials
TRIAL_TYPE = [0:5];
%col = 'jet';
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

A = 23;
E = 1.55;
manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'addPosPlot',true,'addFreq',true,'poscutOff',2,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'plotcolorbar',false,'dotSize',2)

%% Calculate rastermap
% Load files
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*cell_metrics.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file(1).name);

% Generate original spike matrix
dtime = mean(diff(tracking.timestamps));
    
% Excluding the last trial
win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%Only select pyramidal cells
logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
%Only select cells with a rate> 0.1 Hz
rate = sum(spkMat,2)./(length(timestamps)*(1/30));
logVector2 = rate>0.1;
keepCells = logicalVector'& logVector2;
clear spkMat

% Now generate rastermapData
filename = 'IZ48_230714_sess28.rastermapData';
%filename = 'IZ47_230710_sess25.rastermapData';
directory_path  = 'C:\Data\Rastermap\preProcessedData\';
load(strcat(directory_path,filename))

data = zscore(spkMat,[],1);

dataNdArray = py.numpy.array(data);

pyrun("from rastermap import Rastermap") %load interpreter, import main function
rmModel = pyrun("model = Rastermap(n_clusters=50, n_PCs=200, locality=0, time_lag_window=3, grid_upsample = 0).fit(spks)", "model", spks=dataNdArray);
sortIdx = int16(py.memoryview(rmModel.isort.data)) + 1; %back to MATLAB array, 1-indexing

%sortIdx_new = sortIdx;
sortIdx_new = [sortIdx(61:end) sortIdx(1:60)];

%% Plot total rastermap

subplot(numrow,numcol,6:10)

% The sorting is the same but using the output from python because the
% grayscale looks better
rastermapDir = 'C:\Data\Rastermap\AnalyzedData';
A = imread(strcat(rastermapDir,'\Allrastermap.png'));
image(A(150:1350,220:2200,:))
axis off

%% Also extract decoding
decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\IZ48\Final\IZ48_230714_sess28\py_data\theta_decoding_lickLoc_y';
file_name = 'up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
posterior_goal =  ncread(strcat(decodingPath,'\',file_name),'x_position') ;
posterior_pos =  ncread(strcat(decodingPath,'\',file_name),'y_position') ;
post_time =  ncread(strcat(decodingPath,'\',file_name),'time') ;
post_pos =  ncread(strcat(decodingPath,'\',file_name),'y_position_value') ;
post_goal =  ncread(strcat(decodingPath,'\',file_name),'x_position_value') ;

file_name = 'change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';
load(strcat(decodingPath,'\',file_name));

%% Plot a Tone port 6 trial
trialNum = 55;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,numrow,numcol,11,fig2,pos,trialEnd,5798.5,5802.3,behavTrials.toneGain(trialNum))
title('Port 6')
%xlim([5798.5 5802.3])

% Plot manifold
TRIAL_TYPE = [0:5];
A = 23;
E = 1.55;
manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',8,'poscutOff',2,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5798.5 trialEnd])

%Plot position decoding
[minTS,idxDec1] = min(abs(post_time-5798.5));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-5802.3));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end
%Plot position decoding
ax1 = subplot(numrow,numcol,21);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
box off
caxis([0 0.4])
xlim([5798.5 5802.3])

%Plot goal decoding
ax1 = subplot(numrow,numcol,26);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on
plotChangePoints(fig2,numrow,numcol,26,change_point,trial,5798.5,5802.3,post_time,trialNum,posterior_goal)
xlim([5798.5 5802.3])

%% Plot a Tone port 5 trial
trialNum = 61;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,numrow,numcol,12, fig2,pos,trialEnd,5872,5876,behavTrials.toneGain(trialNum))
title('Port 5')

manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',8,'poscutOff',2,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',2,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5872,trialEnd])

%Plot position decoding
[minTS,idxDec1] = min(abs(post_time-5872));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-5876));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end
%Plot position decoding
ax1 = subplot(numrow,numcol,22);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
caxis([0 0.4])
box off
xlim([5872 5876])

%Plot goal decoding
ax1 = subplot(numrow,numcol,27);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on
plotChangePoints(fig2,numrow,numcol,27,change_point,trial,5872,5876,post_time,trialNum,posterior_goal)
xlim([5872 5876])

%% Plot a Tone port 4 trial
trialNum = 60;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,numrow,numcol,13, fig2,pos,trialEnd,5860.5,5863.5,behavTrials.toneGain(trialNum))
title('Port 4')

manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',8,'poscutOff',2,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',3,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5860.5 trialEnd])

%Plot position decoding
[minTS,idxDec1] = min(abs(post_time-5860.5));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-5863.5));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end
%Plot position decoding
ax1 = subplot(numrow,numcol,23);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
caxis([0 0.4])
box off
xlim([5860.5 5863.5])

%Plot goal decoding
ax1 = subplot(numrow,numcol,28);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on
plotChangePoints(fig2,numrow,numcol,28,change_point,trial,5860.5,5863.5,post_time,trialNum,posterior_goal)
xlim([5860.5 5863.5])

%% Plot a Tone port 3 trial
trialNum = 25;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,numrow,numcol,14, fig2,pos,trialEnd,5455.9,5458.9,behavTrials.toneGain(trialNum))
title('Port 3')

manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',8,'poscutOff',2,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',4,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5455.9 trialEnd])

%Plot position decoding
[minTS,idxDec1] = min(abs(post_time-5455.9));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-5458.9));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end
%Plot position decoding
ax1 = subplot(numrow,numcol,24);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
caxis([0 0.4])
box off
xlim([5455.9 5458.9])

%Plot goal decoding
ax1 = subplot(numrow,numcol,29);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on
plotChangePoints(fig2,numrow,numcol,29,change_point,trial,5455.9,5458.9,post_time,trialNum,posterior_goal)
xlim([5455.9 5458.9])

%% Plot a Tone port 2 trial
trialNum = 34;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,numrow,numcol,15,fig2,pos,trialEnd,5551.5,5554.5,behavTrials.toneGain(trialNum))
title('Port 2')

manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',6,'poscutOff',2,'speedThresh',2,...
     'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',5,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5551.5 trialEnd])

%Plot position decoding
[minTS,idxDec1] = min(abs(post_time-5551.5));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-5554.5));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end
%Plot position decoding
ax1 = subplot(numrow,numcol,25);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
caxis([0 0.4])
box off
colorbar
xlim([5551.5 5554.5])

%Plot goal decoding
ax1 = subplot(numrow,numcol,30);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
colorbar
hold on
plotChangePoints(fig2,numrow,numcol,30,change_point,trial,5551.5,5554.5,post_time,trialNum,posterior_goal)
xlim([5551.5 5554.5])

% %% Plot a Tone port 1 trial
% trialNum = 42;
% raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
% time_seg = timestamps(eventVariables.trialNum==trialNum);
% pos = constVariables.y(eventVariables.trialNum==trialNum);
% trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
% plotRaster(raster_seg,time_seg,numrow,numcol,15, fig2,pos,trialEnd,5652,5655, behavTrials.toneGain(trialNum))
% title('Port 1')
% 
% manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',8,...
%     'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',5,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5652 trialEnd])
% 
% 
% %Plot position decoding
% [minTS,idxDec1] = min(abs(post_time-5652));
% if minTS>0.5 %Take the next index
%     idxDec1 = idxDec1+1;
% end
% [minTS,idxDec2] = min(abs(post_time-5655));
% if minTS>0.5 %Take the previous index
%     idxDec2 = idxDec2-1;
% end
% %Plot position decoding
% ax1 = subplot(numrow,numcol,25);
% imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
% set(gca,'YDir','normal')
% colormap(ax1,viridis)
% caxis([0 0.4])
% 
% %Plot goal decoding
% ax1 = subplot(numrow,numcol,30);
% imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
% set(gca,'YDir','normal')
% colormap(ax1,viridis)

% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure4A_rasterUMAP.png'));
saveas(gcf,strcat(expPath,'Figure4A_rasterUMAP.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure4A_rasterUMAP.fig'));

%% Make a new figure with average statistics
%linkDecodingandUMAP

end


function plotRaster(spkMat,time_seg,numrows,numcol,plotloc,fighandle,pos,trialEnd,timeStart,timeEnd,toneGain)

ax1 = subplot(numrows, numcol, plotloc,'Parent',fighandle);
hold on;
[~,idxStart] = min(abs(time_seg-timeStart));
[~,idxEnd] = min(abs(time_seg-timeEnd));

time_seg_trunc = time_seg(idxStart:idxEnd);
[~,idxEndTrial] = min(abs(time_seg_trunc-trialEnd));

% get frequency estimate
gain = [122/10.5, 122/32.6 122/55.53 122/79.62 122/102.79 122/122];
freqExp = log10(25000/1000);
freq = (pos(idxStart:idxEnd)*gain(toneGain+1))/122;
tonepos = 1000*(10.^(freqExp*freq));

% Make colAxis to range according to the tonepos variable 
%colAxis = (1:length(time_seg))./length(time_seg);
spkMat = spkMat(:,idxStart:idxEnd);

for i = 1:size(spkMat, 1) % Loop through each neuron
    spike_times = [];
    spike_times = find(spkMat(i,:)>0); % Find spike times for the current neuron
    spike_times1 = spike_times(spike_times<idxEndTrial);
    y = i * ones(size(spike_times1)); % Y-coordinates for spikes
    scatter(time_seg_trunc(spike_times1), y, 2,tonepos(spike_times1),'filled'); 

    spike_times2 = spike_times(spike_times>=idxEndTrial);
    y = i * ones(size(spike_times2)); % Y-coordinates for spikes
    scatter(time_seg_trunc(spike_times2), y, 2,[0.5 0.5 0.5],'filled'); 
end

plot(time_seg_trunc,pos(idxStart:idxEnd)/130*size(spkMat,1),'Color',[0.6 0.6 0.6],'LineWidth',1)
line([trialEnd trialEnd],[0 size(spkMat,1)],'Color','k','LineWidth',1)
xlabel('Time(s)')
xlim([time_seg(idxStart) time_seg(idxEnd)])
ylim([0 size(spkMat,1)])
caxis([2000 23000])
colormap(ax1,viridis)

end

function plotChangePoints(fighandle,numrow,numcol,plotloc,change_points,trial,tsStart,tsEnd, post_time,trialNum,posterior_goal)

[~,idxstart] = min(abs(post_time-tsStart));
if post_time(idxstart)<tsStart %Take the next index
    idxstart = idxstart+1;
end        
[~,idxend] = min(abs(post_time-tsEnd));
if post_time(idxend)>tsEnd %Take the previous index
    idxend = idxend-1;
end   

[~,decGoal] = max(posterior_goal(:,idxstart:idxend));

[~, goal_dec] = max(posterior_goal, [], 1);

curChanges = change_points{trial==(trialNum-1)};
subplot(numrow,numcol,plotloc,'Parent',fighandle)
for cr = 1:length(curChanges)
    if cr == length(curChanges)
        trialDecGoal = mode(decGoal(curChanges(cr)+1:end));
        line([post_time(idxstart+curChanges(cr)-1) post_time(idxend)],[trialDecGoal trialDecGoal],'Color','m')
    elseif cr == 1
        trialDecGoal = mode(decGoal(curChanges(cr)+1:curChanges(cr+1)+1));
        line([post_time(idxstart) post_time(idxstart+curChanges(cr+1)-1)],[trialDecGoal trialDecGoal],'Color','m')
    else
        trialDecGoal = mode(decGoal(curChanges(cr)+1:curChanges(cr+1)+1));
        line([post_time(idxstart+curChanges(cr)-1) post_time(idxstart+curChanges(cr+1)-1)],[trialDecGoal trialDecGoal],'Color','m')
    end
    line([post_time(idxstart+curChanges(cr)-1) post_time(idxstart+curChanges(cr)-1)],[0 6],'Color','r','LineWidth',1.5);  
end

end
