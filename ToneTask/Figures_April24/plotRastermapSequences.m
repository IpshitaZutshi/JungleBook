function plotRastermapSequences

%% Go to the directory
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';

cd(sessloc)
plotLFP = 1;

if plotLFP
    numrows = 5;
    numcol = 4;
else
    numrows = 3;
    numcol = 4;
end

%% Load files
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

%% Generate original spike matrix
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

%% Now generate rastermapData
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

%% Now identify place cells an tone cells
toneCellLog = zeros(sum(keepCells),1);

count = 0;
for ii = 1:length(keepCells)
    if keepCells(ii) == 1

        count = count+1;
        %% Check if its a tone cell
        dataMatTone = [];
        for jj = 2:7    
            a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end                 
        toneMap = nanmean(dataMatTone,1);        

        [~,idxMax] = max(toneMap);    
        corrTone = []; 
        for pp = 1:6
           for jj = (pp+1):6           
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];
           end
        end     
        toneCorr = nanmean(corrTone);  

        %% Detect fields
        Field_Info = detectFields(toneMap);
        if ~isempty(Field_Info) && (toneCorr > 0.1) && idxMax>30
            toneCellLog(count) = 1;
        end
    end
end

%% Now plot
f1  = figure;
set(f1,'Renderer','painters')
set(f1,'Color','w')
set(f1,'Position',[1 41 1920 970])

subplot(numrows,numcol,1:4)
% imagesc(data(sortIdx_new , :), [0, 1.5]);
% colormap(flipud(gray))
% set(gca,'YDir','normal')
% grid off

% The sorting is the same but using the output from python because the
% grayscale looks better
rastermapDir = 'C:\Data\Rastermap\AnalyzedData';
A = imread(strcat(rastermapDir,'\Allrastermap.png'));
image(A(150:1350,220:2200,:))
axis off

%% Extract data for the UMAP
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ48_230714_sess28\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{8},'.position_behavior_speed_1_smooth_5.mat');
A = 23;
E = 1.55;
TRIAL_TYPE = [0 1 2 3 4 5];

% %% Plot a no tone trial
% trialNum =7;
% raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
% time_seg = timestamps(eventVariables.trialNum==trialNum);
% pos = constVariables.y(eventVariables.trialNum==trialNum);
% trialEnd = eventVariables.trialEnd(eventVariables.trialNum==trialNum);
% plotRaster(raster_seg,time_seg,toneCellLog,3,7,8,f1,pos,trialEnd)
% title('No Tone')
% 
% manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
%     'numrow',3,'numcol',7,'rowloc',3,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',behavTrials.timestamps(trialNum,:))

%% Plot a Tone port 6 trial
trialNum = 55;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,toneCellLog,numrows,numcol,5,f1,pos,trialEnd,5798.5,5802.3,behavTrials.toneGain(trialNum))
title('Port 6')
colorbar
%xlim([5798.5 5802.3])

manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5798.5 trialEnd])
subplot(numrows,numcol,9)

% %% Plot a Tone port 5 trial
% trialNum = 61;
% raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
% time_seg = timestamps(eventVariables.trialNum==trialNum);
% pos = constVariables.y(eventVariables.trialNum==trialNum);
% trialEnd = eventVariables.trialEnd(eventVariables.trialNum==trialNum);
% plotRaster(raster_seg,time_seg,toneCellLog,3,7,10,f1,pos,trialEnd)
% title('Port 5')
% 
% manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
%     'numrow',3,'numcol',7,'rowloc',3,'colloc',3,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',behavTrials.timestamps(trialNum,:))

%% Plot a Tone port 4 trial
trialNum = 60;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,toneCellLog,numrows,numcol,6,f1,pos,trialEnd,5860.5,5863.5,behavTrials.toneGain(trialNum))
title('Port 4')
colorbar

manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',2,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5860.5 trialEnd])
subplot(numrows,numcol,10)

% %% Plot a Tone port 3 trial
% trialNum = 25;
% raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
% time_seg = timestamps(eventVariables.trialNum==trialNum);
% pos = constVariables.y(eventVariables.trialNum==trialNum);
% trialEnd = eventVariables.trialEnd(eventVariables.trialNum==trialNum);
% plotRaster(raster_seg,time_seg,toneCellLog,3,7,12,f1,pos,trialEnd)
% title('Port 3')
% 
% manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
%     'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',5,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',behavTrials.timestamps(trialNum,:))


%% Plot a Tone port 2 trial
trialNum = 34;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,toneCellLog,numrows,numcol,7,f1,pos,trialEnd,5551.5,5554.5,behavTrials.toneGain(trialNum))
title('Port 2')
colorbar

manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
     'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',3,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5551.5 trialEnd])
subplot(numrows,numcol,11)

%% Plot a Tone port 1 trial
trialNum = 42;
raster_seg = spkMat(sortIdx_new,eventVariables.trialNum==trialNum);
time_seg = timestamps(eventVariables.trialNum==trialNum);
pos = constVariables.y(eventVariables.trialNum==trialNum);
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
plotRaster(raster_seg,time_seg,toneCellLog,numrows,numcol,8,f1,pos,trialEnd,5652,5655, behavTrials.toneGain(trialNum))
title('Port 1')
colorbar

manifoldPlot('figHandle',f1,'umap_path',umap_path,'umap_name','behavior_speed_1_smooth_10','behav_file',behav_file,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',4,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'singleTrial',true,'tsWin',[5652 trialEnd])
subplot(numrows,numcol,12)


%% Also plot the LFP and slow gamma/theta ratio
if plotLFP
    channelNum = 7;
    lfp = bz_GetLFP(channelNum,'noPrompts',true);
    plotSpectrum(55,lfp,behavTrials,numrows,numcol,4,1,f1,[5798.5 5802.3])
    plotSpectrum(60,lfp,behavTrials,numrows,numcol,4,2,f1,[5860.5 5863.5])
    plotSpectrum(34,lfp,behavTrials,numrows,numcol,4,3,f1,[5551.5 5554.5])
    plotSpectrum(42,lfp,behavTrials,numrows,numcol,4,4,f1,[5652 5655])
end


%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\RastermaptoUMAP.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\RastermaptoUMAP.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\RastermaptoUMAP.fig'));

end

function plotRaster(spkMat,time_seg,toneCellLog,numrows,numcol,plotloc,fighandle,pos,trialEnd,timeStart,timeEnd,toneGain)

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
    scatter(time_seg_trunc(spike_times1), y, 4,tonepos(spike_times1),'filled'); 

    spike_times2 = spike_times(spike_times>=idxEndTrial);
    y = i * ones(size(spike_times2)); % Y-coordinates for spikes
    scatter(time_seg_trunc(spike_times2), y, 4,[0.5 0.5 0.5],'filled'); 

    % if toneCellLog(i)==1
    %     %plot(time_seg(spike_times), y, '.','Color', [1 0 1], 'MarkerSize', 4); % Plot spikes
    %     scatter(time_seg(spike_times), y, 4,colAxis(spike_times),'filled'); % Plot spikes
    % else
    %     %scatter(time_seg(spike_times), y,'.', 'Color', [0.1 0.1 0.1], 'MarkerSize', 4); % Plot spikes
    %     scatter(time_seg(spike_times), y, 4,colAxis(spike_times),'filled'); 
    % end
end

plot(time_seg_trunc,pos(idxStart:idxEnd)/130*size(spkMat,1),'Color',[0.6 0.6 0.6],'LineWidth',1)
line([trialEnd trialEnd],[0 size(spkMat,1)],'Color','k','LineWidth',1)
xlabel('Time(s)')
xlim([time_seg(idxStart) time_seg(idxEnd)])
ylim([0 size(spkMat,1)])
caxis([2000 23000])
colormap(ax1,plasma)
end

function plotSpectrum(trialNum,lfp,behavTrials,numrows,numcol,rowloc,colloc,fighandle,win)

[~,idxStart] = min(abs(lfp.timestamps-win(1)));
[~,idxEnd] = min(abs(lfp.timestamps-win(2)));

ax1 = subplot(numrows,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle)
nfreqs = 100;
ncyc = 10;
fBand = [2 150];
wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',win);
lfpSpect = log10(abs(wavelet.data));
imagesc(wavelet.timestamps,wavelet.freqs,lfpSpect');
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
clim([2 4.5])
colormap(ax1,"jet")
hold on
line([behavTrials.timestamps(trialNum,2) behavTrials.timestamps(trialNum,2)],[2 150]);
xlim(win)
ylabel('Frequency (log)')

%% Ratio of the power in the 18-35 Hz band and theta band
gamIdx = wavelet.freqs>15 & wavelet.freqs <30;
thetaIdx = wavelet.freqs>6 & wavelet.freqs <12;

powGam = nanmedian(lfpSpect(:,gamIdx),2);
powTheta = nanmedian(lfpSpect(:,thetaIdx),2);

subplot(numrows,numcol,numcol*(rowloc-1)+colloc+numcol,'Parent',fighandle)
plot(wavelet.timestamps,powGam./powTheta)
xlim(win)
hold on
line([behavTrials.timestamps(trialNum,2) behavTrials.timestamps(trialNum,2)],[0.5 1]);
line([0 win(2)],[1 1])
ylabel('Gamma/Theta ratio')
end