function plotThetaSequences

%% Go to the directory
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)

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

%% Now identify tone cells
toneCellLog = zeros(length(keepCells),1);
placeCellLog = zeros(spikes.numcells,1);

for ii = 1:length(keepCells)
    if keepCells(ii) == 1

        %% Check if its a tone cell
        dataMatTone = [];
        dataMat = [];

        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{ii}{jj}]; % Place
            a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end  

        spaceMap = nanmean(dataMat,1);
        toneMap = nanmean(dataMatTone,1);        

        [~,idxMax] = max(toneMap);    
        corrTone = []; 
        corrSpace = [];
        for pp = 1:6
           for jj = (pp+1):6      
               a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];                  
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];
           end
        end     
        toneCorr = nanmean(corrTone);  
        spaceCorr = nanmean(corrSpace);

        %% Detect fields
        Field_Info_tone = detectFields(toneMap);
        Field_Info_space = detectFields(spaceMap);
    
        if ~isempty(Field_Info_tone) && (toneCorr > 0.1) && (idxMax > 30)
            toneCellLog(ii) = 1;       
        end
    
        if ~isempty(Field_Info_space) && (spaceCorr > 0.1)
            placeCellLog(ii) = 1;
        end    
    end
end

%% Now load rastermapData
filename = 'IZ48_230714_sess28.rastermapData';
directory_path  = 'C:\Data\Rastermap\preProcessedData\';
load(strcat(directory_path,filename))
data = zscore(spkMat,[],1);

dataNdArray = py.numpy.array(data);

pyrun("from rastermap import Rastermap") %load interpreter, import main function
rmModel = pyrun("model = Rastermap(n_clusters=50, n_PCs=200, locality=0, time_lag_window=3, grid_upsample = 0).fit(spks)", "model", spks=dataNdArray);
sortIdx = int16(py.memoryview(rmModel.isort.data)) + 1; %back to MATLAB array, 1-indexing

sortIdx_new = [sortIdx(61:end) sortIdx(1:60)];
sortIdx = double(sortIdx_new);

%% Extract theta from channel 63, filter and estimate phase
lfp = bz_GetLFP(63, 'noprompts',true);
passband = [6 12];
[b, a] = butter(3,[passband(1)/(lfp.samplingRate/2) passband(2)/(lfp.samplingRate/2)],'bandpass'); % order 3
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
rawlfp = double(lfp.data(:,1));
power = fastrms(filt,ceil(lfp.samplingRate./passband(1)));  % approximate power is frequency band
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);

%% Start plotting
f1 = figure;
set(gcf,'Color','w')
set(gcf,'Renderer','painters')
set(gcf,"Position",[2100 50 1615 900])

%% Plot the raw data for each trial, first trial 60, which is a port 4 trial
trialNum = 60;
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
twin = [5860.5 5863.5];
[~,idxStart] = min(abs(lfp.timestamps-twin(1)));
[~,idxEnd] = min(abs(lfp.timestamps-twin(2)));
plotRastSpikes(spikes,sortIdx,twin,keepCells,toneCellLog,placeCellLog,6,1,1,f1)

% Rescale filtered LFP between 0 and 1
currLFP = filt(idxStart:idxEnd);
currLFP_scale = (currLFP+abs(min(currLFP)));
maxF = max(currLFP_scale);
subplot(6,1,1)
plot(lfp.timestamps(idxStart:idxEnd),(currLFP_scale./maxF)*length(sortIdx),'Color',[0.3 0.3 0.3])
line([trialEnd trialEnd],[0 length(sortIdx)],'Color','k','LineWidth',1.5)
line([5861.9 5861.9],[0 260],'Color','r','LineWidth',1.5)
xlim([5861.51 5863.51])

% Plot position below
[~,idxStart] = min(abs(tracking.timestamps-twin(1)));
[~,idxEnd] = min(abs(tracking.timestamps-twin(2)));
subplot(6,1,2)
yyaxis left
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.y(idxStart:idxEnd),'Color',[0.2 0.2 0.2])
hold on
ylabel('Position')
yyaxis right
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.v(idxStart:idxEnd))
ylabel('Speed')
xlim([5861.5 5863.5])
xlim([5861.51 5863.51])

%% Next trial
trialNum = 34;
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
twin = [5551.5 5554.5];
[~,idxStart] = min(abs(lfp.timestamps-twin(1)));
[~,idxEnd] = min(abs(lfp.timestamps-twin(2)));

plotRastSpikes(spikes,sortIdx,twin,keepCells,toneCellLog,placeCellLog,6,1,3,f1)

% Rescale filtered LFP between 0 and 1
currLFP = filt(idxStart:idxEnd);
currLFP_scale = (currLFP+abs(min(currLFP)));
maxF = max(currLFP_scale);
subplot(6,1,3)
plot(lfp.timestamps(idxStart:idxEnd),(currLFP_scale./maxF)*length(sortIdx),'Color',[0.3 0.3 0.3])
line([trialEnd trialEnd],[0 length(sortIdx)],'Color','k','LineWidth',1.5)
twin = [5551.5 5554.5];
line([5552.6 5552.6],[0 260],'Color','r','LineWidth',1.5)
xlim([5552.1 5554.1])

% Plot position below
[~,idxStart] = min(abs(tracking.timestamps-twin(1)));
[~,idxEnd] = min(abs(tracking.timestamps-twin(2)));
subplot(6,1,4)
yyaxis left
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.y(idxStart:idxEnd),'Color',[0.2 0.2 0.2])
hold on
ylabel('Position')
yyaxis right
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.v(idxStart:idxEnd))
ylabel('Speed')
xlim([5552.1 5554.1])

%% Next trial
trialNum = 42;
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
twin = [5652 5655];
[~,idxStart] = min(abs(lfp.timestamps-twin(1)));
[~,idxEnd] = min(abs(lfp.timestamps-twin(2)));

plotRastSpikes(spikes,sortIdx,twin,keepCells,toneCellLog,placeCellLog,6,1,5,f1)

% Rescale filtered LFP between 0 and 1
currLFP = filt(idxStart:idxEnd);
currLFP_scale = (currLFP+abs(min(currLFP)));
maxF = max(currLFP_scale);
subplot(6,1,5)
plot(lfp.timestamps(idxStart:idxEnd),(currLFP_scale./maxF)*length(sortIdx),'Color',[0.3 0.3 0.3])
line([trialEnd trialEnd],[0 length(sortIdx)],'Color','k','LineWidth',1.5)
twin = [5652 5655];
line([5653.42 5653.42],[0 260],'Color','r','LineWidth',1.5)
xlim([5652.6 5654.6])

% Plot position below
[~,idxStart] = min(abs(tracking.timestamps-twin(1)));
[~,idxEnd] = min(abs(tracking.timestamps-twin(2)));
subplot(6,1,6)
yyaxis left
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.y(idxStart:idxEnd),'Color',[0.2 0.2 0.2])
hold on
ylabel('Position')
yyaxis right
plot(tracking.timestamps(idxStart:idxEnd),tracking.position.v(idxStart:idxEnd))
ylabel('Speed')
xlim([5652.6 5654.6])

%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure4A_RastermapThetaAssemblies.png'));
saveas(gcf,strcat(expPath,'SupFigure4A_RastermapThetaAssemblies.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure4A_RastermapThetaAssemblies.fig'));
end

function plotRastSpikes(spikes,sortIdx,timeWin,keepCells,toneCellLog,placeCellLog,numrows,numcol,plotloc,fighandle)

ax1 = subplot(numrows, numcol, plotloc,'Parent',fighandle);
hold on;
cellIdx = find(keepCells);
for i =  1:length(cellIdx)% Loop through each neuron
    spike_times = find(spikes.times{cellIdx(i)}>timeWin(1) & spikes.times{cellIdx(i)}<timeWin(2)); % Find spike times for the current neuron        
    %Find its order
    posrast = find(sortIdx==i);
    y = posrast * ones(size(spike_times)); % Y-coordinates for spikes
    if toneCellLog(cellIdx(i))==1
        scatter(spikes.times{cellIdx(i)}(spike_times),y,7,[1 0 1],'filled')         
    elseif placeCellLog(cellIdx(i))==1
        scatter(spikes.times{cellIdx(i)}(spike_times),y,7,[0.1 0.1 0.1],'filled')           
    else
        scatter(spikes.times{cellIdx(i)}(spike_times),y,7,[0.9 0.9 0.9],'filled')        
    end
end
ylim([0 length(cellIdx)])

end