function [spkMat, timestamps, sortIdx, keepCells, toneCellLog, placeCellLog] = calculateRastermap(varargin)

p = inputParser;
addParameter(p,'expPath',pwd);
addParameter(p,'plotfig',false);
addParameter(p,'loadSpecific',false);

parse(p,varargin{:});
expPath = p.Results.expPath;
plotfig = p.Results.plotfig;
loadSpecific = p.Results.loadSpecific;

cd(expPath)
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
sessionInfo = bz_getSessionInfo(pwd,'noPrompts', true);

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
logVector2 = rate>0.2;
keepCells = logicalVector'& logVector2;
spkMat = spkMat(keepCells,:);
cellId = find(keepCells);

if loadSpecific
    clear spkMat
    filename = 'IZ48_230714_sess28.rastermapData';
    %filename = 'IZ47_230710_sess25.rastermapData';
    directory_path  = 'C:\Data\Rastermap\preProcessedData\';
    load(strcat(directory_path,filename))
end

data = zscore(spkMat,[],1);

dataNdArray = py.numpy.array(data);

pyrun("from rastermap import Rastermap") %load interpreter, import main function
rmModel = pyrun("model = Rastermap(n_clusters=50, n_PCs=200, locality=0, time_lag_window=3, grid_upsample = 0).fit(spks)", "model", spks=dataNdArray);
sortIdx = int16(py.memoryview(rmModel.isort.data)) + 1; %back to MATLAB array, 1-indexing

toneCellLog = zeros(spikes.numcells,1);
placeCellLog = zeros(spikes.numcells,1);

for ii = 1:length(keepCells)
    if keepCells(ii) == 1        
        % Check if its a tone/place cell
        dataMatTone = [];
        dataMat = [];
    
        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{ii}{jj}]; % Place
            a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end  
    
        spaceMap = nanmean(dataMat,1);
        toneMap = nanmean(dataMatTone,1);        
    
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
        Field_Info_tone = detectFields(toneMap,'maxFieldSize',40);
        Field_Info_space = detectFields(spaceMap);
    
        if ~isempty(Field_Info_tone) && (toneCorr > 0.1)
            toneCellLog(ii) = 1;       
        end
    
        if ~isempty(Field_Info_space) && (spaceCorr > 0.1)
            placeCellLog(ii) = 1;
        end   
    
        if placeCellLog(ii)==1 && toneCellLog(ii) ==1
            if toneCorr > spaceCorr
                placeCellLog(ii) = 0;
            else
                toneCellLog(ii) = 0;
            end
        end
    end
end

if plotfig
    figure
    imagesc(data(sortIdx , :), [0, 1.5]);
    colormap(flipud(gray))
end

end