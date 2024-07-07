function FigS8_ProbeTrials

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[90 90 1780 914]);

numrows = 3;
numcol = 5;

%% Rastermap example of a probe trial
sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sess)
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
[spkMat, timestamps, sortIdx, keepCells, toneCellLog, placeCellLog] = calculateRastermap('expPath',sess,'loadSpecific',true);
sortIdx_new = [sortIdx(61:end) sortIdx(1:60)];
data = zscore(spkMat,[],1);

trialNum =40;% Correct non-probe trial, lick 5
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
trialStart = behavTrials.timestamps(trialNum,1)-0.033;
twin = [trialStart trialEnd];
ax = subplot(numrows,numcol,1);
imagesc(timestamps,1:260,data(sortIdx_new, :), [0, 1.5]);
colormap(ax, flipud(gray))
xlim(twin)
title('Non-probe, port 5 correct')
set(gca,'YDir','normal')
hold on

%Find the time when the tone stops
freqExp = log10(25000/1000);
[~,idxStart] = min(abs(tracking.timestamps-trialStart));
[~,idxEnd] = min(abs(tracking.timestamps-trialEnd));
pos = tracking.position.y(idxStart:idxEnd)/122;
tonepos = 1000*(10.^(freqExp*pos));
idxProbe = find(tonepos>4000,1,'first');
tProbe = tracking.timestamps(idxStart+idxProbe);
line([tProbe tProbe],[1 260],'Color','r','LineWidth',1.5)
plot(tracking.timestamps,(tracking.position.y/122*250),'Color','b')

trialNum = 22;% correct probe trial, lick 5
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
trialStart = behavTrials.timestamps(trialNum,1)-0.033;
twin = [trialStart trialEnd];
ax = subplot(numrows,numcol,2);
imagesc(timestamps,1:260,data(sortIdx_new, :), [0, 1.5]);
colormap(ax, flipud(gray))
xlim(twin)
title('Probe, port 5 correct')
set(gca,'YDir','normal')
hold on

%Find the time when the tone stops
freqExp = log10(25000/1000);
[~,idxStart] = min(abs(tracking.timestamps-trialStart));
[~,idxEnd] = min(abs(tracking.timestamps-trialEnd));
pos = tracking.position.y(idxStart:idxEnd)/122;
tonepos = 1000*(10.^(freqExp*pos));
idxProbe = find(tonepos>4000,1,'first');
tProbe = tracking.timestamps(idxStart+idxProbe);
line([tProbe tProbe],[1 260],'Color','r','LineWidth',1.5)
plot(tracking.timestamps,(tracking.position.y/122*250),'Color','b')


trialNum = 54;% Incorrect probe trial, lick 5
trialEnd = behavTrials.timestamps(trialNum,2)-0.033;
trialStart = behavTrials.timestamps(trialNum,1)-0.033;
twin = [5790 trialEnd];
ax = subplot(numrows,numcol,3);
imagesc(timestamps,1:260,data(sortIdx_new, :), [0, 1.5]);
colormap(ax, flipud(gray))
xlim(twin)
title('Probe, port 5 incorrect')
set(gca,'YDir','normal')
hold on

%Find the time when the tone stops
freqExp = log10(25000/1000);
[~,idxStart] = min(abs(tracking.timestamps-trialStart));
[~,idxEnd] = min(abs(tracking.timestamps-trialEnd));
pos = (tracking.position.y(idxStart:idxEnd)*122/102.79)/122;
tonepos = 1000*(10.^(freqExp*pos));
idxProbe = find(tonepos>4000,1,'first');
tProbe = tracking.timestamps(idxStart+idxProbe);
line([tProbe tProbe],[1 260],'Color','r','LineWidth',1.5)
plot(tracking.timestamps,(tracking.position.y/122*250),'Color','b')

%% Assemblies during probe trials. 
%[psth_probe,psth_nonprobe,t] = AssembliesProbenoProbe;
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\probeTrialAssemblies.mat');

[~,maxPSTH]  = max(psth_nonprobe,[],2);
[~,sortidx] = sort(maxPSTH,'descend');

RdPu=cbrewer('seq', 'RdPu', 11);

ax1 = subplot(numrows,numcol, 4);
imagesc(t,1:size(psth_nonprobe,1),zscore(psth_nonprobe(sortidx,:),[],2))
set(gca,'YDir','normal')
title('non probe')
caxis([-1 5])
colormap(ax1, RdPu)
hold on
line([0 0],[1 size(psth_nonprobe,1)],'Color','k','LineWidth',1.5)
xlim([-2 0.5])


ax1 = subplot(numrows,numcol, 5);
imagesc(t,1:size(psth_probe,1),zscore(psth_probe(sortidx,:),[],2))
set(gca,'YDir','normal')
title('probe')
hold on
line([0 0],[1 size(psth_probe,1)],'Color','k','LineWidth',1.5)
colormap(ax1, RdPu)
caxis([-1 5])
xlim([-2 0.5])
colorbar


%% UMAP of probe trial
TRIAL_TYPE = [0:5];
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230710_sess25\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
A = 92.3082;
E = 180;

manifoldPlot_probe('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',2,'colloc',1,'probe',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
manifoldPlot_probe('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',1,'probe',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)


umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230628_sess17\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
A = 179.1082;
E = 90;
manifoldPlot_probe('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',2,'colloc',3,'probe',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
manifoldPlot_probe('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',3,'probe',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure8A_ProbeTrials.png'));
saveas(gcf,strcat(expPath,'SupFigure8A_ProbeTrials.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure8A_ProbeTrials.fig'));

%% Make a separate figure with examples of single cells

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[500 150 1050 775]);

numrows = 7;
numcol = 6;

col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)
plotTrialFieldsProbe(155,numrows, numcol, 1, 1, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230628_sess17';
cd(sessloc)
plotTrialFieldsProbe(270,numrows, numcol, 1, 3, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230628_sess17';
cd(sessloc)
plotTrialFieldsProbe(234,numrows, numcol, 1, 5, fig2,col)

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure8B_ProbeTrialsSingleCells.png'));
saveas(gcf,strcat(expPath,'SupFigure8B_ProbeTrialsSingleCells.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure8B_ProbeTrialsSingleCells.fig'));

end

function plotTrialFieldsProbe(cellNum,numrows,numcol,rowloc,colloc,fighandle,col)

file = dir('*.rateMapsAvg.cellinfo.mat');
load(file(1).name);
file = dir(['*.rateMapsAvgLickLocProbe.cellinfo.mat']);
probeMaps = load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);
file = dir(['*spikeData.cellinfo.mat']);
load(file.name);

b = linspace(0,125,50);

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
end


for jj = 1:2
    for kk = 1:6   
        if jj == 1
            idx = find(behavTrials.linTrial(1:(end-1))==0 & ...
            behavTrials.lickLoc(1:(end-1)) ==(kk-1) & behavTrials.probe(1:(end-1))==0);        
        elseif jj == 2
            idx = find(behavTrials.linTrial(1:(end-1))==0 & ...
            behavTrials.lickLoc(1:(end-1)) ==(kk-1) & behavTrials.probe(1:(end-1))==1); 
        else
            idx = find(behavTrials.linTrial(1:(end-1))==0 & ...
            behavTrials.toneGain(1:(end-1)) ==(kk-1) & behavTrials.probe(1:(end-1))==1); 
        end
        
        if isempty(idx)
            continue
        end
        subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+((kk-1)*numcol)+(jj-1)], 'Parent', fighandle);
        for ii = 1:length(idx)
            plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
            hold on
            posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
            scatter(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),5,'r','filled')
            xlim([0 122])
            ylim([0 6])
            axis off
            box off
        end
    
        subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(6*numcol)+(jj-1)], 'Parent', fighandle);
        hold on     
        plot(b,probeMaps.firingMaps.forward.rateMaps{cellNum}{kk+((jj-1)*6)},'Color',col(kk,:),'LineWidth',1.5);   
        xlim([0 122])
        set(gca,'xtick',[])
        box off   
    end
end   

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
        scatter(spikes.times{cellIdx(i)}(spike_times),y,7,[0.7 0.7 0.7],'filled')        
    end
end
ylim([0 length(cellIdx)])

end