function plotProbeTrials_SupFigure4

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[2000 27 788 950]);

BuPu=cbrewer('seq', 'BuPu', 11);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

numrows = 25;
numcol = 9;

%% Panel A1: Example tone tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
cd(sessloc)
cellNum = 86;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgLickLocProbe.cellinfo.mat']);
Maps = load(file.name);
probeMaps = Maps.firingMaps;
plotExampleCellProbeTrials(1,1, spikeData, tracking, behavTrials, firingMaps,probeMaps,fig2,numrows,numcol,cellNum,0)
plotExampleCellProbeTrials(1,3, spikeData, tracking, behavTrials, firingMaps,probeMaps,fig2,numrows,numcol,cellNum,1)

%% Panel A2: Example tone tuned cell 2

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
cd(sessloc)
cellNum = 95;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgLickLocProbe.cellinfo.mat']);
Maps = load(file.name);
probeMaps = Maps.firingMaps;
plotExampleCellProbeTrials(1,5, spikeData, tracking, behavTrials, firingMaps,probeMaps ,fig2,numrows,numcol,cellNum,0)
plotExampleCellProbeTrials(1,7, spikeData, tracking, behavTrials, firingMaps,probeMaps ,fig2,numrows,numcol,cellNum,1)


%% Calculate place cells
Summary = compileProbeTrials('plotfig',false,'savefig',false);
linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);

%% Panel B: Heatmaps of error and no error tuning
idxMaps{1} = Summary.AllcellType & Summary.AllspaceField & Summary.AllspatialCorr>0.1; % Idx for spatial cells
idxMaps{2} = Summary.AllcellType & Summary.AlltoneField & Summary.AlltoneCorr>0.1;  % Idx for tone cells 

for ii = 1:2
   selectedtoneMap = Summary.AlltoneMap(idxMaps{ii},:);
   selectedspaceMap = Summary.AllspaceMap(idxMaps{ii},:);
   selectedtoneMapProbe = Summary.AlltoneMapProbe(idxMaps{ii},:);
   selectedspaceMapProbe = Summary.AllspaceMapProbe(idxMaps{ii},:);
   
   [maxSpace,idxSpace] = max(selectedspaceMap,[],2);        
   [maxTone,idxTone] = max(selectedtoneMap,[],2);

       %Z score normalize
    normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
    normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
    normtoneMapProbe = (selectedtoneMapProbe-nanmean(selectedtoneMapProbe,2))./nanstd(selectedtoneMapProbe,[],2);
    normspaceMapProbe = (selectedspaceMapProbe-nanmean(selectedspaceMapProbe,2))./nanstd(selectedspaceMapProbe,[],2);

    sortidx = [];
    if ii == 1            
        [~,sortidx] = sort(idxSpace,'ascend');
    elseif ii == 2
        [~,sortidx] = sort(idxTone,'ascend');
    end

    ax1 = subplot(numrows,numcol,[(numcol*16)+1+(numcol*5*(ii-1)) (numcol*16)+2+(numcol*5*(ii-1)) (numcol*17)+1+(numcol*5*(ii-1)) (numcol*17)+2+(numcol*5*(ii-1)) ...
        (numcol*18)+1+(numcol*5*(ii-1)) (numcol*18)+2+(numcol*5*(ii-1)) (numcol*19)+1+(numcol*5*(ii-1)) (numcol*19)+2+(numcol*5*(ii-1))]);
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    caxis([-1 3])
    colormap(ax1, YlGnBu) 
    title('Space no probe')

    ax1 = subplot(numrows,numcol,[(numcol*16)+3+(numcol*5*(ii-1)) (numcol*16)+4+(numcol*5*(ii-1)) (numcol*17)+3+(numcol*5*(ii-1)) (numcol*17)+4+(numcol*5*(ii-1)) ...
        (numcol*18)+3+(numcol*5*(ii-1)) (numcol*18)+4+(numcol*5*(ii-1)) (numcol*19)+3+(numcol*5*(ii-1)) (numcol*19)+4+(numcol*5*(ii-1))]);
    imagesc(linPos, 1:length(sortidx),normspaceMapProbe(sortidx,:))
    caxis([-1 3])
    colormap(ax1, YlGnBu)
    title('Space Probe')  
        
    ax1 = subplot(numrows,numcol,[(numcol*16)+5+(numcol*5*(ii-1)) (numcol*16)+6+(numcol*5*(ii-1)) (numcol*17)+5+(numcol*5*(ii-1)) (numcol*17)+6+(numcol*5*(ii-1)) ...
        (numcol*18)+5+(numcol*5*(ii-1)) (numcol*18)+6+(numcol*5*(ii-1)) (numcol*19)+5+(numcol*5*(ii-1)) (numcol*19)+6+(numcol*5*(ii-1))]);
    imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
    caxis([-1 3])
    colormap(ax1, BuPu)
    title('Tone No Probe')  
    
    ax1 = subplot(numrows,numcol,[(numcol*16)+7+(numcol*5*(ii-1)) (numcol*16)+8+(numcol*5*(ii-1)) (numcol*17)+7+(numcol*5*(ii-1)) (numcol*17)+8+(numcol*5*(ii-1)) ...
        (numcol*18)+7+(numcol*5*(ii-1)) (numcol*18)+8+(numcol*5*(ii-1)) (numcol*19)+7+(numcol*5*(ii-1)) (numcol*19)+8+(numcol*5*(ii-1))]);
    imagesc(linTone, 1:length(sortidx),normtoneMapProbe(sortidx,:))
    caxis([-1 3])
    colormap(ax1, BuPu)
    title('Tone Probe') 
end

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure4_probeTrials.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure4_probeTrials.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure4_probeTrials.fig'));

end