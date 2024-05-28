function plotExampleCell(rowloc, colloc,fighandle,numrows,numcol,cellNum)


b = linspace(0,125,50);
c = linspace(0,1.05,50);

freqExp = log10(22000/1000);
for ii = 1:length(c)
    a(ii) = (1000*(10.^(freqExp*c(ii))));
end

file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file.name);

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
   v = tracking.position.v;
   idx2 = v>5 & idx;
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
end

linIdx = find(behavTrials.linTrial==1);
jumpLin = find(diff(linIdx)>1);  
if isempty(jumpLin)
    idx = linIdx;
else
    idx = linIdx(1:jumpLin);
end
   
%% First no tone trials
subplot(numrows, numcol, [(rowloc-1)*numcol+colloc (rowloc-1)*numcol+colloc+1], 'Parent', fighandle);
for ii = 1:length(idx)-1
    plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
    hold on
    posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
    scatter(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),2,'r','filled')
    xlim([0 122])
    ylim([0 6])
    box off
    axis off
end    

if isempty(jumpLin)
    idx = [];
else
    idx = linIdx(jumpLin+1:end);
end

%% Second no tone trials
subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(7*numcol) (rowloc-1)*numcol+colloc+1+(7*numcol)], 'Parent', fighandle);
for ii = 1:length(idx)-1
    plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
    hold on
    posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
    scatter(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),2,'r','filled')
    xlim([0 122])
    ylim([0 6])
    box off
    axis off
end

subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(8*numcol) (rowloc-1)*numcol+colloc+1+(8*numcol)], 'Parent', fighandle);
plot(b,firingMaps.forward.rateMaps{cellNum}{1},'Color',[0 0 0],'LineWidth',1);
hold on
subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(8*numcol) (rowloc-1)*numcol+colloc+1+(8*numcol)], 'Parent', fighandle);
plot(b,firingMaps.forward.rateMaps{cellNum}{28},'Color',[0.5 0.5 0.5],'LineWidth',1);
xlim([0 122])
set(gca,'xtick',[])
box off

%% Tone trials
dataMat = [];
dataMatTone = [];

for kk = 1:6
    idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & ...
        behavTrials.toneGain(1:(end-1)) ==(kk-1) & behavTrials.stim(1:(end-1))==0);
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(kk*numcol) (rowloc-1)*numcol+colloc+1+(kk*numcol)], 'Parent', fighandle);
    for ii = 1:length(idx)
        plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
        hold on
        posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
        scatter(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),2,'r','filled')
        xlim([0 122])
        ylim([0 6])
        axis off
        box off
    end
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(9*numcol) (rowloc-1)*numcol+colloc+1+(9*numcol)], 'Parent', fighandle);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(12*numcol) (rowloc-1)*numcol+colloc+1+(12*numcol)], 'Parent', fighandle);
    hold on
    plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);      
    xlim([1000 22000])
    set(gca,'xtick',[])
    box off   
    xscale log
    
    dataMat = [dataMat;firingMaps.forward.rateMaps{cellNum}{1+kk}];
    atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{1+kk},'linear');
    dataMatTone = [dataMatTone;atm];
end   


YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(10*numcol) (rowloc-1)*numcol+colloc+1+(10*numcol)], 'Parent', fighandle);
imagesc(1,b,nanmean(dataMat,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat,1))))

BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(13*numcol) (rowloc-1)*numcol+colloc+1+(13*numcol)], 'Parent', fighandle);
imagesc(a,1,nanmean(dataMatTone,1))
colormap(ax2, BuPu)
box off
set(gca,'ytick',[])
xscale log
xlabel(num2str(max(nanmean(dataMatTone,1))))

end