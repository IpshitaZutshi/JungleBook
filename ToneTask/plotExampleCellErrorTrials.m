function plotExampleCellErrorTrials(rowloc, colloc,spikeData, tracking, behavTrials, firingMaps,errorMaps, fighandle,numrows,numcol,cellNum,groupStyle)


b = linspace(0,125,50);
c = linspace(0,1.05,50);

freqExp = log10(22000/1000);
for ii = 1:length(c)
    a(ii) = (1000*(10.^(freqExp*c(ii))));
end

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
end


%% Tone trials
dataMat1 = [];
dataMat2 = [];
dataMatTone1 = [];
dataMatTone2 = [];

for kk = 1:6
    if groupStyle==0
        idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==0 & ...
        behavTrials.lickLoc(1:(end-1)) ==(kk-1) & behavTrials.stim(1:(end-1))==0);        
    else
        idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==0 & ...                
        behavTrials.toneGain(1:(end-1)) ==(kk-1) & behavTrials.stim(1:(end-1))==0);
    end
        
    if isempty(idx)
        continue
    end
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+((kk-1)*numcol) (rowloc-1)*numcol+colloc+1+((kk-1)*numcol)], 'Parent', fighandle);
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
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(6*numcol) (rowloc-1)*numcol+colloc+1+(6*numcol)], 'Parent', fighandle);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{7+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off   
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(7*numcol) (rowloc-1)*numcol+colloc+1+(7*numcol)], 'Parent', fighandle);
    hold on        
    plot(b,errorMaps.forward.rateMaps{cellNum}{7+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(11*numcol) (rowloc-1)*numcol+colloc+1+(11*numcol)], 'Parent', fighandle);
    hold on
    plot(a,firingMaps.tone.rateMaps{cellNum}{7+kk},'Color',col(kk,:),'LineWidth',1);      
    xlim([1000 25000])
    xscale log
    set(gca,'xtick',[])
    box off
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(12*numcol) (rowloc-1)*numcol+colloc+1+(12*numcol)], 'Parent', fighandle);
    hold on
    plot(a,errorMaps.tone.rateMaps{cellNum}{7+kk},'Color',col(kk,:),'LineWidth',1);      
    xlim([1000 25000])
    xscale log
    set(gca,'xtick',[])
    box off    
    
    dataMat1 = [dataMat1;firingMaps.forward.rateMaps{cellNum}{7+kk}];
    dataMat2 = [dataMat2;errorMaps.forward.rateMaps{cellNum}{7+kk}];
    
    %atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{7+kk},'linear');
    dataMatTone1 = [dataMatTone1;firingMaps.tone.rateMaps{cellNum}{7+kk}];
    
    %atm = fillmissing(errorMaps.tone.rateMaps{cellNum}{7+kk},'linear');
    dataMatTone2 = [dataMatTone2;errorMaps.tone.rateMaps{cellNum}{7+kk}];
end   


YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(8*numcol) (rowloc-1)*numcol+colloc+1+(8*numcol)], 'Parent', fighandle);
imagesc(b,1,nanmean(dataMat1,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat1,1))))

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(9*numcol) (rowloc-1)*numcol+colloc+1+(9*numcol)], 'Parent', fighandle);
imagesc(b,1,nanmean(dataMat2,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat2,1))))


BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(13*numcol) (rowloc-1)*numcol+colloc+1+(13*numcol)], 'Parent', fighandle);
imagesc(a,1,nanmean(dataMatTone1,1))
colormap(ax2, BuPu)
box off
set(gca,'ytick',[])
xscale log
xlim([1000 25000])
xlabel(num2str(max(nanmean(dataMatTone1,1))))

caxis([0 35])

BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(14*numcol) (rowloc-1)*numcol+colloc+1+(14*numcol)], 'Parent', fighandle);
imagesc(a,1,nanmean(dataMatTone2,1))
colormap(ax2, BuPu)
box off
set(gca,'ytick',[])
xscale log
xlim([1000 25000])
xlabel(num2str(max(nanmean(dataMatTone2,1))))
caxis([0 35])

end