function plotExampleCellProbeTrials(rowloc, colloc,spikeData, tracking, behavTrials, firingMaps,probeMaps, fighandle,numrows,numcol,cellNum,groupStyle)


b = linspace(0,125,50);
a = linspace(2000,22000,50);

%colMap = cbrewer('seq','Blues',18);
%col = [colMap(5,:);colMap(8,:);colMap(10,:);colMap(13,:);colMap(16,:);colMap(18,:)];
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
  %  vy = tracking.position.vy;
  %  idx2 = vy>0.2 & idx;
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
end


%% Tone trials
dataMat1 = [];
dataMat2 = [];
dataMatTone1 = [];
dataMatTone2 = [];

for kk = 1:6

    idx = find(behavTrials.linTrial(1:(end-1))==0 & ...
    behavTrials.lickLoc(1:(end-1)) ==(kk-1) & behavTrials.probe(1:(end-1))==groupStyle);        
       
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
    plot(b,probeMaps.forward.rateMaps{cellNum}{kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off   
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(7*numcol) (rowloc-1)*numcol+colloc+1+(7*numcol)], 'Parent', fighandle);
    hold on        
    plot(b,probeMaps.forward.rateMaps{cellNum}{kk+6},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(11*numcol) (rowloc-1)*numcol+colloc+1+(11*numcol)], 'Parent', fighandle);
    hold on
    plot(a,probeMaps.tone.rateMaps{cellNum}{kk},'Color',col(kk,:),'LineWidth',1);      
    xlim([1000 25000])
    set(gca,'xtick',[])
    box off
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(12*numcol) (rowloc-1)*numcol+colloc+1+(12*numcol)], 'Parent', fighandle);
    hold on
    plot(a,probeMaps.tone.rateMaps{cellNum}{6+kk},'Color',col(kk,:),'LineWidth',1);      
    xlim([1000 25000])
    set(gca,'xtick',[])
    box off    
    
    dataMat1 = [dataMat1;probeMaps.forward.rateMaps{cellNum}{kk}];
    dataMat2 = [dataMat2;probeMaps.forward.rateMaps{cellNum}{6+kk}];
    
    %atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{7+kk},'linear');
    dataMatTone1 = [dataMatTone1;probeMaps.tone.rateMaps{cellNum}{kk}];
    
    %atm = fillmissing(errorMaps.tone.rateMaps{cellNum}{7+kk},'linear');
    dataMatTone2 = [dataMatTone2;probeMaps.tone.rateMaps{cellNum}{6+kk}];
end   


YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(8*numcol) (rowloc-1)*numcol+colloc+1+(8*numcol)], 'Parent', fighandle);
imagesc(1,b,nanmean(dataMat1,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat1,1))))

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(9*numcol) (rowloc-1)*numcol+colloc+1+(9*numcol)], 'Parent', fighandle);
imagesc(1,b,nanmean(dataMat2,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat2,1))))


BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(13*numcol) (rowloc-1)*numcol+colloc+1+(13*numcol)], 'Parent', fighandle);
imagesc(1,a,nanmean(dataMatTone1,1))
colormap(ax2, BuPu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMatTone1,1))))


BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(14*numcol) (rowloc-1)*numcol+colloc+1+(14*numcol)], 'Parent', fighandle);
imagesc(1,a,nanmean(dataMatTone2,1))
colormap(ax2, BuPu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMatTone2,1))))


end