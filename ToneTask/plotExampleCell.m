function plotExampleCell(rowloc, colloc,spikeData, tracking, behavTrials, firingMaps,fighandle,numrows,numcol,cellNum)


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
    xlim([1000 25000])
    set(gca,'xtick',[])
    box off
    
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
imagesc(1,a,nanmean(dataMatTone,1))
colormap(ax2, BuPu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMatTone,1))))
end