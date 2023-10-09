function plotExampleCellAvgMaps(rowloc, colloc,firingMaps,fighandle,numrows,numcol,cellNum)


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

subplot(numrows, numcol, [(rowloc-1)*numcol+colloc (rowloc-1)*numcol+colloc+1], 'Parent', fighandle);
plot(b,firingMaps.forward.rateMaps{cellNum}{1},'Color',[0 0 0],'LineWidth',1);
hold on
subplot(numrows, numcol, [(rowloc-1)*numcol+colloc (rowloc-1)*numcol+colloc+1], 'Parent', fighandle);
plot(b,firingMaps.forward.rateMaps{cellNum}{28},'Color',[0.5 0.5 0.5],'LineWidth',1);
xlim([0 122])
set(gca,'xtick',[])
box off

%% Tone trials
dataMat = [];
dataMatTone = [];

for kk = 1:6
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+numcol (rowloc-1)*numcol+colloc+1+numcol], 'Parent', fighandle);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    
    
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+3*numcol (rowloc-1)*numcol+colloc+1+3*numcol], 'Parent', fighandle);
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
ax1 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+2*numcol (rowloc-1)*numcol+colloc+1+2*numcol], 'Parent', fighandle);
imagesc(1,b,nanmean(dataMat,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat,1))))

BuPu=cbrewer('seq', 'BuPu', 11);
ax2 = subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+4*numcol (rowloc-1)*numcol+colloc+1+4*numcol], 'Parent', fighandle);
imagesc(1,a,nanmean(dataMatTone,1))
colormap(ax2, BuPu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMatTone,1))))
end