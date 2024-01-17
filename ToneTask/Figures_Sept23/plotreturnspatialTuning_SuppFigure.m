function plotreturnspatialTuning_SuppFigure

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 3;
numcol = 8;

%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);

%% Panel C: Place cells heat map - ACgN
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
linPos = linspace(1,122,50);
    
ss = 2;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
idxMaps{1} = idxSess & Summary.AllretFieldlin;
idxMaps{2} = idxSess & Summary.AllretFieldCorrect;
idxMaps{3} = idxSess & Summary.AllretFieldlinEnd;
         
for ii = 1:3            
    selectedlinMap = Summary.AllretMapLinInit(idxMaps{ii},:);
    selectedlinMapEnd = Summary.AllretMapLinEnd(idxMaps{ii},:);
    selectedMapCorrect = Summary.AllretMapCorr(idxMaps{ii},:);
    selectedMapInCorr = Summary.AllretMapIncorr(idxMaps{ii},:);
    
    [maxLin,idxLin] = max(selectedlinMap,[],2);
    [maxSpace,idxSpace] = max(selectedMapCorrect,[],2);    
    [maxLinEnd,idxLinEnd] = max(selectedlinMapEnd,[],2);    
    
    normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
    normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
    normspaceMap = (selectedMapCorrect-nanmean(selectedMapCorrect,2))./nanstd(selectedMapCorrect,[],2);
    normspaceMapInCorr = (selectedMapInCorr-nanmean(selectedMapInCorr,2))./nanstd(selectedMapInCorr,[],2);
      
    if ii ==1            
        [~,sortidx] = sort(idxLin,'ascend');
    elseif ii == 2            
        [~,sortidx] = sort(idxSpace,'ascend');
    elseif ii == 3
        [~,sortidx] = sort(idxLinEnd,'ascend');
    end
         
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+1);
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    caxis([-1 3])
    if ii==1
        title('No tone I')
    end
    
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+2);  
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii==1
        title('Tone Correct')
    end

    ax1 = subplot(numrows,numcol,numcol*(ii-1)+3);  
    imagesc(linPos, 1:length(sortidx),normspaceMapInCorr(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii==1
        title('Tone InCorrect')
    end
    
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+4);  
    imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel('Cell ID')
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii == 1
        title('No tone II')
    end
end   
  
% Control mice
ss = 1;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
idxMaps{1} = idxSess & Summary.AllretFieldlin;
idxMaps{2} = idxSess & Summary.AllretFieldCorrect;
idxMaps{3} = idxSess & Summary.AllretFieldlinEnd;
         
for ii = 1:3            
    selectedlinMap = Summary.AllretMapLinInit(idxMaps{ii},:);
    selectedlinMapEnd = Summary.AllretMapLinEnd(idxMaps{ii},:);
    selectedMapCorrect = Summary.AllretMapCorr(idxMaps{ii},:);
    selectedMapInCorr = Summary.AllretMapIncorr(idxMaps{ii},:);
    
    [maxLin,idxLin] = max(selectedlinMap,[],2);
    [maxSpace,idxSpace] = max(selectedMapCorrect,[],2);    
    [maxLinEnd,idxLinEnd] = max(selectedlinMapEnd,[],2);    
    
    normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
    normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
    normspaceMap = (selectedMapCorrect-nanmean(selectedMapCorrect,2))./nanstd(selectedMapCorrect,[],2);
    normspaceMapInCorr = (selectedMapInCorr-nanmean(selectedMapInCorr,2))./nanstd(selectedMapInCorr,[],2);
      
    if ii ==1            
        [~,sortidx] = sort(idxLin,'ascend');
    elseif ii == 2            
        [~,sortidx] = sort(idxSpace,'ascend');
    elseif ii == 3
        [~,sortidx] = sort(idxLinEnd,'ascend');
    end
         
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+5);
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    caxis([-1 3])
    if ii==1
        title('No tone I')
    end
    
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+6);  
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii==1
        title('Tone Correct')
    end

    ax1 = subplot(numrows,numcol,numcol*(ii-1)+7);  
    imagesc(linPos, 1:length(sortidx),normspaceMapInCorr(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii==1
        title('Tone InCorrect')
    end
    
    ax1 = subplot(numrows,numcol,numcol*(ii-1)+8);  
    imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel('Cell ID')
    xlabel('Position')
    axis off
    caxis([-1 3])
    if ii == 1
        title('No tone II')
    end
end   

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure_returnRun.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure_returnRun.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure_returnRun.fig'));

end