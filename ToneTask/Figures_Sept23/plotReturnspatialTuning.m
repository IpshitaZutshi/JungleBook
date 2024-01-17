function plotReturnspatialTuning

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 18;
numcol = 9;

%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);

%% Panel C: Place cells heat map - ACgN
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
linPos = linspace(1,122,50);
    
ss = 2;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
idxMaps{1} = idxSess & Summary.AlllinField;% & Summary.AlllinCorr>0.1;
idxMaps{2} = idxSess & Summary.AllspaceField;% & Summary.AllspatialCorr>0.1;
         
for ii = 1:2            
    selectedlinMap = Summary.AlllinMapInit(idxMaps{ii},:);
    selectedlinMapEnd = Summary.AlllinMapEnd(idxMaps{ii},:);
    selectedspaceMap = Summary.AllspaceMap(idxMaps{ii},:);
    
    [maxLin,idxLin] = max(selectedlinMap,[],2);
    [maxSpace,idxSpace] = max(selectedspaceMap,[],2);    
    
    normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
    normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
    normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
      
    if ii ==1            
        [~,sortidx] = sort(idxLin,'ascend');
    elseif ii == 2            
        [~,sortidx] = sort(idxSpace,'ascend');
    end
         
    ax1 = subplot(numrows,numcol,[6+(numcol*(3*(ii-1))) 6+numcol+(numcol*(3*(ii-1))) 6+2*numcol+(numcol*(3*(ii-1)))]);    
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    caxis([-1 3.5])
    if ii==1
        title('No tone I')
    end
    
    ax1 = subplot(numrows,numcol,[7+(numcol*(3*(ii-1))) 7+numcol+(numcol*(3*(ii-1))) 7+2*numcol+(numcol*(3*(ii-1)))]);    
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3.5])
    if ii==1
        title('Tone')
    end
    
    ax1 = subplot(numrows,numcol,[8+(numcol*(3*(ii-1))) 8+numcol+(numcol*(3*(ii-1))) 8+2*numcol+(numcol*(3*(ii-1)))]);    
    imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel('Cell ID')
    xlabel('Position')
    axis off
    caxis([-1 3.5])
    if ii == 1
        title('No tone II')
    end
end   
    
%% Panel E: Correlation forward run
subplot(numrows,numcol,[numcol*15+1 numcol*15+2 numcol*15+3 numcol*15+4 numcol*16+1 ...
    numcol*16+2 numcol*16+3 numcol*16+4  numcol*17+1 numcol*17+2 numcol*17+3 numcol*17+4])
col = [52/243 52/243 52/243;...
    180/243 180/243 180/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243];

% Task
idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 ...
    & (Summary.AlllinField | Summary.AllspaceField | Summary.AlllinEndField);
%idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlllinField==1);%1 == 1 | Summary.AlllinField2 == 1);
data{1} = Summary.AlllinCorr(idxSess);

%idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlllinField == 1 | Summary.AllspaceField == 1);
data{3} = Summary.AlltoneNoToneCorr(idxSess);

%idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlllinEndField == 1 | Summary.AllspaceField == 1);
data{4} = Summary.AlltonelinEndCorr(idxSess);

%idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlllinEndField == 1 | Summary.AlllinField == 1);
data{2} = Summary.AlllinlinEndCorr(idxSess);

% Set up a one way ANOVA
datacombined = [data{1};data{2};data{3};data{4}];
group1 = [ones(length(data{1}),1);ones(length(data{2}),1)*2;ones(length(data{3}),1)*3;ones(length(data{4}),1)*4];

Fig6Stats.FwdCorr = groupStats(datacombined,[group1],'inAxis',true,'color',col,'plotType','boxplot','labelSummary',false);
ylim([-1 2])

%% Panel F: Return runs
%Task
idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 ...
    & (Summary.AllretFieldlin | Summary.AllretFieldCorrect | Summary.AllretFieldlinEnd);

data{1} = Summary.AllretlinlinCorr(idxSess);
data{2} = Summary.AllretlinlinEndCorr(idxSess);
data{3} = Summary.AllretlinToneCorr(idxSess);
data{4} = Summary.AllretlinEndToneCorr(idxSess);
%

% Set up a one way ANOVA
datacombined = [data{1};data{2};data{3};data{4}];
group1 = [ones(length(data{1}),1);ones(length(data{2}),1)*2;ones(length(data{3}),1)*3;ones(length(data{4}),1)*4];

subplot(numrows,numcol,[numcol*15+5 numcol*15+6 numcol*15+7 numcol*15+8 ...
    numcol*16+5 numcol*16+6 numcol*16+7 numcol*16+8 numcol*17+5 numcol*17+6 numcol*17+7 numcol*17+8])

Fig6Stats.ReturnCorr = groupStats(datacombined,group1,'inAxis',true,'color',col,'plotType','boxplot','labelSummary',false);
ylim([-1 2])

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure6a_spatialTuning.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure6a_spatialTuning.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure6a_spatialTuning.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\Figure6a_spatialTuning.mat'),'Fig6Stats'); 

end