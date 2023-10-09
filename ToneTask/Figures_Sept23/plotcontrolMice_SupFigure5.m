function plotcontrolMice_SupFigure5

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrows = 19;
numcol = 9;

%% Panel A: Task schematic  - Control 
subplot(numrows,numcol,[1 2 numcol+1 numcol+2])
title('L. Control task schematic','FontName', 'Arial','FontSize', 9)
axis off

%% Panel B: Control Task trajectory example [15], N, Control licks heatmap[16]

subplot(numrows,numcol,[3 4 numcol+3 numcol+4])
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ46\Final\IZ46_230410_sess11';
cd(sessloc)
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);
plot(tracking.position.y,tracking.timestamps, 'Color',[0.5 0.5 0.5])
hold on
for tt  = 1:size(behavTrials.timestamps,1)
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,1)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'k','filled')
    
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,2)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'r','filled')
end
ylim([behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
set(gca,'YDir','reverse')
xlim([0 125])
box off

subplot(numrows,numcol,[5 6 numcol+5 numcol+6])
colMap = cbrewer('seq','Greys',3);
imagesc(behavTrials.numLicks(:,2:7))
colormap(colMap)
caxis([0 1])
xlim([0.5 7.5])
ylabel(strcat('Trials (',num2str(1),'-',num2str(size(behavTrials.timestamps,1)),')'))
%xticklabels({'Lick1','Lick2','Lick3','Lick4','Lick5','Lick6'})
h=gca; h.XAxis.TickLength = [0 0];
hold on
%colMap = cbrewer('seq','Blues',18);
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];
scatter(ones(1,sum(behavTrials.linTrial==1))*6.9,find(behavTrials.linTrial==1),1,[50/243 50/243 50/243],'filled');
scatter(ones(1,sum(behavTrials.toneGain==0&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==0&behavTrials.linTrial==0),1,col(1,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==1&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==1&behavTrials.linTrial==0),1,col(2,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==2&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==2&behavTrials.linTrial==0),1,col(3,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==3&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==3&behavTrials.linTrial==0),1,col(4,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==4&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==4&behavTrials.linTrial==0),1,col(5,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==5&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==5&behavTrials.linTrial==0),1,col(6,:),'filled');
box off


%% Panel C: Example spatial cell task % 13 rows, 2 columns each

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ41\Final\IZ41_220704_sess11';
cd(sessloc)
cellNum = 39;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,7, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);

%% Panel D: Place cells heat map - Control
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
linPos = linspace(1,122,50);
ss = 1;
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
         
    ax1 = subplot(numrows,numcol,[2*numcol+1+(numcol*(3*(ii-1))) 2*numcol+1+numcol+(numcol*(3*(ii-1))) 2*numcol+1+2*numcol+(numcol*(3*(ii-1)))]);    
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    caxis([-1 3.5])
    if ii == 1
        title('No tone I')
    end
    
    ax1 = subplot(numrows,numcol,[2*numcol+2+(numcol*(3*(ii-1))) 2*numcol+2+numcol+(numcol*(3*(ii-1))) 2*numcol+2+2*numcol+(numcol*(3*(ii-1)))]);    
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3.5])
    if ii == 1
        title('Tone')
    end
    
    ax1 = subplot(numrows,numcol,[2*numcol+3+(numcol*(3*(ii-1))) 2*numcol+3+numcol+(numcol*(3*(ii-1))) 2*numcol+3+2*numcol+(numcol*(3*(ii-1)))]);        
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
subplot(numrows,numcol,[numcol*9+1 numcol*9+2 numcol*9+3 numcol*9+4 numcol*10+1 ...
    numcol*10+2 numcol*10+3 numcol*10+4  numcol*11+1 numcol*11+2 numcol*11+3 numcol*11+4])
col = [52/243 52/243 52/243;...
    180/243 180/243 180/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243];

% Control
idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & ...
    (Summary.AlllinField == 1 | Summary.AllspaceField == 1| Summary.AlllinEndField == 1);
data{1} = Summary.AlllinCorr(idxSess);
data{3} = Summary.AlltoneNoToneCorr(idxSess);
data{4} = Summary.AlltonelinEndCorr(idxSess);
data{2} = Summary.AlllinlinEndCorr(idxSess);

% Set up a one way ANOVA
datacombined = [data{1};data{2};data{3};data{4}];
group1 = [ones(length(data{1}),1);ones(length(data{2}),1)*2;ones(length(data{3}),1)*3;ones(length(data{4}),1)*4];

SupFig5Stats.FwdCorr = groupStats(datacombined,group1,'inAxis',true,'color',col,'plotType','boxplot','labelSummary',false);

%% Panel F: Return runs
%Task
idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & ...
    (Summary.AllretFieldlin == 1 | Summary.AllretFieldCorrect == 1| Summary.AllretFieldlinEnd == 1);
data{1} = Summary.AllretlinlinCorr(idxSess);
data{3} = Summary.AllretlinToneCorr(idxSess);
data{4} = Summary.AllretlinEndToneCorr(idxSess);
data{2} = Summary.AllretlinlinEndCorr(idxSess);

% Set up a one way ANOVA
datacombined = [data{1};data{2};data{3};data{4}];
group1 = [ones(length(data{1}),1);ones(length(data{2}),1)*2;ones(length(data{3}),1)*3;ones(length(data{4}),1)*4];

subplot(numrows,numcol,[numcol*13+1 numcol*13+2 numcol*13+3 numcol*13+4 numcol*14+1 ...
    numcol*14+2 numcol*14+3 numcol*14+4  numcol*15+1 numcol*15+2 numcol*15+3 numcol*15+4])

SupFig5Stats.ReturnCorr = groupStats(datacombined,group1,'inAxis',true,'color',col,'plotType','boxplot','labelSummary',false);

%% Control versus Task
idxSess0 = Summary.AllsessType==0 & Summary.AllcellType == 1 & ...
    (Summary.AlllinField == 1 | Summary.AllspaceField == 1| Summary.AlllinEndField == 1);
idxSess1 = Summary.AllsessType==1 & Summary.AllcellType == 1 & ...
    (Summary.AlllinField == 1 | Summary.AllspaceField == 1| Summary.AlllinEndField == 1);

data{1} = Summary.AlllinCorr(idxSess0);
data{3} = Summary.AlltoneNoToneCorr(idxSess0);
data{4} = Summary.AlltonelinEndCorr(idxSess0);
data{2} = Summary.AlllinlinEndCorr(idxSess0);

data{5} = Summary.AlllinCorr(idxSess1);
data{7} = Summary.AlltoneNoToneCorr(idxSess1);
data{8} = Summary.AlltonelinEndCorr(idxSess1);
data{6} = Summary.AlllinlinEndCorr(idxSess1);

col = [52/243 52/243 52/243;...
    52/243 52/243 52/243;...
    180/243 180/243 180/243;...
    180/243 180/243 180/243;...    
    56/243 61/243 150/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243;...       
    80/243 91/243 166/243];
    
% Set up a two way ANOVA
datacombined = [data{1};data{2};data{3};data{4};data{5};data{6};data{7};data{8}];
group1 = [ones(length(data{1}),1);ones(length(data{2}),1)*2;ones(length(data{3}),1)*3;ones(length(data{4}),1)*4;...
    ones(length(data{5}),1);ones(length(data{6}),1)*2;ones(length(data{7}),1)*3;ones(length(data{8}),1)*4];
group2 = [ones(length(data{1}),1);ones(length(data{2}),1)*1;ones(length(data{3}),1)*1;ones(length(data{4}),1)*1;...
    ones(length(data{5}),1)*2;ones(length(data{6}),1)*2;ones(length(data{7}),1)*2;ones(length(data{8}),1)*2];

subplot(numrows,numcol,[numcol*16+1 numcol*16+2 numcol*16+3 numcol*16+4 numcol*17+1 ...
    numcol*17+2 numcol*17+3 numcol*17+4  numcol*18+1 numcol*18+2 numcol*18+3 numcol*18+4])

SupFig5Stats.CtrlvsTask = groupStats(datacombined,[group1 group2],'inAxis',true,'color',col,'plotType','boxplot','labelSummary',false);

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure5_controlMice.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure5_controlMice.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure5_controlMice.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\SupFigure5_controlMice.mat'),'SupFig5Stats'); 

end