function plotPlaceFieldExpansion

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[30 275 1800 550]);

numrows = 3;
numcol = 6;

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

%% Panel A1: Example cells that scale/ do not scale
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
%plotPosPhaseCCG(5, 83, [], numrows, numcol, 1, fig2)
plotTrialFields(79,numrows, numcol, 1, 1, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220919_sess14';
cd(sessloc)
plotTrialFields(21,numrows, numcol, 1, 2, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)
plotTrialFields(155,numrows, numcol, 1, 3, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220828_sess5';
cd(sessloc)
plotTrialFields(36,numrows, numcol, 1, 4, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ40\Final\IZ40_220714_sess18';
cd(sessloc)
plotTrialFields(9,numrows, numcol, 1, 5, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
plotTrialFields(47,numrows, numcol, 1, 6, fig2,col)

%% Compile place field data across cells

%Skip- the first port because it's too short
fieldData = compileFieldExpansion;
trialCOM = fieldData.trialCOM(:,2:end);
ind = sum(~isnan(trialCOM),2)==5;

subplot(numrows, numcol,13)    
Stats.PFPeak = groupStats({fieldData.placefield(ind,2),fieldData.placefield(ind,3),...
   fieldData.placefield(ind,4),fieldData.placefield(ind,5),fieldData.placefield(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Place field peak position') 

subplot(numrows, numcol,14)
Stats.COM = groupStats({fieldData.trialCOM(ind,2),fieldData.trialCOM(ind,3),...
   fieldData.trialCOM(ind,4),fieldData.trialCOM(ind,5),fieldData.trialCOM(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Trial COM') 

% subplot(numrows, numcol,14)
% subMat  = (fieldData.trialCOM(ind,:)-fieldData.avgCOM(ind));
% PlaceFieldExpansionStats.COMdiff = groupStats({subMat(:,2),subMat(:,3),subMat(:,4),subMat(:,5),subMat(:,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
% title('Trial COM - avg COM') 

subplot(numrows, numcol,15)
Stats.PFStart = groupStats({fieldData.fieldStart(ind,2),fieldData.fieldStart(ind,3),...
   fieldData.fieldStart(ind,4),fieldData.fieldStart(ind,5),fieldData.fieldStart(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Place field start') 

subplot(numrows, numcol,16)
Stats.StartRate = groupStats({fieldData.fieldStartRate(ind,2),fieldData.fieldStartRate(ind,3),...
   fieldData.fieldStartRate(ind,4),fieldData.fieldStartRate(ind,5),fieldData.fieldStartRate(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Place field start rate') 

subplot(numrows, numcol,17)
Stats.PFEnd = groupStats({fieldData.fieldEnd(ind,2),fieldData.fieldEnd(ind,3),...
   fieldData.fieldEnd(ind,4),fieldData.fieldEnd(ind,5),fieldData.fieldEnd(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Place field end') 

subplot(numrows, numcol,18)
Stats.EndRate = groupStats({fieldData.fieldEndRate(ind,2),fieldData.fieldEndRate(ind,3),...
   fieldData.fieldEndRate(ind,4),fieldData.fieldEndRate(ind,5),fieldData.fieldEndRate(ind,6)},2:6,'inAxis',true,'color',col(2:end,:),'labelSummary',false,'repeatedMeasures',true);
title('Place field end rate')    

%% Save figure and stats

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure3A_PlaceCellTruncation.png'));
saveas(gcf,strcat(expPath,'SupFigure3A_PlaceCellTruncation.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure3A_PlaceCellTruncation.fig'));
save(strcat(expPath,'SupFigure3A_PlaceCellTruncation.mat'),'Stats'); 


end

function plotTrialFields(cellNum,numrows, numcol, rowloc, colloc, fighandle,col)

file = dir('*.rateMapsAvg.cellinfo.mat');
load(file(1).name);
file = dir(['*.rateMapsTrial.cellinfo.mat']);
trialMap = load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);

YlGnBu=cbrewer('seq', 'YlGnBu', 11);

idxTrials = find(behavTrials.linTrial==0);

dataMat = [];
for kk = 1:(length(idxTrials)-1)
    dataMat(kk,:) = trialMap.firingMaps.forward.rateMaps{cellNum}{idxTrials(kk)};
end  
spaceMap = nanmean(dataMat,1);

bPos = linspace(0,125,50);
ax1 = subplot(numrows,numcol,(numcol*(rowloc-1))+colloc,'Parent',fighandle);

%% Detect fields
rm = [];
for kk = 1:6      
    hold on
    rm = [rm;firingMaps.forward.rateMaps{cellNum}{1+kk}];
    plot(bPos,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5)  
    xlim([0 122])
    set(gca,'xtick',[])
    box off
    field_Info = detectFields(firingMaps.forward.rateMaps{cellNum}{1+kk},'minFieldSize',2,'maxFieldSize',35,'maxRate',2);
    if ~isempty(field_Info)
        [~,idxMax] = max(field_Info(:,3)-field_Info(:,2)); % Take the location that has the largest field
        extractedField = firingMaps.forward.rateMaps{cellNum}{1+kk}(field_Info(idxMax,2):field_Info(idxMax,3));
        last_nonNaN_index = find(~isnan(extractedField), 1, 'last');
        extractedField = extractedField(1:last_nonNaN_index);
        trialCOM(kk) = (((1:length(extractedField))*extractedField')./sum(extractedField))+(field_Info(idxMax,2)-1);
        %[maxFR,idx] = max(firingMaps.forward.rateMaps{cellNum}{kk+1});

        %line([bPos(idx) bPos(idx)],[0 maxFR],'Color',col(kk,:),'LineWidth',1)
        %line([trialCOM(kk)*125/50 trialCOM(kk)*125/50],[0 field_Info(idxMax,1)],'Color',col(kk,:),'LineWidth',1.5)
        scatter(trialCOM(kk)*125/50,-1,30,col(kk,:),"filled")
    else
        trialCOM(kk) = nan;
    end
end

% h = imagesc(bPos,1:6,rm);
% set(h, 'AlphaData', ~isnan(rm))
% scatter(trialCOM*125/50,1:6,25,'magenta','filled')
%set(gca,'Ydir','normal')
%ylim([0.5 6.5])
% xlim([0 122])
% %set(gca,'xtick',[])
% colormap(ax1,YlGnBu)
% box off
% colorbar
       
% Plot average map
subplot(numrows,numcol,(numcol*(rowloc-1))+colloc+numcol,'Parent',fighandle);

Field_Info_space = detectFields(spaceMap,'minFieldSize',4,'maxFieldSize',35);
if ~isempty(Field_Info_space) 
    [~,idxMax] = max(Field_Info_space(:,3)-Field_Info_space(:,2)); % Take the location that has the largest field
    extractedField = spaceMap(Field_Info_space(idxMax,2):Field_Info_space(idxMax,3));
    trialCOM = (((1:length(extractedField))*extractedField')./sum(extractedField))+(Field_Info_space(idxMax,2)-1);

    plot(bPos, spaceMap,'Color',[0.2 0.2 0.2],'LineWidth',1.5)
    line([trialCOM*125/50 trialCOM*125/50],[0 Field_Info_space(idxMax,1)],'Color',[0.2 0.2 0.2],'LineWidth',1.5)
end
xlim([0 122])
%set(gca,'xtick',[])
box off
end