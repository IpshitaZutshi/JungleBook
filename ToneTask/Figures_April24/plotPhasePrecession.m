function plotPhasePrecession

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[100 40 1140 950]);

numrows =10;
numcol = 6;

col = [ 241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

%% Panel A1: Place cell phase precession

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
plotTrialFields(47,numrows, numcol, 1, 1, fig2)
plotPhasePlots(47,numrows, numcol, 2, 1, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotTrialFields(5,numrows, numcol, 1, 2, fig2)
plotPhasePlots(5,numrows, numcol, 2, 2, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)
plotTrialFields(155,numrows, numcol, 1, 3, fig2)
plotPhasePlots(155,numrows, numcol, 2, 3, fig2,col)


%% Panel A2: Tone cell phase precession
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotTrialFields(103,numrows, numcol, 1, 4, fig2)
plotPhasePlots(103,numrows, numcol, 2, 4, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
plotTrialFields(28,numrows, numcol, 1, 5, fig2)
plotPhasePlots(28,numrows, numcol, 2, 5, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
plotTrialFields(167,numrows, numcol, 1, 6, fig2)
plotPhasePlots(167,numrows, numcol, 2, 6, fig2,col)

%% Now do the same analysis as for place fields, but comparing phase precession slope and intercept
% First for place cells 
precessData = compilePhasePrecession('tonecell',false,'plotfig',false);

for ii = 1:5
    subplot(numrows,numcol,6*numcol+ii)
    histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
    hold on
    line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
    line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
    idxSig{ii} = precessData.sigP(:,ii)<0.05;
    fractSig(ii) = sum(idxSig{ii})./sum(~isnan(precessData.slope(:,ii)));
    title(strcat('Place Port', num2str(ii),' fract:',num2str(fractSig(ii))))
end

subplot(numrows,numcol,7*numcol+1)
PhasePrecessionStats.Place.slope_all = groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
   precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope all') 

subplot(numrows,numcol,7*numcol+2)
PhasePrecessionStats.Place.intercept_all = groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
   precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept all') 

subplot(numrows,numcol,7*numcol+3)
PhasePrecessionStats.Place.slope_sig = groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
   precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope sig') 

subplot(numrows,numcol,7*numcol+4)
PhasePrecessionStats.Place.intercept_sig = groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
   precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept sig') 

PhasePrecessionStats.Place.fractSig = fractSig;

%% Next for tone cells 

precessData = compilePhasePrecession('tonecell',true,'plotfig',false);
idxSig = [];
fractSig = [];

for ii = 1:5
    subplot(numrows,numcol,8*numcol+ii)
    histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
    hold on
    line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
    line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
    idxSig{ii} = precessData.sigP(:,ii)<0.05;
    fractSig(ii) = sum(idxSig{ii})./sum(~isnan(precessData.slope(:,ii)));
    title(strcat('Tone Port', num2str(ii),' fract:',num2str(fractSig(ii))))
end

subplot(numrows,numcol,9*numcol+1)
PhasePrecessionStats.Tone.slope_all = groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
   precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope all') 

subplot(numrows,numcol,9*numcol+2)
PhasePrecessionStats.Tone.intercept_all = groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
   precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept all') 

subplot(numrows,numcol,9*numcol+3)
PhasePrecessionStats.Tone.slope_sig = groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
   precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope sig') 

subplot(numrows,numcol,9*numcol+4)
PhasePrecessionStats.Tone.intercept_sig = groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
   precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept sig') 

PhasePrecessionStats.Tone.fractSig = fractSig;


%% Save figure and stats

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure3B_PhasePrecession.png'));
saveas(gcf,strcat(expPath,'SupFigure3B_PhasePrecession.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure3B_PhasePrecession.fig'));
save(strcat(expPath,'SupFigure3B_PhasePrecession.mat'),'PhasePrecessionStats'); 

end

function plotPhasePlots(cellNum, numrows, numcol, rowloc, colloc, fighandle, col)

    % Load necessary data files
    file = dir('*.rateMapsAvg.cellinfo.mat');
    load(file(1).name);
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file.name);
    file = dir(['*spikeData.cellinfo.mat']);
    load(file.name);
    file = dir(['*.spikes.cellinfo.mat']);
    load(file.name);

    % Extract spike position and phase data and running speed
    tSp = spikeData.pos{cellNum};
    pSp = rad2deg(spikeData.phase{cellNum}); 
    vSp = tracking.position.vy(spikeData.posIdx{cellNum});

    %% Get phase within tracking intervals
    [idx] = InIntervals(spikes.times{cellNum},[tracking.timestamps(1) tracking.timestamps(end)]); 
    tsBehav1 = spikes.times{cellNum}(idx);

    shift = 1;
    speedThresh = 5;

    %% Plot spike phase relationships for each trial type
    for ii = 2:6
        subplot(numrows, numcol, colloc + ((ii - 1) * numcol), 'Parent', fighandle);
        
        % Select relevant trials
        idxTrial = behavTrials.lickLoc == (ii - 1) & behavTrials.linTrial == 0 & behavTrials.correct == 1;
        intervals = behavTrials.timestamps(idxTrial,:) - 0.033;

        % Select spikes within trial intervals
        bools1 = InIntervals(tsBehav1, intervals);
        bools_1 = bools1 & vSp > speedThresh;

        % Plot scatter plot
        scat_dat1 = [tSp(bools_1); tSp(bools_1)];
        phasedat = pSp(bools_1);
        phasedeg = phasedat - 180;
        phasedat = deg2rad(phasedeg);
        phase = (phasedat + pi) / (2 * pi);        
        scat_dat2 = [phase; phase + shift];
        scatter(scat_dat1, scat_dat2, 2, 'k', 'filled');
        hold on;
        xlim([0 122]);
        set(gca, 'xtick', []);
        box off;   

        %% Get the phase precession
        % Define place fields
        field_Info = detectFields(firingMaps.forward.rateMaps{cellNum}{1 + ii}, 'minFieldSize', 2, 'maxRate', 5, 'maxFieldSize', 35,'percentRate',0.2);
        if  ~isempty(field_Info)
            [~, idxMax] = max(field_Info(:, 3) - field_Info(:, 2));
            fieldStart = field_Info(idxMax, 2) * 125 / 50;
            fieldEnd = field_Info(idxMax, 3) * 125 / 50;
            fieldSize = fieldEnd - fieldStart;
    
            % Select spikes within place fields
            bools_1 = bools1 & tSp > fieldStart & tSp < fieldEnd & vSp> speedThresh; % Skip early spikes that are not accurate 
            lindat = (tSp(bools_1) - fieldStart) ./ fieldSize;
            phidat = pSp(bools_1);
            phasedeg = phidat - 180;
            phidat = deg2rad(phasedeg);

            % Calculate circular-linear regression
            [slope, offset, ~, R, p] = rccc(lindat, phidat, [-3 3], 0.05); 
            offset = (offset+pi)/(2*pi);
            reg = offset + slope * lindat;

            % Plot regression line
            hold on;        
            plot([tSp(bools_1) tSp(bools_1)], [reg reg+shift], 'r-', 'LineWidth', 1.5);            
            ylim([0 2])
            title(strcat('Slope: ', num2str(slope), ' intercept: ', num2str(offset)))
        end

    end     
end

function plotTrialFields(cellNum,numrows, numcol, rowloc, colloc, fighandle)

file = dir('*.rateMapsAvg.cellinfo.mat');
load(file(1).name);
file = dir(['*.rateMapsTrial.cellinfo.mat']);
trialMap = load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);

YlGnBu=cbrewer('seq', 'YlGnBu', 11);

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

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
end
