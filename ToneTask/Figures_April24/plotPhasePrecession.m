function plotPhasePrecession

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[30 275 1800 550]);

numrows =9;
numcol = 6;

col = [ 241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

%% Panel A1: Place cell phase precession

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
plotPhasePlots(47,numrows, numcol, 1, 1, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotPhasePlots(5,numrows, numcol, 1, 2, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)
plotPhasePlots(155,numrows, numcol, 1, 3, fig2,col)


%% Panel A2: Tone cell phase precession
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotPhasePlots(103,numrows, numcol, 1, 4, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
plotPhasePlots(28,numrows, numcol, 1, 5, fig2,col)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
plotPhasePlots(167,numrows, numcol, 1, 6, fig2,col)

%% Now do the same analysis as for place fields, but comparing phase precession slope and intercept
% First for place cells 
precessData = compilePhasePrecession('tonecell',false,'plotfig',false);

for ii = 1:5
    subplot(numrows,numcol,5*numcol+ii)
    histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
    hold on
    line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
    line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
    idxSig{ii} = precessData.sigP(:,ii)<0.05;
    fractSig(ii) = sum(idxSig{ii})./sum(~isnan(precessData.slope(:,ii)));
    title(strcat('Place Port', num2str(ii),' fract:',num2str(fractSig(ii))))
end

subplot(numrows,numcol,6*numcol+1)
PhasePrecessionStats.Place.slope_all = groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
   precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope all') 

subplot(numrows,numcol,6*numcol+2)
PhasePrecessionStats.Place.intercept_all = groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
   precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept all') 

subplot(numrows,numcol,6*numcol+3)
PhasePrecessionStats.Place.slope_sig = groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
   precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope sig') 

subplot(numrows,numcol,6*numcol+4)
PhasePrecessionStats.Place.intercept_sig = groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
   precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept sig') 

PhasePrecessionStats.Place.fractSig = fractSig;

%% Next for tone cells 

precessData = compilePhasePrecession('tonecell',true,'plotfig',false);
idxSig = [];
fractSig = [];

for ii = 1:5
    subplot(numrows,numcol,7*numcol+ii)
    histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
    hold on
    line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
    line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
    idxSig{ii} = precessData.sigP(:,ii)<0.05;
    fractSig(ii) = sum(idxSig{ii})./sum(~isnan(precessData.slope(:,ii)));
    title(strcat('Tone Port', num2str(ii),' fract:',num2str(fractSig(ii))))
end

subplot(numrows,numcol,8*numcol+1)
PhasePrecessionStats.Tone.slope_all = groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
   precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope all') 

subplot(numrows,numcol,8*numcol+2)
PhasePrecessionStats.Tone.intercept_all = groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
   precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept all') 

subplot(numrows,numcol,8*numcol+3)
PhasePrecessionStats.Tone.slope_sig = groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
   precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Slope sig') 

subplot(numrows,numcol,8*numcol+4)
PhasePrecessionStats.Tone.intercept_sig = groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
   precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col,'labelSummary',false);
title('Intercept sig') 

PhasePrecessionStats.Tone.fractSig = fractSig;


%% Save figure and stats

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\PhasePrecession.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\PhasePrecession.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_April2024\PhasePrecession.fig'));
save(strcat(expPath,'\Compiled\Figures_April2024\PhasePrecession.mat'),'PhasePrecessionStats'); 


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
        subplot(numrows, numcol, colloc + ((ii - 2) * numcol), 'Parent', fighandle);
        
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
        scatter(scat_dat1, scat_dat2, 5, 'k', 'filled');
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
