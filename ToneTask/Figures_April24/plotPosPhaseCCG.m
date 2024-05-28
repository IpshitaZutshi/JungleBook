function plotPosPhaseCCG(cell1, cell2, cell3,numrows, numcol, rowloc, fighandle)

%% Load files
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*cell_metrics.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file(1).name);

file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);


tSp{1} = spikeData.pos{cell1}; tSp{2} = spikeData.pos{cell2}; 
pSp{1} = rad2deg(spikeData.phase{cell1}); pSp{2} = rad2deg(spikeData.phase{cell2}); 

if ~isempty(cell3)
    tSp{3} = spikeData.pos{cell3};
    pSp{3} = rad2deg(spikeData.phase{cell3});
end

%% Get phase
[idx] = InIntervals(spikes.times{cell1},[tracking.timestamps(1) tracking.timestamps(end)]); 
tsBehav1 = spikes.times{cell1}(idx);

[idx] = InIntervals(spikes.times{cell2},[tracking.timestamps(1) tracking.timestamps(end)]); 
tsBehav2 = spikes.times{cell2}(idx);

if ~isempty(cell3)
    [idx] = InIntervals(spikes.times{cell3},[tracking.timestamps(1) tracking.timestamps(end)]); 
    tsBehav3 = spikes.times{cell3}(idx);
end

bPos = linspace(0,125,50);
binSize = 0.005;

%% Now, for each trial type, plot the spike phase relationships
for ii = 1:6

    idxTrial = behavTrials.lickLoc == (ii-1) & behavTrials.linTrial == 0 & behavTrials.correct == 1;
    intervals = behavTrials.timestamps(idxTrial,:)-0.033;

    %Select spikes within the relevant interval
    bools1 = InIntervals(tsBehav1,intervals);
    bools2 = InIntervals(tsBehav2,intervals);

    bools_1 = bools1 & tSp{1}>5;% & tSp{1}<70 ;
    bools_2 = bools2 & tSp{2}>5;% & tSp{2}<70 ;
    
    if ~isempty(cell3)
        bools3 = InIntervals(tsBehav3,intervals);
        bools_3 = bools3 & tSp{3}>5;
    end

    subplot(numrows, numcol,((rowloc-1)*numcol)+ii,'Parent',fighandle)
    scat_dat1 = [tSp{1}(bools_1);tSp{1}(bools_1)];
    scat_dat2 = [pSp{1}(bools_1);pSp{1}(bools_1)+360];
    scatter(scat_dat1 ,scat_dat2,4,[0/255 0/255 0/255],'filled')
    hold on
    scat_dat1 = [tSp{2}(bools_2);tSp{2}(bools_2)];
    scat_dat2 = [pSp{2}(bools_2);pSp{2}(bools_2)+360];    
    scatter(scat_dat1,scat_dat2,4,[8/255 133/255 161/255],'filled')

    if ~isempty(cell3)
        scat_dat1 = [tSp{3}(bools_3);tSp{3}(bools_3)];
        scat_dat2 = [pSp{3}(bools_3);pSp{3}(bools_3)+360];    
        scatter(scat_dat1,scat_dat2,4,'m','filled')%[224/255 163/255 46/255]
    end

    title(strcat('Port ',num2str(ii)))
    xlim([0 125])
    ylim([0 720])

    % Also plot the mean firing rates for those cells
    subplot(numrows, numcol,((rowloc-1)*numcol)+numcol+ii,'Parent',fighandle)
    plot(bPos,firingMaps.forward.rateMaps{cell1}{ii+1},'Color',[0/255 0/255 0/255],'Linewidth',1.5)
    hold on
    plot(bPos,firingMaps.forward.rateMaps{cell2}{ii+1},'Color',[8/255 133/255 161/255],'Linewidth',1.5) 
    if ~isempty(cell3)
        plot(bPos,firingMaps.forward.rateMaps{cell3}{ii+1},'Color','m','Linewidth',1.5)%[224/255 163/255 46/255]   
    end
    xlim([0 125])
    FieldInfo1 = detectFields(firingMaps.forward.rateMaps{cell1}{ii+1},'minFieldSize',2,'maxFieldSize',35);
    FieldInfo2 = detectFields(firingMaps.forward.rateMaps{cell2}{ii+1},'minFieldSize',2,'maxFieldSize',35);
    ylim([0 40])

    if ~isempty(FieldInfo1)
        line([bPos(FieldInfo1(1,4)) bPos(FieldInfo1(1,4))],[0 FieldInfo1(1,1)],'Color',[0/255 0/255 0/255],'LineWidth',1)
    end
    if ~isempty(FieldInfo2)
        line([bPos(FieldInfo2(1,4)) bPos(FieldInfo2(1,4))],[0 FieldInfo2(1,1)],'Color',[8/255 133/255 161/255],'LineWidth',1)
    end
    if ~isempty(FieldInfo1) & ~isempty(FieldInfo2)
        title(strcat('Dist: ',num2str(bPos(FieldInfo2(1,4))-bPos(FieldInfo1(1,4)))))
    end
    box off
    
    if sum(bools1)>0 && sum(bools2)>0
        % Finally, plot the ccgs between the cells for spikes 
        data1{1} = tsBehav1(bools_1 & tSp{1}>5);
        data1{2} = tsBehav2(bools_2 & tSp{2}>5);
        [cor, lag] = CCG(data1,[],'binSize',binSize,'duration',0.5);
        [cor2, lag2, smooth] = SmoothCor(cor(:,1,2), lag, binSize);
        subplot(numrows, numcol,((rowloc-1)*numcol)+2*numcol+ii,'Parent',fighandle)
        bar(lag2,cor(:,1,2),'EdgeColor','none','FaceColor',[0.6 0.6 0.6])
        hold on
        plot(lag2,smooth,'Color','k','LineWidth',1.5)
        a = bandpass(smooth,[6 12],1/binSize);
        plot(lag,a,'Color','r','LineWidth',1.5)
        xlim([-0.2 0.2])
    end
end

end