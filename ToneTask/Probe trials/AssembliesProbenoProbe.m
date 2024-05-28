function [psth_probe,psth_nonprobe,t] = AssembliesProbenoProbe

sess = {'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};    

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
plotfig = 0;
saveMat = 1;

psth_nonprobe = [];
psth_probe = [];

for ss = 1:length(sess)
    try
        %% Load files
        cd(strcat(expPath,sess{ss})) 
    
        file = dir(['*spikes.cellinfo.mat']);
        load(file(1).name);
        file = dir(['*cell_metrics.cellinfo.mat']);
        load(file.name);
        file = dir(['*rateMapsAvg.cellinfo.mat']);
        load(file(1).name);
        file = dir(['*TrialBehavior.Behavior.mat']);
        load(file(1).name);
        file = dir(['*Tracking.Behavior.mat']);
        load(file(1).name);
        [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
        
        dtime = mean(diff(tracking.timestamps));
        spkData = bz_SpktToSpkmat(spikes.times,'dt',dtime, 'win',[behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
        spikeMat = spkData.data';
        
        %Only select pyramidal cells
        logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
        keepCells = logicalVector;
        spikeMat = spikeMat(keepCells,:);
        cellId = find(keepCells);
        
        SpikeCount = zscore(spikeMat,[],1);
        AssemblyTemplates = assembly_patterns(SpikeCount);
        
        Vectors = [];
        assembliesID = [];
    
        for aa = 1:size(AssemblyTemplates,2)
            assembliesVector = AssemblyTemplates(:,aa);
            % flip if negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            Vectors(:,aa) = assembliesVector;
        
            %For each assembly, define significant cells
            sigCells = find(assembliesVector> 2.5*std(assembliesVector));
            assembliesID{aa} = cellId(sigCells);
        end  
        
        Activities = assembly_activity(Vectors,SpikeCount);

        % Calculate PSTH of assemblies for probe and non-rpobe trials
        % separately
        assemblyTimes = [];
        for aa = 1:size(Activities,1)
            idx = find(Activities(aa,:) > 1.5*std(Activities(aa,:)));
            assemblyTimes{aa} = spkData.timestamps(idx);
            [stccg, t] = CCG({assemblyTimes{aa} behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.probe==0,2)},[],'binSize',0.1,'duration',5,'norm','rate');
            psth_nonprobe = [psth_nonprobe; stccg(:,2,1)'];    

            [stccg, t] = CCG({assemblyTimes{aa} behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.probe==1,2)},[],'binSize',0.1,'duration',5,'norm','rate');
            psth_probe = [psth_probe; stccg(:,2,1)'];  
        end
    end 
end

if saveMat
    save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\probeTrialAssemblies.mat','psth_probe','psth_nonprobe','t')
end

if plotfig
    [~,maxPSTH]  = max(psth_nonprobe,[],2);
    [~,sortidx] = sort(maxPSTH,'descend');
    
    figure
    RdPu=cbrewer('seq', 'RdPu', 11);
    ax1 = subplot(1,2,1);
    imagesc(t,1:size(psth_nonprobe,1),zscore(psth_nonprobe(sortidx,:),[],2))
    set(gca,'YDir','normal')
    title('non probe')
    caxis([-2 5])
    colormap(ax1, RdPu)
    hold on
    line([0 0],[1 size(psth_nonprobe,1)],'Color','k','LineWidth',1.5)
    xlim([-2 0.5])
    
    ax1 = subplot(1,2,2);
    imagesc(t,1:size(psth_probe,1),zscore(psth_probe(sortidx,:),[],2))
    set(gca,'YDir','normal')
    title('probe')
    hold on
    line([0 0],[1 size(psth_probe,1)],'Color','k','LineWidth',1.5)
    colormap(ax1, RdPu)
    caxis([-2 5])
    xlim([-2 0.5])
end

end