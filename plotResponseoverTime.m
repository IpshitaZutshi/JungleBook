function plotResponseoverTime


[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
load([sessionInfo.FileName '.MergePoints.events.mat'])
% 
% if exist([sessionInfo.FileName '.region.mat'],'file') 
%     load([sessionInfo.FileName '.region.mat']);
%     pyrCh =39;%region.CA1sp; 24 25 0 38 (CA1, s.r., s.l.m, DG?)25,26,19,30,4,61,63,56,38,43,41(IZ28) 39 21(sr) 43 17 (slm) 47 16 46 (CA3?) 19 45
% else
%     disp('First calculate .region file to identify ripple channel! Skipping');
%     return
% end
pyrCh =81;
analyzeLFP = 1;
analyzeSpikes = 0;
fixNoise = 0;

%% First plot LFP spectogram, and theta, high, mid, low gamma over time
if analyzeLFP
    figure
    lfp = bz_GetLFP(pyrCh,'noPrompts', true);
    if fixNoise
        lfpNoise = bz_GetLFP(29,'noPrompts', true);
        lfp.data = lfp.data-lfpNoise.data;
    end
    params.Fs = lfp.samplingRate; params.fpass = [2 300]; params.tapers = [3 5]; params.pad = 1;
    [S,t,f] = mtspecgramc(single(lfp.data),[2 1],params);
    S1 = log10(S); % in Db
    S_det= bsxfun(@minus,S1,polyval(polyfit(f,mean(S1,1),2),f)); % detrending

    subplot(5,1,1)
    imagesc(t,f,S_det',[-1.5 1.5]);
    set(gca,'YDir','normal','XTick',[]); ylabel('Freqs');
    hold on
    for kk = 1:length(MergePoints.timestamps(:,2))
        xline(MergePoints.timestamps(kk,2),'w','LineWidth',1)
    end

    subplot(5,1,2)
    fIdx = f>=6 & f <=12;
    plot(t,mean(S_det(:,fIdx),2));
    ylabel('Power')
    xlabel('Time')
    title('theta power')
    xlim([0 max(t)])

    subplot(5,1,3)
    fIdx = f>=20 & f <=45;
    plot(t,mean(S_det(:,fIdx),2));
    ylabel('Power')
    xlabel('Time')
    title('20-45 Hz power')
    xlim([0 max(t)])

    subplot(5,1,4)
    fIdx = f>=50 & f <=90;
    plot(t,mean(S_det(:,fIdx),2));
    ylabel('Power')
    xlabel('Time')
    title('50-90 Hz power')
    xlim([0 max(t)])

    subplot(5,1,5)
    fIdx = f>=95 & f <=150;
    plot(t,mean(S_det(:,fIdx),2));
    ylabel('Power')
    xlabel('Time')
    title('95 - 150 Hz power')
    xlim([0 max(t)])
    
end

%% Plot spike responses

if analyzeSpikes
    % Load spikes
    if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
       spikes = bz_LoadPhy('noPrompts',true);
    end
    [cell_metrics] = bz_CellMetricsSimpleIZ(spikes,sessionInfo.rates.wideband);

    spkRate = [];
    region = [];
    win = [300 150];
    
    spikemat = bz_SpktToSpkmat(spikes,'dt',10,'binsize',30,'units','rate');
    for jj = 1:size(spikes.UID,2)
        %spkRate(jj,:) = calculateSpkIFR(spikes.times{jj},MergePoints.timestamps(end,2),win);
        maxChan = spikes.maxWaveformCh(jj);
        idxChan = [];
        for aa = 1:(size(sessionInfo.AnatGrps,2)-1)
            if isempty(idxChan)
                idxChan = find(sessionInfo.AnatGrps(aa).Channels==maxChan);
            end
        end  
        if ~isempty(idxChan)
            spkChan(jj) = idxChan;
        else
            spkChan(jj) = nan;
        end
    end

    [~,sortIdx] = sort(spkChan);
    %sortIdx = cell_metrics.FR<5 & spkChan' >41;
    spkRateNorm = zscore(spikemat.data,[],1);

    timebins = 1:win(2):(MergePoints.timestamps(end,2)-(win(1)/2));
    
    figure
    subplot(3,1,1)
    imagesc(spikemat.timestamps,1:1:size(spkRateNorm(:,sortIdx),2),spkRateNorm(:,sortIdx)')
    hold on
    for kk = 1:length(MergePoints.timestamps(:,2))
        xline(MergePoints.timestamps(kk,2),'Color','w','LineWidth',1)
    end

    subplot(3,1,2)
    imagesc(spikemat.timestamps,1:1:sum(sortIdx),spikemat.data(:,sortIdx)')
    hold on
    for kk = 1:length(MergePoints.timestamps(:,2))
        xline(MergePoints.timestamps(kk,2),'Color','w','LineWidth',1)
    end
    
    subplot(3,1,3)
    plot(spikemat.timestamps, mean(spkRateNorm(:,sortIdx)',1));
    hold on
    for kk = 1:length(MergePoints.timestamps(:,2))
        xline(MergePoints.timestamps(kk,2),'LineWidth',1)
    end    
end

end

function spkRate = calculateSpkIFR(ts,t,win)

timebins = 1:win(2):(t-(win(1)/2));

for tt = 1:length(timebins)
    spkRate(tt) = sum(ts>=timebins(tt) & ts < (timebins(tt)+ (win(1)/2)))/ win(1);
end

end

