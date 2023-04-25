function plotStimCA1CortexSingleCells

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};
ctxCh = {{19 54},{19 54},{19 54},{27 49},{8 56},{0 54},{31 57},{27 40},{30 49},{5 53},{4 55},{19 58},{4 37 104},{4 37 104},{0 50 99},{40},{63 104},{48 99}};
analogCh = 2;
 
YlGnBu=cbrewer('div', 'RdYlBu', 100);
YlGnBu = YlGnBu(end:-1:1,:);
YlGnBu(YlGnBu>1) = 1;

PSTHCtxUD = [];
PSTHCtxStim = [];
PSTHCA1UD = [];
PSTHCA1Stim = [];

for i=1:size(pathToSessionsAll,2)

    cd(strcat(expPath,pathToSessionsAll{i}))
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

    load([sessionInfo.FileName '.spikes.cellinfo.mat'])
    load([sessionInfo.FileName '.ripples.events.mat'])

    ripTS = ripples.peaks;

    load([sessionInfo.FileName '.SlowWaves.events.mat'])
    load([sessionInfo.FileName '.pulses.events.mat']);

    if exist([sessionInfo.FileName '.SleepState.states.mat'],'file') 
       load([sessionInfo.FileName '.SleepState.states.mat']);
    else
       disp('No sleep states associated with this session');
    end

    if i< 16
        pulTr = (pulses.stimComb==analogCh);
    else
        pulTr = (pulses.stimPerID'==1 & (pulses.stimComb==2 | pulses.stimComb==3)); % for the CA3 mice, both 2 and 3 are the same
    end

    % Only take pulses that happened during NREM sleep
    events = pulses.intsPeriods(1,pulTr)';         
    stStim = events(InIntervals(events,SleepState.ints.NREMstate));
      
    % UP-DOWN transitions
    st = SlowWaves.ints.DOWN(:,1);
    
    % Identify channels
    Ctxchans = [];
    for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
        idx = find(sessionInfo.AnatGrps(jj).Channels==ctxCh{i}{jj});
        Ctxchans = [Ctxchans sessionInfo.AnatGrps(jj).Channels(1:idx)];
    end
        
    CA1chans = [];
    for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
        idx = find(sessionInfo.AnatGrps(jj).Channels==ctxCh{i}{jj});
        CA1chans = [CA1chans sessionInfo.AnatGrps(jj).Channels(idx+1:end)];
    end
        
    % Combine spike timing - Only select CA1 cells
    for jj = 1:spikes.numcells
        if ~isempty(find(CA1chans == spikes.maxWaveformCh(jj)))
            [stccg, t] = CCG({spikes.times{jj} stStim},[],'binSize',0.01,'duration',2,'norm','rate');
            PSTHCA1Stim = [PSTHCA1Stim; stccg(:,2,1)'];
            
            [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.01,'duration',2,'norm','rate');
            PSTHCA1UD = [PSTHCA1UD; stccg(:,2,1)'];
            
        elseif ~isempty(find(Ctxchans == spikes.maxWaveformCh(jj)))
            
            [stccg, t] = CCG({spikes.times{jj} stStim},[],'binSize',0.01,'duration',2,'norm','rate');
            PSTHCtxStim = [PSTHCtxStim; stccg(:,2,1)'];
            
            [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.01,'duration',2,'norm','rate');
            PSTHCtxUD = [PSTHCtxUD; stccg(:,2,1)'];            
        end
    end
    
end

figure
set(gcf,'Color','w')
set(gcf, 'Renderer','painters')
set(gcf,'Position',[1 41 1920 970])

subplot(3,2,1)
imagesc(t,1:1:size(PSTHCtxUD,1),PSTHCtxUD)
colormap(YlGnBu)
colorbar    
caxis([0 5])
title('Cortex UD')

subplot(3,2,3)
imagesc(t,1:1:size(PSTHCtxStim,1),PSTHCtxStim)
colormap(YlGnBu)
colorbar    
caxis([0 5])
title('Cortex Stim')

subplot(3,2,5)
meanpsth = nanmean(PSTHCtxUD,1);
stdpsth = nanstd(PSTHCtxUD,1)./sqrt(size(PSTHCtxUD,1));               
hold on
fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);

meanpsth = nanmean(PSTHCtxStim,1);
stdpsth = nanstd(PSTHCtxStim,1)./sqrt(size(PSTHCtxStim,1));               
hold on
fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[175/243 54/243 60/243],'linestyle','none','FaceAlpha',0.5);                    
hi = line(t,meanpsth,'LineWidth',1.5,'Color',[175/243 54/243 60/243]);
title('PSTH Spiking')

subplot(3,2,2)
imagesc(t,1:1:size(PSTHCA1UD,1),PSTHCA1UD)
colormap(YlGnBu)
colorbar    
caxis([0 5])
title('CA1 UD')

subplot(3,2,4)
imagesc(t,1:1:size(PSTHCA1Stim,1),PSTHCA1Stim)
colormap(YlGnBu)
colorbar    
caxis([0 5])
title('CA1 Stim')

subplot(3,2,6)
meanpsth = nanmean(PSTHCA1UD,1);
stdpsth = nanstd(PSTHCA1UD,1)./sqrt(size(PSTHCA1UD,1));               
hold on
fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);

meanpsth = nanmean(PSTHCA1Stim,1);
stdpsth = nanstd(PSTHCA1Stim,1)./sqrt(size(PSTHCA1Stim,1));               
hold on
fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[175/243 54/243 60/243],'linestyle','none','FaceAlpha',0.5);                    
hi = line(t,meanpsth,'LineWidth',1.5,'Color',[175/243 54/243 60/243]);
title('PSTH Spiking')


saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1SingleCells.png'));
saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1SingleCells.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1SingleCells.fig'));

end