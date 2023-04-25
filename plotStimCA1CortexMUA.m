function plotStimCA1CortexMUA

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cond = 0;
saveMat = 1;
force = 0;
if cond == 0
    pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
    'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};%'IZ26\IZ26_0um_201003_sess1',{4 37 104},
    ctxCh = {{19 54},{19 54},{19 54},{27 49},{8 56},{0 54},{31 57},{27 40},{30 49},{5 53},{4 55},{19 58},{4 37 104},{0 50 99},{40},{63 104},{48 99}};
    analogCh = 2;
elseif cond == 1
    pathToSessionsAll = {'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3'};
    analogCh = 1;
    ctxCh = {{4 55},{19 58},{4 37 104},{4 37 104},{0 50 99}};
    
end
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);
downUP = 0;

if exist(strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\StimCA1UnitData',num2str(cond),'.mat'),'file') && ~force
    load(strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\StimCA1UnitData',num2str(cond),'.mat'))
else
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
        st = events(InIntervals(events,SleepState.ints.NREMstate));
        
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
        MUAts = [];
        numMUAts = 0;
        for jj = 1:spikes.numcells
            if ~isempty(find(CA1chans == spikes.maxWaveformCh(jj)))
                MUAts = [MUAts; spikes.times{jj}];
                numMUAts = numMUAts+1;
            end
        end

        MUActx = [];
        numctxMUA = 0;
        for jj = 1:spikes.numcells
            if ~isempty(find(Ctxchans == spikes.maxWaveformCh(jj)))
                MUActx = [MUActx; spikes.times{jj}];
                numctxMUA = numctxMUA+1;
            end
        end
        
        %% PSTHs across the session

        % Firing rate of cortical cells around stimulation
        [stccg, tUD] = CCG({MUActx st},[],'binSize',0.02,'duration',15,'norm','rate');
        distCtxUnitUD(i,:) = stccg(:,2,1)./numctxMUA;
        
        % Firing rate of CA1 cells around stimulation
        [stccg, tUD] = CCG({MUAts st},[],'binSize',0.02,'duration',15,'norm','rate');
        distUnitUD(i,:) = stccg(:,2,1)./numMUAts;

    end

    if saveMat
        save(['Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\StimCA1UnitData',num2str(cond),'.mat'], ...
            'distCtxUnitUD','distUnitUD','tUD');
    end
end

figure
set(gcf,'Color','w')
set(gcf, 'Renderer','painters')
set(gcf,'Position',[1 41 1920 970])

subplot(3,1,1)
imagesc(tUD,1:1:size(distUnitUD,1),distUnitUD)
colormap(YlGnBu)
colorbar    
%caxis([-2 2])
xlim([-2 8])

subplot(3,1,2)
imagesc(tUD,1:1:size(distCtxUnitUD,1),distCtxUnitUD)
colormap(YlGnBu)
colorbar    
caxis([-2 2])
xlim([-2 8])

subplot(3,1,3)
meanpsth = nanmean(distUnitUD,1);
stdpsth = nanstd(distUnitUD,1)./sqrt(size(distUnitUD,1));               
hold on
fill([tUD; flipud(tUD)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(tUD,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);

meanpsth = nanmean(distCtxUnitUD,1);
stdpsth = nanstd(distCtxUnitUD,1)./sqrt(size(distCtxUnitUD,1));               
hold on
fill([tUD; flipud(tUD)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[175/243 54/243 60/243],'linestyle','none','FaceAlpha',0.5);                    
hi = line(tUD,meanpsth,'LineWidth',1.5,'Color',[175/243 54/243 60/243]);
title('Spiking around stimulation')
xlim([-2 8])

saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1Spiking_',num2str(cond),'.png'));
saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1Spiking_',num2str(cond),'.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Ripples\CortexCA1Comparison\StimCA1Spiking_',num2str(cond),'.fig'));

end