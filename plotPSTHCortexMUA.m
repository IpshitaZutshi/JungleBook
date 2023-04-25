function plotPSTHCortexMUA

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cond = 1;
force = 0;
if cond == 0
    pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
    'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};
    ctxCh = {{19 54},{19 54},{19 54},{27 49},{8 56},{0 54},{31 57},{27 40},{30 49},{5 53},{4 55},{19 58},{4 37 104},{4 37 104},{0 50 99},{40},{63 104},{48 99}};
    analogCh = 2;
elseif cond == 1
    pathToSessionsAll = {'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3'};
    analogCh = [1 2 3];
    ctxCh = {{4 55},{19 58},{4 37 104},{4 37 104},{0 50 99}};
    
end
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);
downUP = 0;

if exist(strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\PSTHData',num2str(cond),'.mat'),'file') && ~force
    load(strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\PSTHData',num2str(cond),'.mat'))
else
    for ii = analogCh
        distPre{ii} = zeros(length(pathToSessionsAll),121);
    end

    for i=1:size(pathToSessionsAll,2)

        cd(strcat(expPath,pathToSessionsAll{i}))
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

        load([sessionInfo.FileName '.spikes.cellinfo.mat'])
        load([sessionInfo.FileName '.ripples.events.mat'])
        load([sessionInfo.FileName '.SleepState.states.mat']);

        ripTS = ripples.peaks;

        load([sessionInfo.FileName '.SlowWaves.events.mat'])
        % Load pulses
        disp('Getting analog-in inputs...');
        [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);


        % Identify channels
        CTXchans = [];
        for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
            idx = find(sessionInfo.AnatGrps(jj).Channels==ctxCh{i}{jj});
            CTXchans = [CTXchans sessionInfo.AnatGrps(jj).Channels(1:idx)];
        end

        if downUP ==0
            st = SlowWaves.ints.DOWN(:,1);
        elseif downUP == 1
            st = SlowWaves.ints.DOWN(:,2);
        end

        % Combine spike timing - Only select cortical cells
        MUAts = [];
        for jj = 1:spikes.numcells
            if ~isempty(find(CTXchans == spikes.maxWaveformCh(jj)))
                MUAts = [MUAts; spikes.times{jj}];
            end
        end

        %% PSTHs during stim
        % Identify ripples during stimulation
        for aa = analogCh
            keepIdx = [];
            if cond ==1
                if exist('pulses')
                    if aa<=2
                        pulTr = (pulses.stimComb==aa);
                    else
                        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==aa);
                    end
                end   
            elseif cond==0
                if i< 16
                    if aa<=2
                        pulTr = (pulses.stimComb==aa);
                    else
                        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==aa);
                    end
                else
                    pulTr = (pulses.stimPerID'==1 & (pulses.stimComb==2 | pulses.stimComb==3)); % for the CA3 mice, both 2 and 3 are the same
                end
            end
            events = pulses.intsPeriods(1,pulTr)';  

            % Only keep events that overlap with NREM sleep
            keepIdx = InIntervals(events,SleepState.ints.NREMstate);

            keepIdx = logical(keepIdx);

            if sum(keepIdx)>2 
                [stccg, t] = CCG({MUAts events(keepIdx)},[],'binSize',0.1,'duration',12,'norm','rate');            
                distPre{aa}(i,:) = stccg(:,2,1);
            else
                distPre{aa}(i,:) = nan;     
            end
        end
    end
    
    save(['Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Ripples\DownState\PSTHData',num2str(cond),'.mat'], ...
            'distPre','t');
end


figure
set(gcf,'Color','w')
set(gcf, 'Renderer','painters')
set(gcf,'Position',[1 41 1920 970])

for ii = analogCh
    subplot(2,3,ii)
    aa = ~isnan(distPre{ii}(:,1));
    RipPre = distPre{ii}(aa,:);
    imagesc(t,1:1:size(distPre{ii},1),zscore(RipPre,[],2))
    colormap(YlGnBu)
    colorbar    
    %caxis([-2 2])
    
    subplot(2,3,ii+3)
    meanpsth = nanmean(zscore(distPre{ii},[],2),1);
    stdpsth = nanstd(zscore(distPre{ii},[],2),1)./sqrt(size(distPre{ii},1));               
    hold on
    fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0 0 0],'linestyle','none','FaceAlpha',0.5);                    
    hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0 0 0]);    


%     subplot(5,4,ii+17) % Plot Stats
%     idxTS = t<-0.02 & t>-0.08;
%     temp = zscore(distRipPre{ii},[],2);
%     data{1} = nanmean(temp(:,idxTS),2);
%     temp = zscore(distRipStim{ii},[],2);
%     data{2} = nanmean(temp(:,idxTS),2);
%     temp = zscore(distRipPost{ii},[],2);
%     data{3} = nanmean(temp(:,idxTS),2);
%     stats{ii} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','inAxis',true);
    
end
% 
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\PSTHSpiking_',num2str(cond),'.png'));
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\PSTHSpiking_',num2str(cond),'.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\PSTHSpiking_',num2str(cond),'.fig'));
% save(strcat(expPath,'Compiled\Ripples\DownState\UDSpiking_',num2str(cond),'.mat'),'stats');  

end