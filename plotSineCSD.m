function plotSineCSD

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

[pulses] = bz_getAnalogPulsesSine('analogCh',65);

freqRange = [50, 10, 20, 0];%unique(pulses.freq);

figure
set(gcf,'Position',[100 100 1500 900])
for ff = 1:length(freqRange)

    for jj = 1%:(size(sessionInfo.AnatGrps,2)-1)
        freqidx = pulses.freq == freqRange(ff);
        events = pulses.intsPeriods(1,freqidx);
        lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:numel(sessionInfo.AnatGrps(jj).Channels)-mod(numel(sessionInfo.AnatGrps(jj).Channels),8)),'noPrompts', true);
        twin1 = 5;
        twin2 = 5;
        [csd,lfpAvg] = bz_eventCSD(lfp,events','twin',[twin1 twin2],'plotLFP',false,'plotCSD',false);
        taxis = linspace(-twin1,twin2,size(csd.data,1));
        cmax = max(max(csd.data)); 
        subplot(1,4,ff);
        %subplot((size(sessionInfo.AnatGrps,2)-1),1,jj);
        contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
        set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('STIMULATION, Shank #',num2str(jj)),'FontWeight','normal'); 
        colormap jet; try caxis([-100 100]); end
        title(['Stim Frequency = ' num2str(freqRange(ff)) ' Hz'])
        hold on
        for kk = 1:size(lfpAvg.data,2)
            plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k','LineWidth',0.6)
        end
    end
end


end