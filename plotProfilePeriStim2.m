function plotProfilePeriStim2

analogEv = 64;
numAnalog = 2;
pre = 1.4;
post = 1.4;
fixChannels = 1;
%[colormap] = cbrewer('seq','PuBuGn',90);
%[colormap] = copper(32);
for ii = 1:numAnalog
    analogCh(ii) = (analogEv-1)+ii;
end

[pulses] = bz_getAnalogPulses('analogCh',analogCh);
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);


for jj = 1 %1:(size(sessionInfo.AnatGrps,2)-1)
    lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels,'noPrompts', true);
    if fixChannels
        lfp = bz_interpolateLFP(lfp);
    end
    
    [colormap] = cbrewer('seq','PuBuGn',length(sessionInfo.AnatGrps(jj).Channels)+40);
    colormap(colormap<0) = 0;
    data = lfp.data;
    timestamps = lfp.timestamps;
%         
%         [Lia] = ismember(sessionInfo.AnatGrps(ii).Channels, channels);
%         nC = 1:length(sessionInfo.AnatGrps(ii).Channels);
%         nC = nC(Lia)';
    for i = 2%:(numAnalog+1)
        if i<=numAnalog
            pulTr = (pulses.stimComb==i);
        else
            pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
        end
        events = pulses.intsPeriods(:,pulTr);
        events = round(events*1250);
        events = events(:,(events(1,:) + (15*1250) <= size(data,1)) & (events(1,:) - (5*1250) > 0));

        for pp = [36 38]%1:length(events(1,:))
            figure
            set(gcf,'renderer','Painters')
            hold on
            for kk = 1:(length(sessionInfo.AnatGrps(jj).Channels))
                plot(timestamps((events(1,pp)-(1250*pre)):(events(1,pp)+(1250*post))),0.75*(data((events(1,pp)-(1250*pre)):(events(1,pp)+(1250*post)),kk))-(kk-1)*400,'Color',colormap(kk+22,:))
            end
            title(strcat('AnalogCh = ',num2str(i)))
            xlim([timestamps((events(1,pp)-(1250*pre))) timestamps((events(1,pp)+(1250*post)))])
            line([timestamps(events(1,pp)) timestamps(events(1,pp))],[2500 (-2*10^4)],'Color','red')
            line([timestamps(events(2,pp)) timestamps(events(2,pp))],[2500 (-2*10^4)],'Color','red')
        end


    end
end