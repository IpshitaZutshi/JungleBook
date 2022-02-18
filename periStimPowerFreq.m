lfpCh = [51 57 39 37];%[21 17 16 19]; %[18 31 39 34]for IZ4, [59 63 53 37] for IZ5

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
frange = [1 250];
params.Fs = sessionInfo.rates.lfp; params.fpass = frange; 
params.tapers = [3 5]; params.pad = 1;

if ~exist('analogEv','var')
    analogEv = 64;
    numAnalog = 2;

    if ~isempty(analogEv)
        for ii = 1:numAnalog
            analogCh(ii) = (analogEv-1)+ii;
        end
    end
end


figure    
for l = 1:length(lfpCh)
    lfp = bz_GetLFP(lfpCh(l),'noPrompts', true);
    [S1,t,f] = mtspecgramc_fast(single(lfp.data),[2 1],params);
     S = S1;%correctOneOverf(S1,f);

    [pulses] = bz_getAnalogPulses('analogCh',analogCh);
   
    for i =1:3
        if i<=numAnalog
            pulTr = (pulses.stimComb==i);
        else
            pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
        end
        events = pulses.intsPeriods(1,pulTr);
        events = events(((events + 5) <= max(t)) & ((events - 5) > 0));
        thetaProfile_Pre = [];
        thetaProfile_Post = [];
        for pp = 1:length(events)
            events_tmp = abs(t-events(pp));
            [~,idx_temp] = min(events_tmp);
            thetaProfile_Pre(:,:,pp) = S((idx_temp-1),:)';
            thetaProfile_Post(:,:,pp) = S((idx_temp+1),:)';
        end
        subplot(length(lfpCh),3,3*(l-1)+i)
        thetaProfile_Post(:,:,5) = nan;
        plot(f,nanmean(thetaProfile_Pre,3),'k')
        hold on
        plot(f,nanmean(thetaProfile_Post,3),'b')
        ylabel('Power')
        xlabel('Frequency')
        xlim([frange(1) frange(2)])
        if i ==2
           title(strcat('Channel ',num2str(lfpCh(l))))
        end
        clear thetaProfile_Pre thetaProfile_Post
    end
end


