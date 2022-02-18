lfpCh = [8 26 62 6];%[21 17 16 19]; %[18 31 39 34]for IZ4, [59 63 53 37] for IZ5

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
frange = [2 300];


if ~exist('analogEv','var')
    analogEv = 64;
    numAnalog = 2;

    if ~isempty(analogEv)
        for ii = 1:numAnalog
            analogCh(ii) = (analogEv-1)+ii;
        end
    end
end

[pulses] = bz_getAnalogPulses('analogCh',analogCh);
 
figure    
for l = 1:length(lfpCh)
    lfp = bz_GetLFP(lfpCh(l),'noPrompts', true);
    for i =1:3
        if i<=numAnalog
            pulTr = (pulses.stimComb==i);
        else
            pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
        end
        events = pulses.intsPeriods(:,pulTr)';
        events(1,:) =  events(1,:)+1;
        events(2,:) =  events(2,:)-0.5;
        events_pre = events-6;
        
        [specslope_stim,spec_stim] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[4 150],'spectype','fft','ints',events);
        [specslope,spec] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[4 150],'spectype','fft','ints',events_pre);

        subplot(length(lfpCh),3,3*(l-1)+i)
        plot(specslope.freqs,mean(specslope.resid,1),'k')
        hold on
        plot(specslope_stim.freqs,mean(specslope_stim.resid,1),'b')
        ylabel('Corrected power')
        xlabel('Frequency')
        if i ==2
           title(strcat('Channel ',num2str(lfpCh(l))))
        end
    end
end


