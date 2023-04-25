function sessPeriStimSpectogram

%Plots the average spectogram around a stimulation pulse - during NonREM
%Plots the average spectogram around an UP-DOWN transition

wavelet = 1;
expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cond = 0; % 0 is mEC only, 1 is mEC bilateral
if cond == 0
    pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
    'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};
    analogCh = 2;
elseif cond == 1
    pathToSessionsAll = {'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3'};
    analogCh = [1 2 3];
end

for ii=1:size(pathToSessionsAll,2)
        
    cd(strcat(expPath,pathToSessionsAll{ii}));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);   
    file = dir('*.SlowWaves.events.mat');
    load(file.name);    
    
    if exist([sessionInfo.FileName '.SleepState.states.mat'],'file') 
       load([sessionInfo.FileName '.SleepState.states.mat']);
    else
       disp('No sleep states associated with this session');
       return
    end
    
    if exist([sessionInfo.FileName '.ripples.events.mat'],'file') 
       load([sessionInfo.FileName '.ripples.events.mat']);
    else
       disp('First calculate .ripple file! Skipping');
       return
    end
    
    %Load ripple LFP
    pyrCh = ripples.detectorinfo.detectionchannel;
    lfp = bz_GetLFP(pyrCh);
    
    % Timestamps of UP DOWN transitions 
    stUD = SlowWaves.ints.DOWN(:,1);
    stDU = SlowWaves.ints.DOWN(:,2);
    
    if length(stUD)>250
        stUD = stUD(1:250);
        stDU = stDU(1:250);
    end
    % Timestamps of mEC pulses
    [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
    
    % Event triggered spectogram
    [wavAvg, lfpAvg] = bz_eventWavelet(lfp, stUD, 'twin',[2 2],'plotWave',false,'plotLFP',false);
    freqs = wavAvg.freqs;
    timest = wavAvg.timestamps;
    wavAvgstUD(:,:,ii) = wavAvg.data;
    lfpAvgstUD(:,:,ii) = lfpAvg.data;
    
    [wavAvg, lfpAvg] = bz_eventWavelet(lfp, stDU, 'twin',[2 2],'plotWave',false,'plotLFP',false); 
    wavAvgstDU(:,:,ii) = wavAvg.data;
    lfpAvgstDU(:,:,ii) = lfpAvg.data;
    
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
            if ii< 16
                if aa<=2
                    pulTr = (pulses.stimComb==aa);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==aa);
                end
            else
                pulTr = (pulses.stimPerID'==1 & (pulses.stimComb==2 | pulses.stimComb==3)); % for the CA3 mice, both 2 and 3 are the same
            end
        end
       
        % Only take pulses that happened during NREM sleep
        events = pulses.intsPeriods(1,pulTr)';         
        eventsUD = events(InIntervals(events,SleepState.ints.NREMstate));
        
        events = pulses.intsPeriods(2,pulTr)';  
        eventsDU  = events(InIntervals(events,SleepState.ints.NREMstate));   
        
        if ~isempty(eventsUD)
            [wavAvg, lfpAvg] = bz_eventWavelet(lfp, eventsUD, 'twin',[2 2],'plotWave',false,'plotLFP',false); 
            wavAvgeventUD(:,:,ii) = wavAvg.data;
            lfpAvgeventUD(:,:,ii) = lfpAvg.data;

            [wavAvg, lfpAvg] = bz_eventWavelet(lfp, eventsDU, 'twin',[2 2],'plotWave',false,'plotLFP',false); 
            wavAvgeventDU(:,:,ii) = wavAvg.data;
            lfpAvgeventDU(:,:,ii) = lfpAvg.data;          
        end
    end
end

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
subplot(2,2,1)
contourf(timest*1000,freqs,nanmean(wavAvgstUD,3)',30,'LineColor','none');hold on;
set(gca,'YScale','log');
ylim([freqs(1) freqs(end)]);
colormap jet;
xlabel('time (ms)'); ylabel('frequency (Hz)');
fmin =  freqs(round(length(freqs)/4)); 
fmax =  fmin*2;
lfp_avg = nanmean(lfpAvgstUD,3);
lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;
caxis([30 150])
xlim([-1000 1000])

subplot(2,2,2)
contourf(timest*1000,freqs,nanmean(wavAvgstDU,3)',30,'LineColor','none');hold on;
set(gca,'YScale','log');
ylim([freqs(1) freqs(end)]);
colormap jet;
xlabel('time (ms)'); ylabel('frequency (Hz)');
fmin =  freqs(round(length(freqs)/4)); 
fmax =  fmin*2;
lfp_avg = nanmean(lfpAvgstDU,3);
lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;
caxis([30 150])
xlim([-1000 1000])

subplot(2,2,3)
contourf(timest*1000,freqs,nanmean(wavAvgeventUD,3)',30,'LineColor','none');hold on;
set(gca,'YScale','log');
ylim([freqs(1) freqs(end)]);
colormap jet;
xlabel('time (ms)'); ylabel('frequency (Hz)');
fmin =  freqs(round(length(freqs)/4)); 
fmax =  fmin*2;
lfp_avg = nanmean(lfpAvgeventUD,3);
lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;
caxis([30 150])
xlim([-1000 1000])

subplot(2,2,4)
contourf(timest*1000,freqs,nanmean(wavAvgeventDU,3)',30,'LineColor','none');hold on;
set(gca,'YScale','log');
ylim([freqs(1) freqs(end)]);
colormap jet;
xlabel('time (ms)'); ylabel('frequency (Hz)');
fmin =  freqs(round(length(freqs)/4)); 
fmax =  fmin*2;
lfp_avg = nanmean(lfpAvgeventDU,3);
lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;
caxis([30 150])
xlim([-1000 1000])

saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleSpectogram.png'));
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleSpectogram.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleSpectogram..fig'));
end
