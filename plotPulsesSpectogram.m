function plotPulsesSpectogram

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

if exist([sessionInfo.FileName '.region.mat'],'file') 
    load([sessionInfo.FileName '.region.mat']);
    pyrCh = 9;%region.CA1sp;
else
    disp('First calculate .region file to identify ripple channel! Skipping');
    return
end

lfp = bz_GetLFP(pyrCh,'noPrompts', true);
[pulses] = bz_getAnalogPulsesSine('analogCh',[64 65]);

params.Fs = lfp.samplingRate; params.fpass = [2 100]; params.tapers = [3 5]; params.pad = 1;

[S,t,f] = mtspecgramc(single(lfp.data),[2 1],params);
S1 = log10(S); % in Db
S_det= bsxfun(@minus,S1,polyval(polyfit(f,mean(S1,1),2),f)); % detrending

nf = 1;

figure(nf)
k= 0;
for i = 1:size(pulses.intsPeriodsSess,2)
    k= k+1;
    %subplot(3,4,k)
    subplot(3,size(pulses.intsPeriodsSess,2),i)
    [~,idxStart] = min(abs(t-pulses.intsPeriodsSess(1,i)));
    [~,idxEnd] = min(abs(t-pulses.intsPeriodsSess(2,i)));
    imagesc(t((idxStart-3):(idxEnd+3)),f,S_det((idxStart-3):(idxEnd+3),:)',[-1.5 1.5]);
    set(gca,'YDir','normal'); ylabel('Freqs');

    %subplot(3,4,k+4)
    subplot(3,size(pulses.intsPeriodsSess,2),i+size(pulses.intsPeriodsSess,2))
    idxstart = find(pulses.intsPeriods(1,:)==pulses.intsPeriodsSess(1,i));
    idxend = find(pulses.intsPeriods(2,:)==pulses.intsPeriodsSess(2,i));
    
    startTS = pulses.intsPeriods(1,idxstart:idxend);
    endTS = pulses.intsPeriods(2,idxstart:idxend);
    hold on
    for ll = 1:length(startTS)
        amp = pulses.amplitude(pulses.timestamps(:,1)==startTS(ll));
        [~,idxTemp] = min(abs(t-startTS(ll)));
        profilePre(:,:,ll) = S((idxTemp-1),:)';
        profilePost(:,:,ll) = S((idxTemp+1),:)';
        plot([startTS(ll) startTS(ll)], [0 amp],'b')
        plot([endTS(ll) endTS(ll)], [0 amp],'k')
    end
    xlim([t(idxStart-3) t(idxEnd+3)]);
    
    %subplot(3,4,k+8)
    subplot(3,size(pulses.intsPeriodsSess,2),i+2*(size(pulses.intsPeriodsSess,2)))    
    plot(f,nanmean(profilePre,3),'k')
    hold on
    plot(f,nanmean(profilePost,3),'b')
    ylabel('Power')
    xlabel('Frequency') 
    
    clear profilePre profilePost
    
end

saveas(figure(nf),strcat('summ\Spectogram.png'));
saveas(figure(nf),strcat('summ\Spectogram.eps'));
saveas(figure(nf),strcat('summ\Spectogram.fig'));
nf = nf+1;

figure(nf)
events = pulses.intsPeriods(1,:);
events = events(((events + 5) <= max(t)) & ((events - 5) > 0));
thetaProfile_Pre = [];
thetaProfile_Post = [];

for pp = 1:length(events)
    events_tmp = abs(t-events(pp));
    [~,idx_temp] = min(events_tmp);
    thetaProfile_Pre(:,:,pp) = S((idx_temp-1),:)';
    thetaProfile_Post(:,:,pp) = S((idx_temp+1),:)';
end

freqRange = unique(pulses.freq);
for ff = 1:length(freqRange)
    
    subplot(1,length(freqRange),ff) 
    freqidx = pulses.freq == freqRange(ff);
    plot(f,nanmean(thetaProfile_Pre(:,:,freqidx),3),'k')
    hold on
    plot(f,nanmean(thetaProfile_Post(:,:,freqidx),3),'b')
    ylabel('Power')
    xlabel('Frequency')
    title(['Stim Frequency = ' num2str(freqRange(ff)) ' Hz'])
    
%     subplot(2,length(freqRange),ff+length(freqRange)) 
%     
%     events = pulses.intsPeriods(:,freqidx)';
%     events(1,:) =  events(1,:)+1;
%     events(2,:) =  events(2,:)-0.5;
%     events_pre = events-6;
% 
%     [specslope_stim,~] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[2 100],'spectype','fft','ints',events);
%     [specslope,~] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[2 100],'spectype','fft','ints',events_pre);
% 
%     plot(specslope.freqs,mean(specslope.resid,1),'k')
%     hold on
%     plot(specslope_stim.freqs,mean(specslope_stim.resid,1),'b')
%     ylabel('Corrected power')
%     xlabel('Frequency')

end

% saveas(figure,strcat('summ\PowerVsFreq.png'));
% saveas(figure,strcat('summ\PowerVsFreq..eps'));
% saveas(figure,strcat('summ\PowerVsFreq..fig'));

end