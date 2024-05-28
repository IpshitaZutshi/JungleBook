function plotTrialSpectogram

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[30 275 1800 550]);

%% Plot an example trial     
%sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220919_sess14';
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sessloc)
channelNum = 7;
lfp = bz_GetLFP(channelNum,'noPrompts',true);
file = dir('*.TrialBehavior.Behavior.mat');
load(file.name);
% 
% plotSpectrum(27,lfp,behavTrials,6,1,1,fig2)
% plotSpectrum(29,lfp,behavTrials,6,1,3,fig2)
% plotSpectrum(33,lfp,behavTrials,6,1,5,fig2)

plotSpectrum(55,lfp,behavTrials,6,2,1,1,fig2)
plotSpectrum(60,lfp,behavTrials,6,2,1,2,fig2)
plotSpectrum(34,lfp,behavTrials,6,2,4,1,fig2)
plotSpectrum(42,lfp,behavTrials,6,2,4,2,fig2)


% 
% for ii = 1:6
%     idx = behavTrials.linTrial==0 & behavTrials.stim==0 & behavTrials.correct==1 &behavTrials.lickLoc==(ii-1);
%     lickTS = behavTrials.timestamps(idx,2);
%     twin = 3;
%     nfreqs = 100;
%     ncyc = 7;
%     fBand = [2 150];
% 
%     lfpSpect = [];
%     for iRip = 1:length(lickTS)
%         [~,idx] = min(abs(lfp.timestamps-lickTS(iRip)));
%         intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(1*1250))];
%         wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
%         lfpSpect(:,:,iRip) = log10(abs(wavelet.data));   
%         freqs = wavelet.freqs;                
%     end
%     ts = wavelet.timestamps-wavelet.timestamps(1);
% 
%     subplot(6,1,ii)
%     imagesc(ts,freqs,nanmedian(lfpSpect,3)')
%     set(gca,'YDir','normal')
%     set(gca,'YScale','log')
%     ylim([2 150])
%     caxis([1 4])
% end

end

function plotSpectrum(trialNum,lfp,behavTrials,numrows,numcol,rowloc,colloc,fighandle)

[~,idxStart] = min(abs(lfp.timestamps-behavTrials.timestamps(trialNum,1)));
[~,idxEnd] = min(abs(lfp.timestamps-behavTrials.timestamps(trialNum,2)));

subplot(numrows,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle)
plot(lfp.timestamps(idxEnd-(2*1250):idxEnd+625),lfp.data(idxEnd-(2*1250):idxEnd+625))
hold on
line([behavTrials.timestamps(trialNum,2) behavTrials.timestamps(trialNum,2)],[-2000 3000],'LineWidth',1.5);
xlim([behavTrials.timestamps(trialNum,2)-2 behavTrials.timestamps(trialNum,2)+0.5])
title(strcat('Port :',num2str(behavTrials.toneGain(trialNum)+1),'    TrialNum:',num2str(trialNum)))

subplot(numrows,numcol,numcol*(rowloc-1)+colloc+numcol,'Parent',fighandle)
nfreqs = 100;
ncyc = 15;
fBand = [2 150];
wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',[behavTrials.timestamps(trialNum,2)-2 behavTrials.timestamps(trialNum,2)+0.5]);
lfpSpect = log10(abs(wavelet.data));
imagesc(wavelet.timestamps,wavelet.freqs,lfpSpect');
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
caxis([2 4])
hold on
line([behavTrials.timestamps(trialNum,2) behavTrials.timestamps(trialNum,2)],[2 150],'Color','w','LineWidth',1.5);
xlim([behavTrials.timestamps(trialNum,2)-2 behavTrials.timestamps(trialNum,2)+0.5])

%% Ratio of the power in the 18-35 Hz band and theta band
gamIdx = wavelet.freqs>18 & wavelet.freqs <30;
thetaIdx = wavelet.freqs>6 & wavelet.freqs <12;

powGam = nanmean(lfpSpect(:,gamIdx),2);
powTheta = nanmean(lfpSpect(:,thetaIdx),2);

subplot(numrows,numcol,numcol*(rowloc-1)+colloc+(numcol*2),'Parent',fighandle)
hold on
plot(wavelet.timestamps,powGam./powTheta)
xlim([behavTrials.timestamps(trialNum,2)-2 behavTrials.timestamps(trialNum,2)+0.5])
end