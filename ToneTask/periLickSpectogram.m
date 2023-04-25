channelNum = 79;
lfp = bz_GetLFP(channelNum,'noPrompts',true);

file = dir('*.TrialBehavior.Behavior.mat');
load(file.name);

idx = behavTrials.linTrial==0 & behavTrials.stim==0 ;
lickTS = behavTrials.timestamps(idx,2);
twin = 3;
nfreqs = 100;
ncyc = 7;
fBand = [2 150];

lfpSpect = [];
for iRip = 1:length(lickTS)
    [~,idx] = min(abs(lfp.timestamps-lickTS(iRip)));
    intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(1*1250))];
    wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
    lfpSpect(:,:,iRip) = log10(abs(wavelet.data));   
    freqs = wavelet.freqs;                
end
ts = wavelet.timestamps-wavelet.timestamps(1);

figure
set(gcf,'Renderer','painters')
set(gcf,'color','w')
subplot(1,2,1)
imagesc(ts,freqs,nanmedian(lfpSpect,3)')
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
caxis([1 4])


idx = behavTrials.linTrial==0 & behavTrials.stim==1 ;
lickTS = behavTrials.timestamps(idx,2);
twin = 3;
nfreqs = 100;
ncyc = 7;
fBand = [2 150];

lfpSpect = [];
for iRip = 1:length(lickTS)
    [~,idx] = min(abs(lfp.timestamps-lickTS(iRip)));
    intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(0.5*1250))];
    wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
    lfpSpect(:,:,iRip) = log10(abs(wavelet.data)); 
    freqs = wavelet.freqs;                
end

subplot(1,2,2)
imagesc(1:1:size(lfpSpect,1),freqs,nanmedian(lfpSpect,3)')
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
caxis([1 4])