function periLickSpectogram

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[30 275 1800 550]);

channelNum = 7;
lfp = bz_GetLFP(channelNum,'noPrompts',true);

file = dir('*.TrialBehavior.Behavior.mat');
load(file.name);

for ii = 1:6
    for jj = 1:2
        
        idx = behavTrials.linTrial==0 & behavTrials.stim==0 & behavTrials.correct==1 &behavTrials.lickLoc==(ii-1);
        if jj == 1
            lickTS = behavTrials.timestamps(idx,2);
        else
            lickTS = behavTrials.timestamps(idx,1);
        end
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
            
        subplot(6,2,2*(ii-1)+jj)
        imagesc(ts,freqs,nanmedian(lfpSpect,3)')
        set(gca,'YDir','normal')
        set(gca,'YScale','log')
        ylim([2 150])
        caxis([2 4.5])
        colormap jet
        hold on
        line([3 3],[2 150],'Color','w','LineWidth',1.5);
        if jj == 1 && ii==1
            title('Forward')
        elseif jj == 2 && ii == 1
            title('Return')
        end
    end
end

end