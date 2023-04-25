function compileSineWavePSTH

pathtoSess = {'E:\IZ12\Home cage\IZ12_792um_200213_sess15','E:\IZ17\Home Cage\IZ17_488um_200810_sess12',...
    'E:\IZ18\Home cage\IZ18_360um_200810_sess13','F:\IZ21\Homecage\IZ21_0um_200827_sess11','G:\IZ23\Home cage\IZ23_360um_201006_sess13',...
  'F:\IZ30\IZ30_792um_201123_sess8','F:\IZ31\IZ31_0um_201123_sess8'};%'F:\IZ24\IZ24_288um_200929_sess9''F:\IZ25\IZ25_0um_201008_sess9',  'F:\IZ25\IZ25_0um_201009_sess10',...
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
freqs = [0 10 20 50];
for ff = 1:length(freqs)
    psth{ff} = [];    
end
cellClass = [];

for ss = 1:length(pathtoSess)
    
    cd(pathtoSess{ss})
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    
    % If a pulses file does not exist, make one
    if exist([sessionInfo.FileName '.pulses.events.mat'])
        load([sessionInfo.FileName '.pulses.events.mat']);
    else 
        [pulses] = bz_getAnalogPulsesSine;
    end
    
    if ~isempty(dir('*Kilosort*')) % if kilosort has been run
        spikes = loadSpikes;
    else
        disp(strcat('First run Kilosort for ',pathtoSess{ss}))
    end
    
    [cell_metrics] = bz_CellMetricsSimpleIZ(spikes,sessionInfo.rates.wideband);
    cellClass = [cellClass; cell_metrics.putativeClass];
    
    %% Extract the different types of stim
    
    for ff =  1:length(freqs)
        st = pulses.intsPeriods(1,pulses.freq == freqs(ff))';
        
        for unit = 1:length(spikes.ids)
            
            [stccg, t] = CCG({spikes.times{unit} st},[],'binSize',0.2,'duration',21,'normtype','rate');
            psth{ff} = [psth{ff}; stccg(:,2,1)'];
        end
    end       

end

figure
set(gcf,'Renderer','painters')
col = [85/243 85/243 85/243;...
    224/243 163/243 46/243;...
    175/243 54/243 60/243;...
    56/243 61/243 150/243];
    
for cc = 1:2
    idx = cellClass==cc;
    for ff = 1:length(freqs)
        
        temppsth = psth{ff}(idx,:);
        subplot(2,5,5*(cc-1)+ff)
        imagesc(1,1:1:size(temppsth,1),zscore(temppsth,[],2))
        caxis([-2 2])

        subplot(2,5,5*(cc-1)+5)
        plot(t,nanmean(zscore(temppsth,[],2)),'color',col(ff,:))
        hold on
        xlim([t(1) t(end)]);
    end
end

saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\sineWavePSTH.png'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\sineWavePSTH.fig'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\sineWavePSTH.eps'),'epsc');
end
