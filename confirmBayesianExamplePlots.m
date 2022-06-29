function confirmBayesianExamplePlots(bayes_allSpikes,cell_metrics,ripples, firingMaps)


sess = 0; %0 for baseline, 1 for stim
stem = round((100/175)*110);
lfp = bz_GetLFP(ripples.detectorinfo.detectionchannel);

rdBlu=cbrewer('div', 'Spectral', (length(cell_metrics.UID)));
rdBlu(rdBlu>1) = 1;

sigRips = 1:1:20;%bayes_allSpikes.all{1}.significant_Imax;

rateMaps = [];
for ii = 1:length(cell_metrics.UID)
    rateMaps(ii,:) = firingMaps.rateMaps{ii}{sess+1}(stem:95);
end
[~,b] = max(rateMaps,[],2);
[~,sortidx] = sort(b);
    
for ii = 1:20
    
    figure
    subplot(3,2,[1 3])
    imagesc(zscore(rateMaps(sortidx,:),[],2))
    xlabel('Position on center arm')
    ylabel('Neuron ID')
    caxis([-2 3])
    
    TS = [bayes_allSpikes.all{1}.ripple_HSE.timestamps(sigRips(ii),1)-0.05 ...
        bayes_allSpikes.all{1}.ripple_HSE.timestamps(sigRips(ii),2)+0.05];
    TStrain = round(TS(1)*1250):1:round(TS(2)*1250);
    
    
    %% Buildspikeraster
    Spikerast = [];
    Spikerast(length(cell_metrics.UID),length(TStrain)) = 0;
    for jj = 1:length(cell_metrics.UID)
        tSp = cell_metrics.spikes.times{jj}(cell_metrics.spikes.times{jj}>=TS(1) & cell_metrics.spikes.times{jj}<=TS(2));
        tSp = round(tSp*1250);
        [~,idxRast] = ismember(tSp,TStrain);
        Spikerast(jj,idxRast) = 1;
    end
    
    subplot(3,2,2)
    [~,startIdx] = min(abs(lfp.timestamps-TS(1)));
    [~,endIdx] = min(abs(lfp.timestamps-TS(2)));
    plot(lfp.timestamps(startIdx:endIdx),lfp.data(startIdx:endIdx))
    xlim([TS(1) TS(2)])
    
    subplot(3,2,4)
    sortloc = 1;
    for jj=sortidx'  
        if sum(Spikerast(jj,:))==0
           sortloc = sortloc +1;
        else
            plot_raster(find(Spikerast(jj,:)),sortloc,1,rdBlu(sortloc,:),2)
            sortloc = sortloc +1;
            hold on
        end
    end
    set(gca,'YDir','rev')
    xlim([0 size(Spikerast,2)])
    
    subplot(3,2,5)
    for kk= 1:size(bayes_allSpikes.all{1}.ripple_score_super,2)
        if (bayes_allSpikes.all{1}.ripple_score_super(kk).eventID == sigRips(ii))
            imagesc(bayes_allSpikes.all{1}.ripple_score_super(kk).prob)
        end
    end
    
    subplot(3,2,6)
    imagesc(zscore(bayes_allSpikes.all{1}.template_beh(sortidx,:),[],2))
end
end