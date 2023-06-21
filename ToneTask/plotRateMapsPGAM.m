function plotRateMapsPGAM(resultsPGAM)

relDistCells = resultsPGAM.tuned_relDistStop_mi;
lickCells = resultsPGAM.kernelStrength_licks;

cellID  = union(relDistCells,lickCells);

file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file(1).name);

file = dir(['*TrialBehavior.Behavior.mat']);
load(file(1).name);

for zz = 1:3
    
    if zz ==1 
        idCells = relDistCells;
    elseif zz == 2
        idCells = lickCells;
    else
        idCells = cellID;
    end
    spaceMap{zz} = [];
    toneMap{zz} = [];
    
    for ii = 1:length(idCells)

        % Average spatial and tone ratemaps
        dataMat = [];
        dataMatTone = [];
        for jj = 2:7           
            dataMat = [dataMat;firingMaps.forward.rateMaps{idCells(ii)}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{idCells(ii)}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end 

        spaceMap{zz}(ii,:) = nanmean(dataMat,1);        
        toneMap{zz}(ii,:) = nanmean(dataMatTone,1);     
    end

end

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);

figure
for ii = 1:3
    
    [maxTone,idxTone] = max(toneMap{ii},[],2);
    sortidx =1:length(idxTone);
    %[~,sortidx] = 1:length(idxTone);%sort(idxTone,'ascend');
    normtoneMap = (toneMap{ii}-nanmean(toneMap{ii},2))./nanstd(toneMap{ii},[],2);
    normspaceMap = (spaceMap{ii}-nanmean(spaceMap{ii},2))./nanstd(spaceMap{ii},[],2);
        
    colormap(YlGnBu)
    subplot(2,3,ii)
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    caxis([-1 4])
    subplot(2,3,ii+3)
    imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
    caxis([-1 4])
end
end

