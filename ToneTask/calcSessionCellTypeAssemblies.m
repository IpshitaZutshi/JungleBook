function [AllPlaceWt, AllToneWt, tBins] = calcSessionCellTypeAssemblies

usePGAM = 1;
plotfig = 0;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};     
%sess = {'IZ48\Final\IZ48_230703_sess21'};

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

AllPlaceWt(length(sess),22) = nan;
AllToneWt(length(sess),22) = nan;

for ss = 1:length(sess)
    try
        %% Load files
        cd(strcat(expPath,sess{ss})) 
    
        file = dir(['*spikes.cellinfo.mat']);
        load(file(1).name);
        file = dir(['*cell_metrics.cellinfo.mat']);
        load(file.name);
        file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
        load(file(1).name);
        file = dir(['*TrialBehavior.Behavior.mat']);
        load(file(1).name);
        file = dir(['*Tracking.Behavior.mat']);
        load(file(1).name);
        [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
        
        dtime = mean(diff(tracking.timestamps));
        spkData = bz_SpktToSpkmat(spikes.times,'dt',dtime, 'win',[behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
        spikeMat = spkData.data';
        
        %Only select pyramidal cells
        logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
        keepCells = logicalVector;
        spikeMat = spikeMat(keepCells,:);
        cellId = find(keepCells);
        
        SpikeCount = zscore(spikeMat,[],1);
        AssemblyTemplates = assembly_patterns(SpikeCount);
        
        Vectors = [];
        assembliesID = [];
        assembliesIDNeg = [];
        for aa = 1:size(AssemblyTemplates,2)
            assembliesVector = AssemblyTemplates(:,aa);
            % flip if negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            Vectors(:,aa) = assembliesVector;
        
            %For each assembly, define significant cells
            sigCells = find(assembliesVector> 2*std(assembliesVector));
            assembliesID{aa} = cellId(sigCells);

            sigCells = find(assembliesVector< -2*std(assembliesVector));
            assembliesIDNeg{aa} = cellId(sigCells);

        end  
        
        Activities = assembly_activity(Vectors,SpikeCount);
        
        %% Also classify cells as tone or place
        
        toneCellLog = zeros(spikes.numcells,1);
        placeCellLog = zeros(spikes.numcells,1);
        
        if usePGAM == 0
            for ii = 1:length(keepCells)
                if keepCells(ii) == 1
            
                    %% Check if its a tone cell
                    dataMatTone = [];
                    dataMat = [];
            
                    for jj = 2:7    
                        dataMat = [dataMat;firingMaps.forward.rateMaps{ii}{jj}]; % Place
                        a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear');
                        dataMatTone = [dataMatTone;a];
                    end  
            
                    spaceMap = nanmean(dataMat,1);
                    toneMap = nanmean(dataMatTone,1);        
              
                    corrTone = []; 
                    corrSpace = [];
                    for pp = 1:6
                       for jj = (pp+1):6      
                           a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
                           corrSpace = [corrSpace a(1,2)];                  
                           a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
                           corrTone = [corrTone a(1,2)];
                       end
                    end     
                    toneCorr = nanmean(corrTone);  
                    spaceCorr = nanmean(corrSpace);
            
                    %% Detect fields
                    Field_Info_tone = detectFields(toneMap,'maxFieldSize',40);
                    Field_Info_space = detectFields(spaceMap);
                
                    if ~isempty(Field_Info_tone) && (toneCorr > 0.1)
                        toneCellLog(ii) = 1;       
                    end
                
                    if ~isempty(Field_Info_space) && (spaceCorr > 0.1)
                        placeCellLog(ii) = 1;
                    end   

                    % if placeCellLog(ii)==1 && toneCellLog(ii) ==1
                    %     if toneCorr > spaceCorr
                    %         placeCellLog(ii) = 0;
                    %     else
                    %         toneCellLog(ii) = 0;
                    %     end
                    % end
                end
            end
        else
            PGAMpath = 'C:\Data\PGAMAnalysis\processedData\';
            fileloc = strcat(PGAMpath,sessionInfo.FileName,'\results_struct.mat');
            load(fileloc)
            % Get neuronID 
            neuronID = [];
            for aa = 1:size(results,2)
                neuronID(aa) = results(aa).neuron+1;
            end
        
            numCells = length(cell_metrics.UID);
        
            for aa = 1:numCells
                sigMat(1:2) = 0;
                mutInfo(1:2) = nan;
                if ~strcmp(cell_metrics.putativeCellType(aa),'Pyramidal Cell')
                    continue
                end        
                idx = find(neuronID == aa);
                if isempty(idx)
                    continue
                end
                pyrID(aa) = 1;
                for id = [1 4]
                    if (results(idx(id)).pval < (10^-5) && ~isnan(results(idx(id)).mutual_info))
                        if id == 1 % Place
                            sigMat(1) = 1;
                            mutInfo(1) = results(idx(id)).mutual_info;
                        else % Frequency
                            sigMat(2) = 1;
                            mutInfo(2) = results(idx(id)).mutual_info;
                        end
                    end
                end
                if sigMat(1)==1 && sigMat(2)==0
                    placeCellLog(aa) = 1;        
                elseif sigMat(1)==0 && sigMat(2)==1
                    toneCellLog(aa) = 1;
                elseif sigMat(1)==1 && sigMat(2)==1 % Both are significant
                    if mutInfo(1)> mutInfo(2)
                        placeCellLog(aa) = 1;
                    elseif mutInfo(2)> mutInfo(1)
                        toneCellLog(aa) = 1;
                    end
                end
            end   
        end
        
        %% For each assembly, calculate the proportion of cells that are tone vs place
        toneCell = toneCellLog(keepCells);
        placeCell = placeCellLog(keepCells);
        fractTonePlace = [];
    
        for aa = 1:length(assembliesID)              
            fractTonePlace(aa,1) = sum(Vectors(logical(toneCell),aa));
            fractTonePlace(aa,2) = sum(Vectors(logical(placeCell),aa));   
        end
        
        %% Sort assemblies around lick time and plot the place versus tone distribution by sorting assembly by latency
        assemblyTimes = [];
        psth1 = [];
        for aa = 1:size(Activities,1)
            idx = find(Activities(aa,:) > 1.5*std(Activities(aa,:)));
            assemblyTimes{aa} = spkData.timestamps(idx);
            [stccg, t] = CCG({assemblyTimes{aa} behavTrials.timestamps(behavTrials.linTrial==0,2)},[],'binSize',0.1,'duration',6,'norm','rate');
            psth1(aa,:) = stccg(:,2,1);        
        end
        
        %find time of lick =
        idxLickEnd = find(t>0.1,1,'first');
        idxLickStart = find(t>-2,1,'first');
        [~,maxPSTH]  = max(psth1,[],2);
        
        %Only keep assemblies that have a peak between -1.5 and 0.2 s around the
        %lick
        toKeep = find(maxPSTH>=idxLickStart & maxPSTH<=idxLickEnd);
        tBins = t(idxLickStart:idxLickEnd);
        
        placeWt = [];
        toneWt = [];
        placeWt(1:length(tBins)) = nan;
        toneWt(1:length(tBins)) = nan;
        
        % Now for each assembly,find the peak location, and the weighted sum of
        % place/ tone cells.
        for tt = 1:length(toKeep)
            maxLoc = maxPSTH(toKeep(tt));
            tLoc = t(maxLoc);
            idxLoc = tBins==tLoc;
        
            %Weight of assemblies
            pWt = sum(Vectors(logical(placeCell),toKeep(tt)));  
            tWt = sum(Vectors(logical(toneCell),toKeep(tt)));  
            placeWt(idxLoc) = nanmean([pWt placeWt(idxLoc)]);
            toneWt(idxLoc) = nanmean([tWt toneWt(idxLoc)]);
        end
        
        AllPlaceWt(ss,:) = placeWt;
        AllToneWt(ss,:) = toneWt;
    
        if plotfig
            [~,sortidx] = sort(maxPSTH,'descend');
            sortedtoKeep = sortidx;
        
            figure
            subplot(1,3,[1 2])
            imagesc(zscore(psth1(sortedtoKeep,:),[],2))
            set(gca,'YDir','normal')
            
            subplot(1,3,3)
            plot(fractTonePlace(sortedtoKeep,1),1:length(sortedtoKeep),'Color','m','LineWidth',1.5)
            hold on
            plot(fractTonePlace(sortedtoKeep,2),1:length(sortedtoKeep),'Color','k','LineWidth',1.5)
            ylim([1 length(sortedtoKeep)])
            
          %  save('sessAssemblyPlot.mat','Vectors','psth1','sortedtoKeep','toneCellLog','placeCellLog','keepCells','assembliesID','Activities')
        end
    catch
    end
end

%save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Assemblycontribution.mat','AllPlaceWt','AllToneWt','tBins')
end