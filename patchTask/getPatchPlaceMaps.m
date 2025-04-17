function getPatchPlaceMaps(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveLoc',[],@isstr);
addParameter(p,'forceCalculate',false,@isstr);
addParameter(p,'plotfig',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveLoc = p.Results.saveLoc;
forceCalculate = p.Results.forceCalculate;
plotfig = p.Results.plotfig;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

if isempty(saveLoc)
    saveLoc = strcat(basepath,'\Maps');
    if ~isfolder('Maps')
        mkdir('Maps')
    end    
end

fprintf('Computing place fields\n'); 

portloc = [2 16 37 66 86 110 122.5]/2.5;

%% Assign spike position to each spike
if ~isempty(dir([basepath filesep '*.spikeData.cellinfo.mat']))
    file = dir([basepath filesep '*.spikeData.cellinfo.mat']);
    load(file.name);
else
    for unit = 1:length(spikes.UID)
        spikeData.posIdx{unit} = [];
        spikeData.pos{unit} = [];
        [idx] = InIntervals(spikes.times{unit},[tracking.timestamps(1) tracking.timestamps(end)]); 
        tsBehav = spikes.times{unit}(idx);
        if isempty(tsBehav)
            spikeData.posIdx{unit} = [];
        else
            for tt = 1:length(tsBehav)
                [~,closestIndex] = min(abs(tracking.timestamps-tsBehav(tt)));
                spikeData.posIdx{unit}(tt) = closestIndex;
            end
        end
        spikeData.pos{unit} = tracking.position.y(spikeData.posIdx{unit});
    end
    save([sessionInfo.FileName '.spikeData.cellinfo.mat'],'spikeData'); 
end

%% Group spiking by trials
if ~isempty(dir([basepath filesep '*.rateMapsTrial.cellinfo.mat'])) && ~forceCalculate
    file = dir([basepath filesep '*.rateMapsTrial.cellinfo.mat']);
    load(file.name);
else
    for pf = 1:(size(behavTrials.timestamps,1)-1)    
        [idx] = InIntervals(tracking.timestamps,[behavTrials.timestamps(pf) behavTrials.timestamps(pf+1)]);
        portRun(pf,:) = [behavTrials.port(pf) behavTrials.port(pf+1)];
        rewarded(pf,:) = [behavTrials.reward_outcome(pf) behavTrials.reward_outcome(pf+1)];

        % current high patch
        highProbPatch = mean(behavTrials.ports_probability(pf,1:3)) > mean(behavTrials.ports_probability(pf,5:7));    

        % Determine whether the next lick was in a high-probability patch
        highPatch(pf) = ismember(behavTrials.port(pf+1), 1:3) & highProbPatch | ...
                                ismember(behavTrials.port(pf+1), 5:7) & ~highProbPatch;

        %Approximate running direction    
        direction(pf) = portRun(pf,2)>portRun(pf,1);

        if sum(idx)>0
            positions{pf} = [tracking.timestamps(idx) tracking.position.x(idx) tracking.position.y(idx)]; 
        else
            positions{pf} = [];
        end
    end

    firingMaps = bz_getRateMaps(positions,spikes,'xRange',[0 6],'yRange',[0 125], 'binSize',2.5,'saveMat',false);
    firingMaps.portRun = portRun;
    firingMaps.rewarded = rewarded;
    firingMaps.highPatch = highPatch;
    firingMaps.direction = direction;
    
    save([sessionInfo.FileName '.rateMapsTrial.cellinfo.mat'],'firingMaps'); 
end

col1 = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    0.5 0.5 0.5;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

if plotfig
    
    plotSpikeData = 0;
    plotTrialMaps = 1;
    
    %% First plot the spike data plots
    if plotSpikeData
        for cellNum = 29%1:length(spikeData.pos)  
            figure
            set(gcf,'Color','w')
            set(gcf,'Position',[182 166 1585 762])
            plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1.2)
            hold on
            scatter(tracking.timestamps(spikeData.posIdx{cellNum}),tracking.position.y(spikeData.posIdx{cellNum}),8,'r','filled')

            for kk = 1:7
                idx = find(behavTrials.port==kk);
                for ii= 1:length(idx)
                    lickTS = behavTrials.timestamps(idx(ii));
                    [~,ixTS] = min(abs(tracking.timestamps-lickTS));    
                    %scatter(tracking.timestamps(ixTS),tracking.position.y(ixTS),10,col1(kk,:),'filled')                        
                    scatter(tracking.timestamps(ixTS),tracking.position.y(ixTS),10,'k','filled')                        
                end
            end

            hold on
            plot(behavTrials.timestamps,behavTrials.patch_number*120,'b')
            ylim([0 125])
            xlabel('Time(s)')
            ylabel('Position on track (cm)')
            box off 

            % saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.png'],'png');
            % saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.fig'],'fig');
            % saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.eps'],'epsc');
            % close all;
        end
    end
    
    if plotTrialMaps
        % Figure out unique port combinations
        data = unique(firingMaps.portRun, 'rows');

        % Sort 
        ascending = data(data(:,1) < data(:,2), :);
        descending = data(data(:,1) > data(:,2), :);
        
        % Sort each section
        ascending = sortrows(ascending);
        descending = sortrows(descending, [1, 2], 'descend');
        
        % Concatenate ascending first, then descending
        sorted_data = [ascending; descending];

        % Stack the ratemaps according to this sorting
        idx = cell(size(sorted_data, 1), 1);
        idxSorted = cell(size(sorted_data, 1), 1);
        idxHighPatch = cell(size(sorted_data, 1), 1);
        idxLowPatch = cell(size(sorted_data, 1), 1);

        highPatchLick = cell(size(sorted_data, 1), 1);
        for ss = 1:size(sorted_data, 1)
            idx{ss} = find(ismember(firingMaps.portRun, sorted_data(ss, :), 'rows'));
            highPatchLick{ss} = firingMaps.highPatch(idx{ss});
            % Sort highPatchLick descending, so that trials are separated
            % by high versus low patch
            [~,idxPatch] = sort(highPatchLick{ss});
            idxSorted{ss} = idx{ss}(idxPatch);
            idxHighPatch{ss} = idx{ss}(highPatchLick{ss}==0);
            idxLowPatch{ss} = idx{ss}(highPatchLick{ss}==1);
        end

        numCells = length(firingMaps.rateMaps);

        for cellNum = 157
            sortedRateMaps{cellNum} = [];
            sortedAvgRateMapsHighPatch{cellNum} = [];
            sortedAvgRateMapsLowPatch{cellNum} = [];
            for ss = 1:size(sorted_data, 1)
                if ~isempty(idxSorted{ss})
                    for trialNum = 1:length(idxSorted{ss})
                        sortedRateMaps{cellNum} = [sortedRateMaps{cellNum}; firingMaps.rateMaps{cellNum}{idxSorted{ss}(trialNum)}];
                    end
                else
                    ratemap = nan(1,50);
                    sortedRateMaps{cellNum} = [sortedRateMaps{cellNum};ratemap]; 
                end

                % Make average maps
                avgRM = [];
                if ~isempty(idxHighPatch{ss})
                    for trialNum = 1:length(idxHighPatch{ss})
                        avgRM = [avgRM; firingMaps.rateMaps{cellNum}{idxHighPatch{ss}(trialNum)}];
                    end

                    % Find the fraction of non-NaN values in each column
                    fracNonNan = sum(~isnan(avgRM), 1) ./ size(avgRM, 1);
                    
                    % Keep only columns where at least 10% of the values are non-NaN
                    avgRM_filtered = avgRM(:, fracNonNan >= 0.20);
                    
                   % Compute the average over the filtered columns
                    if ~isempty(avgRM_filtered)
                        % Take the average of filtered data
                        avg_filtered = nanmean(avgRM_filtered, 1);
                        % Pad with NaNs to match the original size
                        padded_avg = nan(1, size(avgRM, 2));
                        padded_avg(fracNonNan >= 0.20) = avg_filtered;
                    else
                        padded_avg = nan(1, size(avgRM, 2));
                    end
                    sortedAvgRateMapsHighPatch{cellNum} = [sortedAvgRateMapsHighPatch{cellNum}; padded_avg];
                else
                    ratemap = nan(1,50);
                    sortedAvgRateMapsHighPatch{cellNum} = [sortedAvgRateMapsHighPatch{cellNum};ratemap];
                end

                % Make average maps
                avgRM = [];
                if ~isempty(idxLowPatch{ss})
                    for trialNum = 1:length(idxLowPatch{ss})
                        avgRM = [avgRM; firingMaps.rateMaps{cellNum}{idxLowPatch{ss}(trialNum)}];
                    end
                     % Find the fraction of non-NaN values in each column
                    fracNonNan = sum(~isnan(avgRM), 1) ./ size(avgRM, 1);
                    
                    % Keep only columns where at least 10% of the values are non-NaN
                    avgRM_filtered = avgRM(:, fracNonNan >= 0.20);
                    
                   % Compute the average over the filtered columns
                    if ~isempty(avgRM_filtered)
                        % Take the average of filtered data
                        avg_filtered = nanmean(avgRM_filtered, 1);
                        % Pad with NaNs to match the original size
                        padded_avg = nan(1, size(avgRM, 2));
                        padded_avg(fracNonNan >= 0.20) = avg_filtered;
                    else
                        padded_avg = nan(1, size(avgRM, 2));
                    end
                    sortedAvgRateMapsLowPatch{cellNum} = [sortedAvgRateMapsLowPatch{cellNum}; padded_avg];
                else
                    ratemap = nan(1,50);
                    sortedAvgRateMapsLowPatch{cellNum} = [sortedAvgRateMapsLowPatch{cellNum};ratemap];
                end                
            end

            fig2 = figure;            
            set(gcf,'Color','w')
            subplot(2,3,[1 4])
            h = imagesc(sortedRateMaps{cellNum});
            set(h, 'AlphaData', ~isnan(sortedRateMaps{cellNum}))
            %caxis([0 20])

            subplot(2,3,2)
            h = imagesc(sortedAvgRateMapsHighPatch{cellNum});
            set(h, 'AlphaData', ~isnan(sortedAvgRateMapsHighPatch{cellNum}))
            caxis([0 max(max(sortedAvgRateMapsHighPatch{cellNum}))])
            title('High probability','Color','r')
            for ss = 1:size(sorted_data,1)
                line([portloc(sorted_data(ss,1)) portloc(sorted_data(ss,1))],[ss-0.5 ss+0.5],'LineWidth',2,'Color','k')
                hold on
                line([portloc(sorted_data(ss,2)) portloc(sorted_data(ss,2))],[ss-0.5 ss+0.5],'LineWidth',2,'Color','k')
            end
            colorbar

            
            subplot(2,3,3)
            h = imagesc(sortedAvgRateMapsLowPatch{cellNum});
            set(h, 'AlphaData', ~isnan(sortedAvgRateMapsLowPatch{cellNum}))
            caxis([0 max(max(sortedAvgRateMapsLowPatch{cellNum}))])
            title('Low probability','Color','b')
            for ss = 1:size(sorted_data,1)
                line([portloc(sorted_data(ss,1)) portloc(sorted_data(ss,1))],[ss-0.5 ss+0.5],'LineWidth',2,'Color','k')
                hold on
                line([portloc(sorted_data(ss,2)) portloc(sorted_data(ss,2))],[ss-0.5 ss+0.5],'LineWidth',2,'Color','k')
            end
            colorbar

            subplot(2,3,6);
            timeaxis = 1:50;
            if cellNum == 157
                data1 = sortedAvgRateMapsHighPatch{cellNum}([7 8 9 10],:);
                data2 = sortedAvgRateMapsLowPatch{cellNum}([7 8 9 10],:); 
                %data3 = sortedAvgRateMapsHighPatch{cellNum}([14 17 19],:);
            else
                data1 = sortedAvgRateMapsHighPatch{cellNum}([15 18 20 21],:);
                data2 = sortedAvgRateMapsLowPatch{cellNum}([15 18 20 21],:); 
                data3 = sortedAvgRateMapsHighPatch{cellNum}([14 17 19],:);                
            end
            plotAvgStd(data1,2,3,6,fig2,timeaxis','r',0)
            hold on                       
            plotAvgStd(data2,2,3,6,fig2,timeaxis','b',0)            
            %plotAvgStd(data3,2,3,6,fig2,timeaxis','k',0)
            ylabel('Average firing rate')
            xlabel('Position on track')
            title(strcat('Cell num',num2str(cellNum)))


%          saveas(gcf,['FiringMap',filesep ,'cell_' num2str(cellNum) '.png'],'png');
%          saveas(gcf,['FiringMap',filesep ,'cell_' num2str(cellNum) '.fig'],'fig');
%          saveas(gcf,['FiringMap',filesep ,'cell_' num2str(cellNum) '.eps'],'epsc');
%          close all;
        end        

    end  
end
end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col, useMedian)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);
    if ~useMedian
        meanpsth = nanmean(array,1);
        stdpsth = nanstd(array,1)./sqrt(size(array,1));
        lArr  = meanpsth-stdpsth;
        uArr = meanpsth+stdpsth;
    else
        meanpsth = nanmedian(array,1);
        % Bootstrapping to estimate variability of the median
        n_bootstraps = 1000;
        bootstrap_medians = zeros(n_bootstraps, size(array, 2));
        
        for i = 1:n_bootstraps
            resample_indices = randi([1, size(array, 1)], size(array, 1), 1);  % Generate random indices with replacement
            bootstrap_sample = array(resample_indices, :);  % Create bootstrap sample
            bootstrap_medians(i, :) = nanmedian(bootstrap_sample);  % Compute median of the bootstrap sample
        end
        
        % Compute the 2.5th and 97.5th percentiles for the bounds
        lArr = prctile(bootstrap_medians, 2.5);
        uArr = prctile(bootstrap_medians, 97.5);

    end

    lArr(isnan(lArr)) = 0;
    uArr(isnan(uArr)) = 0;
    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);
    %yscale log

end