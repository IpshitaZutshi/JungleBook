function plotaccelerationbyPort

sesstoAnalyze = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28',... 
    }; 

plotfig = 1;
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';

col = [0 0 0;...
    0.5 0.5 0.5;...
    119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name);
    
    file = dir('*.Tracking.Behavior.mat');
    load(file.name);

    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.correct),1) = 0;
    end

    if ~isfield(behavTrials,'stim')
        behavTrials.stim(1:length(behavTrials.correct),1) = 0;
    end
    
    if strcmp(sesstoAnalyze{ii}(1:4),'IZ39')==1
        Summary.mouseID(ii) = 1;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ40')==1
        Summary.mouseID(ii) = 2;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ43')==1
        Summary.mouseID(ii) = 3;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ44')==1
        Summary.mouseID(ii) = 4;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ47')==1
        Summary.mouseID(ii) = 5;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ48')==1
        Summary.mouseID(ii) = 6;
    end

    %Get the index for different conditions
    linIdx = find(behavTrials.linTrial==1);
    jumpLin = find(diff(linIdx)>1);
    
    % Group trials
    if isempty(jumpLin)
        idx{1} = linIdx; % Idx of initial runs
        idx{2} = [];
    else
        idx{1} = linIdx(1:jumpLin); % Idx of initial runs
        idx{2} = linIdx(jumpLin:end); % Idx of later runs
    end
    idx{3} = find(behavTrials.linTrial==0 & behavTrials.probe==0 & behavTrials.stim==0);
   
    % Extract a trial by trial acceleration for each position bin,
    % separated for NT1, NT2 and TONE trials
    for tt = 3:8
        sessSum{ii}{tt}=[];
    end
    % 
    % if plotfig
    %     figure
    % end

    for tt = 1:3
        for trial = 1:length(idx{tt})
            lickloc = behavTrials.lickLoc(idx{tt}(trial),:);

            [idxPos] = InIntervals(tracking.timestamps,behavTrials.timestamps(idx{tt}(trial),:));            
            yPos = tracking.position.y(idxPos);
            vPos = tracking.position.v(idxPos);
            accPos = gradient(vPos)./(1/tracking.samplingRate);

            %Only consider positions below 8,  
            binEdges = 5:1:122;

            % Use histcounts to bin the yPos data
            [binCounts,~, binIndices] = histcounts(yPos, binEdges);
            % Initialize arrays to store the average acceleration in each bin
            avgAccel = zeros(size(binEdges)-1);
            
            % Calculate average acceleration in each bin
            for i = 1:(length(binEdges)-1)
                indicesInBin = binIndices == i;
                avgAccel(i) = mean(accPos(indicesInBin));
            end     

            if tt == 3
                sessSum{ii}{lickloc+3} = [sessSum{ii}{lickloc+3}; avgAccel];
            else
                sessSum{ii}{tt}(trial,:) = avgAccel;
            end

            % if plotfig
            %     if tt ==3
            %         idloc = lickloc+3;
            %     else
            %         idloc = tt;
            %     end
            %     subplot(1,8,idloc)
            %     plot(binEdges(2:end),avgAccel,'color',col(idloc,:))
            %     hold on
            %     xlim([binEdges(2) binEdges(end)])
            % end            
        end
    end

    for tt = 1:8
        if ~isempty(sessSum{ii}{tt})
            sessTally{tt}(ii,:) = nanmedian(sessSum{ii}{tt},1);
        else
            sessTally{tt}(ii,:) = nan(1,length(binEdges)-1);
        end

    end
    if plotfig
        figure
        for tt = 1:8
            if ~isempty(sessSum{ii}{tt})
                plot(binEdges(2:end),nanmedian(sessSum{ii}{tt},1),'Color',col(tt,:))
                hold on
                xlim([binEdges(2) binEdges(end)])
            end
        end
    end
end
end