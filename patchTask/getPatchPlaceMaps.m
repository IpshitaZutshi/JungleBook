function getPatchPlaceMaps(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveLoc',[],@isstr);
addParameter(p,'forceCalculate',[],@isstr);
addParameter(p,'plotfig',false,@islogical);

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
    firingMaps.portRun = rewarded;
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
    plotTrialMaps = 0;
    
    %% First plot the spike data plots
    if plotSpikeData
        for cellNum = 1:length(spikeData.pos)  
            figure
            set(gcf,'Color','w')
            set(gcf,'Position',[182 166 1585 762])
            plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1.2)
            hold on
            scatter(tracking.timestamps(spikeData.posIdx{cellNum}),tracking.position.y(spikeData.posIdx{cellNum}),5,'m','filled')

            for kk = 1:7
                idx = find(behavTrials.port==kk);
                for ii= 1:length(idx)
                    lickTS = behavTrials.timestamps(idx(ii));
                    [~,ixTS] = min(abs(tracking.timestamps-lickTS));    
                    scatter(tracking.timestamps(ixTS),tracking.position.y(ixTS),10,col1(kk,:),'filled')                        
                end
            end

            ylim([0 125])
            xlabel('Time(s)')
            ylabel('Position on track (cm)')
            box off 

            saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.png'],'png');
            saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.fig'],'fig');
            saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.eps'],'epsc');
            close all;
        end
    end
    
    if plotTrialMaps
    end
    
%     for pf = 1:length(firingMaps.forward.rateMaps)
%          %  First get lin track
%          figure
%          set(gcf,'Renderer','painters')
%          set(gcf,'Position',[2200 200 1185 712])
%          for ll = 1:length(labels)
%              for ii = 1:length(idx)
%                 datamat = [];
%                 numtrials = find(idx{ii}==1);
%                 for kk = 1:length(numtrials)
%                     datamat = [datamat;firingMaps.(labels{ll}).rateMaps{pf}{numtrials(kk)}];
%                 end
%                 subplot(3,13,13*(ll-1)+ii) % linear trials      
%                 if ~isempty(datamat)
%                     h = imagesc(datamat);
%                     set(h, 'AlphaData', ~isnan(datamat))
%                     
%                 else
%                     axis off
%                 end
%              end                       
%          end                  
% 
%         saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.png'],'png');
%         saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.fig'],'fig');
%         close all;
%     end    
end
end