function FigS5_AssemblySequences

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

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[226 259 1652 647]);

for ii = 1:14
    PSTHAssemblies{ii} = [];
end

% for ss = 1:length(sess)
%     try
%         %% Load files
%         cd(strcat(expPath,sess{ss})) 
% 
%         file = dir(['*spikes.cellinfo.mat']);
%         load(file(1).name);
%         file = dir(['*cell_metrics.cellinfo.mat']);
%         load(file.name);
%         file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
%         load(file(1).name);
%         file = dir(['*TrialBehavior.Behavior.mat']);
%         load(file(1).name);
%         file = dir(['*Tracking.Behavior.mat']);
%         load(file(1).name);
%         [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
% 
%         dtime = mean(diff(tracking.timestamps));
%         spkData = bz_SpktToSpkmat(spikes.times,'dt',dtime, 'win',[behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
%         spikeMat = spkData.data';
% 
%         %Only select pyramidal cells
%         logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
%         keepCells = logicalVector;
%         spikeMat = spikeMat(keepCells,:);
%         cellId = find(keepCells);
% 
%         SpikeCount = zscore(spikeMat,[],1);
%         AssemblyTemplates = assembly_patterns(SpikeCount);
% 
%         Vectors = [];
%         assembliesID = [];
%         assembliesIDNeg = [];
%         for aa = 1:size(AssemblyTemplates,2)
%             assembliesVector = AssemblyTemplates(:,aa);
%             % flip if negative weights:
%             assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
%                 /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
%             Vectors(:,aa) = assembliesVector;
% 
%             %For each assembly, define significant cells
%             sigCells = find(assembliesVector> 2*std(assembliesVector));
%             assembliesID{aa} = cellId(sigCells);
% 
%             sigCells = find(assembliesVector< -2*std(assembliesVector));
%             assembliesIDNeg{aa} = cellId(sigCells);
% 
%         end  
% 
%         Activities = assembly_activity(Vectors,SpikeCount);
% 
%         %% Calculate PSTH of the assemblies
%         for ii = 1:14
%             assemblyTimes = [];
%             psth1 = [];
%             for aa = 1:size(Activities,1)
%                 idx = find(Activities(aa,:) > 1.5*std(Activities(aa,:)));
%                 assemblyTimes{aa} = spkData.timestamps(idx);
%                 if ii< 7
%                     st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.lickLoc==(ii-1) & behavTrials.correct==1,2);                   
%                 elseif ii == 7
%                     st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.correct==1,2);
%                 elseif ii>7 && ii <14
%                     st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.lickLoc==(ii-8) & behavTrials.correct==0,2);
%                 else
%                     st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.correct==0,2);
%                 end
% 
%                 if ~isempty(st)
%                     [stccg, t] = CCG({assemblyTimes{aa} st},[],'binSize',0.1,'duration',10,'norm','rate'); 
%                     psth1(aa,:) = stccg(:,2,1);  
%                 else
%                     fillArr = nan(1,101);
%                     psth1(aa,:) = fillArr;
%                 end
%             end
%             PSTHAssemblies{ii} = [PSTHAssemblies{ii}; psth1];
%         end
% 
%     catch
%     end
% end

%save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Assemblysequences.mat','PSTHAssemblies')

load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Assemblysequences.mat')

%% Plot figure
t = linspace(-5,5,101);

[~,maxRateIdx] = max(PSTHAssemblies{7},[],2);
[maxRate,idxmax] = sort(maxRateIdx,'descend'); 

RdPu=cbrewer('seq', 'Greys', 11);

for ii = 1:14
    ax1 = subplot(3,7,ii);
    norm = zscore(PSTHAssemblies{ii},[],2);
    imagesc(t,1:length(idxmax),norm(idxmax,:))
    set(gca,'YDir','normal')
    colormap(ax1,RdPu)
    caxis([-1 4])
    colorbar
    hold on
end

col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];


for ii = 1:6
    latencies = [];
    curPSTH = PSTHAssemblies{ii}(idxmax,:);

    [~, peakIdx] = max(curPSTH,[],2);
    selidx = abs(peakIdx-maxRate)>30;

    latencies = peakIdx;
    latencies(selidx) = nan;

    X = (1:length(latencies))';
    y = latencies;
    b = regress(X, [ones(length(latencies),1) y]);
    %b = robustfit(y,X);
    
    % Store the slope
    slopes(ii) = b(2);  % b(2) is the slope

    hold off
    subplot(3,7,14+ii)
    plot(y, X, '.', 'Color', col(ii, :));
    hold on
    %lsline
    plot(y,b(1) + b(2)*y, '-', 'Color', col(ii, :), 'LineWidth', 2); 
    box off
    title(num2str(slopes(ii)))
end

subplot(3,7,21)
bar(slopes,'EdgeColor','none')
box off
%ylim([-0.12 0])


%Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure5A_AssemblySequences.png'));
saveas(gcf,strcat(expPath,'SupFigure5A_AssemblySequences.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure5A_AssemblySequences.fig'));
end