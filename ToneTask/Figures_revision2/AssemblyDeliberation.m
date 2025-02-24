function AssemblyDeliberation

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

% AssemblyDelib.Tone.Choice = [];
% AssemblyDelib.Tone.Delib = [];
% AssemblyDelib.Other.Choice = [];
% AssemblyDelib.Other.Delib = [];
% 
% for ss = 1:length(sess)
%     try
%         %% Load files
%         cd(strcat(expPath,sess{ss})) 
% 
%         file = dir(['*spikes.cellinfo.mat']);
%         load(file(1).name);
%         file = dir(['*cell_metrics.cellinfo.mat']);
%         load(file.name);
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
%         % Find deceleration points
%         Dec = findDecelerationPoints('plotfig',false);
% 
%         %% Calculate PSTH of the assemblies
%         assemblyTimes = [];
%         psth1 = [];
%         for aa = 1:size(Activities,1)
%             idx = find(Activities(aa,:) > 1.5*std(Activities(aa,:)));
%             assemblyTimes{aa} = spkData.timestamps(idx);
%             st = behavTrials.timestamps(behavTrials.linTrial==0,2);
%             if ~isempty(st)
%                 [stccg, t] = CCG({assemblyTimes{aa} st},[],'binSize',0.1,'duration',10,'norm','rate'); 
%                 psth1(aa,:) = stccg(:,2,1);  
%             else
%                 fillArr = nan(1,101);
%                 psth1(aa,:) = fillArr;
%             end
%         end
% 
%         % Now, if the peak of the assembly is between 500 ms, or between 2500 and 2000ms of the lick,
%         % calculate its psth around choice and deliberation
%         [~,maxassembly] = max(psth1,[],2); 
%         t_max = t(maxassembly);
%         % Tone cell assemblies
%         assemblyIdx_Tone = find(t_max>=-0.3 & t_max<=0);
% 
%         % Other assemblies    
%         assemblyIdx_Other = find(t_max>=-1 & t_max<=-0.7);
% 
%         %% Now, for these assemblies, calculate the PSTH around deliberations
%         delib_tone_choice = [];
%         delib_tone_delib = [];
%         for aa = 1:length(assemblyIdx_Tone)
%             st = Dec.ts(Dec.decType==1); % Choice decelerations
%             [stccg, t_delib] = CCG({assemblyTimes{assemblyIdx_Tone(aa)} st},[],'binSize',0.1,'duration',2,'norm','rate'); 
%             delib_tone_choice(aa,:) = stccg(:,2,1); 
% 
%             st = Dec.ts(Dec.decType==2); % Delib decelerations
%             [stccg, t_delib] = CCG({assemblyTimes{assemblyIdx_Tone(aa)} st},[],'binSize',0.1,'duration',2,'norm','rate'); 
%             delib_tone_delib(aa,:) = stccg(:,2,1); 
%         end
% 
%         delib_other_choice = [];
%         delib_other_delib = [];
%         for aa = 1:length(assemblyIdx_Other)
%             st = Dec.ts(Dec.decType==1); % Choice decelerations
%             [stccg, t_delib] = CCG({assemblyTimes{assemblyIdx_Other(aa)} st},[],'binSize',0.1,'duration',2,'norm','rate'); 
%             delib_other_choice(aa,:) = stccg(:,2,1); 
% 
%             st = Dec.ts(Dec.decType==2); % Delib decelerations
%             [stccg, t_delib] = CCG({assemblyTimes{assemblyIdx_Other(aa)} st},[],'binSize',0.1,'duration',2,'norm','rate'); 
%             delib_other_delib(aa,:) = stccg(:,2,1); 
%         end
% 
%         AssemblyDelib.Tone.Choice = [AssemblyDelib.Tone.Choice;delib_tone_choice];
%         AssemblyDelib.Tone.Delib = [AssemblyDelib.Tone.Delib;delib_tone_delib];
%         AssemblyDelib.Other.Choice = [AssemblyDelib.Other.Choice;delib_other_choice];
%         AssemblyDelib.Other.Delib = [AssemblyDelib.Other.Delib;delib_other_delib];
%     end
% end

%save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\AssemblyDelib.mat','AssemblyDelib')

load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\AssemblyDelib.mat')

xAxis  = linspace(-1, 1,21)';
plotAvgStd(zscore(AssemblyDelib.Tone.Choice,[],2),1,2,1,fig2,xAxis,'k')
plotAvgStd(zscore(AssemblyDelib.Other.Choice,[],2),1,2,1,fig2,xAxis,'b')

plotAvgStd(zscore(AssemblyDelib.Tone.Delib,[],2),1,2,2,fig2,xAxis,'k')
plotAvgStd(zscore(AssemblyDelib.Other.Delib,[],2),1,2,2,fig2,xAxis,'b')

% Quantify t = -0.2 ms
data1 = AssemblyDelib.Tone.Choice(:,10);
data2 = AssemblyDelib.Other.Choice(:,10);
Summary.choice = groupStats([{data1} {data2}],[],'doPlot',false);

data1 = AssemblyDelib.Tone.Delib(:,10);
data2 = AssemblyDelib.Other.Delib(:,10);
Summary.delib = groupStats([{data1} {data2}],[],'doPlot',false);

%Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\Revision2\';
saveas(gcf,strcat(expPath,'AssemblyDeliberation.png'));
saveas(gcf,strcat(expPath,'AssemblyDeliberation.eps'),'epsc');
saveas(gcf,strcat(expPath,'AssemblyDeliberation.fig'));
save(strcat(expPath,'AssemblyDeliberation.mat'),'Summary');
end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    stdpsth = nanstd(array,1)./sqrt(size(array,1));
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end