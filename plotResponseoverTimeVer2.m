function plotResponseoverTimeVer2(varargin)


p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'refChannel',[],@isnumeric);
addParameter(p,'pyrChPlus',[0 5 10 15 20 25 30],@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
refChannel = p.Results.refChannel;
pyrChPlus = p.Results.pyrChPlus;
%[0 1 2 3 4 5 6 7 8 9 10]
%[0 6 12 18 24 30 36 42 45 50]
%%
if ~exist('expPath') || isempty(expPath)
    expPath = pwd; % select folder
end

%% Load necessary files
cd(expPath);

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
load([sessionInfo.FileName '.MergePoints.events.mat'])
% 
% if exist([sessionInfo.FileName '.region.mat'],'file') 
%     load([sessionInfo.FileName '.region.mat']);
%     pyrCh = region.CA1sp;
% else
%     disp('First calculate .region file to identify ripple channel! Skipping');
%     return
% end
pyrCh = 74;

figure
set(gcf,'renderer','painters')
set(gcf,'Position',[100 100 900 1300])

for kk = 1:length(pyrChPlus)  
    if isempty(refChannel)
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus(kk));   
            end
        end
        lfp = bz_GetLFP(pyrChNew,'noPrompts', true);    
    else
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus(kk));   
            end
        end    
        lfp = bz_GetLFP([pyrChNew refChannel],'noPrompts', true);
        lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        lfp.data = lfp.data(:,1);
    end

    %% Calculate spectogram
    params.Fs = lfp.samplingRate; params.fpass = [2 100]; params.tapers = [3 5]; params.pad = 1;
    [S,t,f] = mtspecgramc(single(lfp.data),[2 1],params);
    S1 = log10(S); % in Db
    S_det= bsxfun(@minus,S1,polyval(polyfit(f,mean(S1,1),2),f)); % detrending

    subplot(length(pyrChPlus),1,kk)
    imagesc(t,f,S_det',[-1.5 1.5]);
    %imagesc(t,f,S1');
    colorbar
    set(gca,'YDir','normal'); ylabel('Freqs');
    xlabel('Time (s)')
    hold on
    for mm = 1:length(MergePoints.timestamps(:,2))
        xline(MergePoints.timestamps(mm,2),'w','LineWidth',1);
    end
    title(strcat('Channel number: ',num2str(pyrChNew)))
end

saveas(gcf,strcat('summ\ResponseOverTime.png'));
% saveas(gcf,strcat('summ\ResponseOverTime.eps'),'epsc');
% saveas(gcf,strcat('summ\ResponseOverTime.fig'));
end
