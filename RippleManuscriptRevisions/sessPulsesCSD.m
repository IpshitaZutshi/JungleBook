function sessPulsesCSD(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',pwd,@isfolder);
parse(p,varargin{:});
expPath = p.Results.expPath;

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

%% Start collecting data
for ii = 1:size(allSess,1)
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
    
    if exist([sessionInfo.FileName '.SleepState.states.mat'],'file') 
          load([sessionInfo.FileName '.SleepState.states.mat']);
       else
          disp('No sleep states associated with this session');
    end
    load([sessionInfo.FileName '.session.mat']);
    load([sessionInfo.FileName '.pulses.events.mat']);
    
    pulTr = pulses.stimComb==2;
    homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
    homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
    pulTr = pulTr & homeCagePulseidx;
    events = pulses.intsPeriods(1,pulTr)';
    
    eventsState = events(InIntervals(events,SleepState.ints.NREMstate));
    
    file = dir(('*.hippocampalLayers.channelinfo.mat'));
    load(file.name);  
    pyrCh = hippocampalLayers.pyramidal;    
    % Only select channels from the shank that includes the pyramidal
    % channel
    % Get channels from pyr-11 ch to pyr+38
    for ch = 1:size(sessionInfo.AnatGrps,2)
%         if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
%             startCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)-11); % changed
%             endCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)+38); % changed
%             channel_order = sessionInfo.AnatGrps(ch).Channels(startCh:endCh); % changed
%         end
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                channel_order = sessionInfo.AnatGrps(ch).Channels;
            end
    end               

    % If there's a reference channel, define
    if isfield(session.channelTags,'RippleNoise')
        refChannel = session.channelTags.RippleNoise.channels-1;
    else
        refChannel = [];
    end      

    %% Load LFP
    lfp = bz_GetLFP(channel_order,'noPrompts', true);
    % Correct noise and interpolate broken channels
    lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
                
%    lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:numel(sessionInfo.AnatGrps(jj).Channels)-mod(numel(sessionInfo.AnatGrps(jj).Channels),8)),'noPrompts', true);
    twin = 0.4;
    [csd,lfpAvg] = bz_eventCSD(lfp,eventsState,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
      
    csdAv(:,:,ii) = csd.data;
    lfpAv(:,:,ii) = lfpAvg.data;
end

acrossSessCSD = nanmean(csdAv,3);
acrossSessLFP = nanmean(lfpAv,3);
figure
taxis = linspace(-twin,twin,size(csd.data,1));
cmax = max(max(acrossSessCSD)); 
contourf(taxis,1:size(acrossSessCSD,2),acrossSessCSD',40,'LineColor','none');hold on;
set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel');  
colormap jet; try caxis([-cmax cmax]); end
hold on
for kk = 1:size(acrossSessLFP,2)
    plot(taxis,(acrossSessLFP(:,kk)/1000)+kk-1,'k','LineWidth',2)
end

set(gcf,'Color','w')
set(gcf,'renderer','painters')
saveas(gcf,strcat(allSess(1).folder,'\Summ\StimNREMCSD.png'));
saveas(gcf,strcat(allSess(1).folder,'\Summ\StimNREMCSD.eps'),'epsc');
saveas(gcf,strcat(allSess(1).folder,'\Summ\StimNREMCSD.fig'));
end