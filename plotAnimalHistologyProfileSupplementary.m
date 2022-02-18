function plotAnimalHistologyProfileSupplementary

animalID = {'IZ12','IZ13','IZ18','IZ24','IZ25','IZ26','IZ31',...
     'IZ21','IZ27','IZ28','IZ33','IZ20',...
     'IZ15','IZ29','IZ30','IZ32'};%,

sessionFolder = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Histology\Track\Profile\';
expDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
histfolder = 'C:\Users\ipshi\Dropbox (NYU Langone Health)\ThetaSourceFigures\Histology_grayscale\';

%% Define main figure
fig1 = figure(1);
set(gcf,'Position',[70 350 1100 800])
set(gcf,'Renderer','painters')
    
for m = 1:length(animalID)
    
    file = dir(strcat(sessionFolder,animalID{m},'_', '*.hippocampalLayers.channelInfo.mat'));
    sessionLoc = file.name;
    sessCells = strsplit(sessionLoc,'.');
    
    if strcmp(animalID{m},'IZ27')==1 || strcmp(animalID{m},'IZ28')==1 || strcmp(animalID{m},'IZ29')==1 || ...
            strcmp(animalID{m},'IZ33')==1 
        pathDir = strcat(expDir,animalID{m},'\Saline\',sessCells{1},'\');
    else
        pathDir = strcat(expDir,animalID{m},'\Final\',sessCells{1},'\');
    end
    cd(pathDir);
    load(file.name);
    sessionInfo = bz_getSessionInfo;
    file = dir(('*.session.mat'));
    load(file.name);     
    clear hippocampalLayers
    load([sessionInfo.FileName '.hippocampalLayersCSD.channelInfo.mat'])
    % Only select channels from the shank that includes the pyramidal channel
    for ch = 1:size(sessionInfo.AnatGrps,2)
        if ismember(hippocampalLayers.pyramidal, sessionInfo.AnatGrps(ch).Channels)
            channel_order = sessionInfo.AnatGrps(ch).Channels;
        end
    end               
    %Silly exception for these animals because shank was broken
    if strcmp(session.animal.name,'IZ11') || strcmp(session.animal.name,'IZ15') 
        channel_order = channel_order(1:11);
    end
    % If there's a reference channel, define
    if isfield(session.channelTags,'RippleNoise')
        refChannel = session.channelTags.RippleNoise.channels-1;
    else
        refChannel = [];
    end         
        

    %% First plot histology figure
    subplot(8,8,[(floor(m/4)*16)+(2*(rem(m,4))+1) (floor(m/4)*16)+(2*(rem(m,4))+2)]);
    currImg = imread(strcat(histfolder,animalID{m},'.png'));
    imshow(currImg)
    title(animalID{m})
    
    
    %% Get LFP
    lfp = bz_GetLFP(channel_order,'noPrompts', true);
    % Correct noise and interpolate broken channels
    lfp = bz_interpolateLFP(lfp,'refChan',refChannel);     
    nC = 1:1:length(channel_order);

    %% Plot ripple CSD   
    ripples = bz_FindRipples(pwd,hippocampalLayers.pyramidal,'noise',refChannel);
    twin = 0.1;
    [evCsd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
    
    subplot(8,8,(floor(m/4)*16)+8+2*(mod(m,4))+1);
    contourf(evCsd.timestamps,(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
    box off; colormap(jet); caxis([-60 60]);%caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);
    hold on
    for kk = 1:size(lfpAvg.data,2)
        plot(lfpAvg.timestamps,(lfpAvg.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
    end
    xs = [evCsd.timestamps(1) evCsd.timestamps(end)];
    plot(xs, [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
    %text(xs(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);

    plot(xs,[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
    %text(xs(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);

    plot(xs,[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
    %text(xs(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
    ylim([min(nC) max(nC)]);
    set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');
    axis off
%     subplot(6,15,(floor(m/5)*30)+15+3*(rem(m,5))+3);
%     csdRippleProfile = [0 mean(evCsd.data(find(evCsd.timestamps > -10 & evCsd.timestamps < 10),:)) 0];
%     plot(csdRippleProfile,nC,'color','b');
%     hold on
%     plot(xs, [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
%     text(xs(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);
% 
%     plot(xs,[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
%     text(xs(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);
% 
%     plot(xs,[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
%     text(xs(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
%     ylim([min(nC) max(nC)]);
%     set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');
%     
    %% Plot theta CSD   
    subplot(8,8,(floor(m/4)*16)+8+2*(rem(m,4))+2);
    file = dir(('*.SessionPulses.Events.mat'));
    load(file.name);    
    file = dir(('*.SessionArmChoice.Events.mat'));
    load(file.name);      
    efields = fieldnames(sessionPulses); 
    rewardTS = sessionArmChoice.(efields{1}).timestamps;
    startDelay = sessionArmChoice.(efields{1}).delay.timestamps(1,:)';     
    endDelay = sessionArmChoice.(efields{1}).delay.timestamps(2,:)';      
    startTS = endDelay(sessionPulses.(efields{1}).stim(1:(end-1))==0);        
    endTS = rewardTS(find(sessionPulses.(efields{1}).stim(1:(end-1))==0)+1);
    events = [startTS';endTS'];
    [csd, lfpAvg] = bz_eventThetaCSD2(lfp,events,'twin',[0.2 0.2],'plotLFP',false,'plotCSD',false,'thetaCh',hippocampalLayers.pyramidal);
    
    %plot the figure
    taxis = linspace(-0.2,0.2,size(csd.data,1));
    contourf(taxis,(nC(2:end-1)),csd.data',40,'LineColor','none');hold on;
    box off; set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); 
    cmax = max(max(abs(csd.data)));
    colormap jet; caxis([-60 60]);%try caxis([-cmax cmax]); end
    hold on
    for kk = 1:size(lfpAvg.data,2)
        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
    end
    ylim([min(nC) max(nC)]);
    
    ax = axis;
    plot(ax(1:2), [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
    %text(ax(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);

    plot(ax(1:2),[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
    %text(ax(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);

    plot(ax(1:2),[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
    %text(ax(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
    axis off
    
%     subplot(6,15,(floor(m/5)*30)+15+3*(rem(m,5))+3);
%     csdThetaProfile = [0 mean(csd.data(find(csd.timestamps > -10 & csd.timestamps < 10),:)) 0];
%     plot(csdThetaProfile,nC,'color','m');
%     hold on
%     plot(xs, [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
%     text(xs(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);
% 
%     plot(xs,[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
%     text(xs(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);
% 
%     plot(xs,[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
%     text(xs(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
%     ylim([min(nC) max(nC)]);
%     set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');   
%     
    
end
end