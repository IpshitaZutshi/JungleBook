function plotOptomECResponse

animalID = {'IZ20','IZ21','IZ23','IZ24','IZ25','IZ26',...
    'IZ27','IZ28','IZ31','IZ32','IZ33','IZ34'};
lightOut =[6.9 4.8 9.8 9.4 12.4 5 12.2 8.5 11.6 13.5 11.5 2];
sessionFolder = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Histology\Track\Profile\';
expDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

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
        
    %% Define main figure
    fig1 = figure(1);
%     tcl = tiledlayout(1,5);
%     tcl.TileSpacing = 'compact';
%     tcl.Padding = 'compact';
    set(gcf,'Position',[70 350 1700 400])
    set(gcf,'Renderer','painters')

    %% First load histology figure
    f1 = openfig(strcat(histologyLoc,animalID{m},'_new.fig'));
    ax = gca;
    figure(1)
    s1 = subplot(1,6,[1 2],'parent',fig1);
    axcp = copyobj(ax,fig1);
    set(axcp,'Position',get(s1,'position'));
    axis off
    delete(s1);
    close(f1)
    ax = axis;
    text(ax(2)/2,-100,animalID{m},'color',[.1 .1 .1],'FontWeight','bold');
    
    %% Get LFP
    lfp = bz_GetLFP(channel_order,'noPrompts', true);
    % Correct noise and interpolate broken channels
    lfp = bz_interpolateLFP(lfp,'refChan',refChannel);     
    
    %% Plot theta power/gamma/hfo
    figure(1)
    subplot(1,6,3)
    powerProfile_theta = bz_PowerSpectrumProfile([5 12],'channels',0:(length(channel_order)-1),'lfp',lfp,'showfig',false,'saveMat',false); 
    powerProfile_mg = bz_PowerSpectrumProfile([45 120],'channels',0:(length(channel_order)-1),'lfp',lfp,'showfig',false,'saveMat',false,'forceDetect',true);
    powerProfile_hfo = bz_PowerSpectrumProfile([120 240],'channels',0:(length(channel_order)-1),'lfp',lfp,'showfig',false,'saveMat',false); 
    
    hold on
    dev1 = powerProfile_theta.mean - powerProfile_theta.std;
    dev2 = powerProfile_theta.mean + powerProfile_theta.std;   
    nC = 1:1:length(powerProfile_theta.channels);
    hold on
    fill([dev1 flip(dev2)],[nC flip(nC)],[.8 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
    p1 = plot(powerProfile_theta.mean,nC,'color',[.8 .2 .2]);

    dev1 = powerProfile_hfo.mean - powerProfile_hfo.std;
    dev2 = powerProfile_hfo.mean + powerProfile_hfo.std;
    fill([dev1 flip(dev2)]+ 10,[nC flip(nC)],[.2 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
    p2 = plot(powerProfile_hfo.mean + 10,nC,'color',[.2 .2 .2]);
    ylim([min(nC) max(nC)]);

    dev1 = powerProfile_mg.mean - powerProfile_mg.std;
    dev2 = powerProfile_mg.mean + powerProfile_mg.std;
    fill([dev1 flip(dev2)]+ 10,[nC flip(nC)],[.2 .2 .6],'FaceAlpha',.2,'EdgeColor','none');
    p3 = plot(powerProfile_mg.mean + 10,nC,'color',[.2 .2 .6]);
    ylim([min(nC) max(nC)]);
    
    ax = axis;
    plot([10 ax(2)],[find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1]);
    text(ax(2),find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);
    
    plot([10 ax(2)],[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1]);
    text(ax(2),find(channel_order==hippocampalLayers.radiatum),'Rad','color',[.5 .5 .1]);

    plot([10 ax(2)],[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1]);
    text(ax(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);

    legend([p1 p2 p3], '[5-12Hz]', '[120-240Hz]','[45-120Hz]','Location','southwest');
    set(gca,'TickDir','out'); set(gca,'YDir','rev'); xlabel('dB'); ylabel('Channels');

    
    %% Plot ripple CSD   
    ripples = bz_FindRipples(pwd,hippocampalLayers.pyramidal,'noise',refChannel);
    twin = 0.1;
    [evCsd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
    subplot(1,6,4)
    contourf(evCsd.timestamps,(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
    box off; colormap(jet); caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);
    hold on
    for kk = 1:size(lfpAvg.data,2)
        plot(lfpAvg.timestamps,(lfpAvg.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
    end
    xs = [evCsd.timestamps(1) evCsd.timestamps(end)];
    plot(xs, [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
    text(xs(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);

    plot(xs,[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
    text(xs(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);

    plot(xs,[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
    text(xs(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
    ylim([min(nC) max(nC)]);
    set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');
    
    %% Plot theta CSD   
    subplot(1,6,5)
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
    colormap jet; try caxis([-cmax cmax]); end
    hold on
    for kk = 1:size(lfpAvg.data,2)
        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
    end
    ylim([min(nC) max(nC)]);
    
    ax = axis;
    plot(ax(1:2), [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
    text(ax(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);

    plot(ax(1:2),[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
    text(ax(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);

    plot(ax(1:2),[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
    text(ax(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
        
    %% Calculate mec stim CSD
    try
        subplot(1,6,6)
        [pulses] = bz_getAnalogPulsesSine;
        pulTr = (pulses.stimComb==2);
        [csd,lfpAvg] = bz_eventCSD(lfp,pulses.intsPeriods(1,pulTr)','twin',[0.5 0.5],'plotLFP',false,'plotCSD',false);
        taxis = linspace(-0.5,0.5,size(csd.data,1));
        cmax = max(max(abs(csd.data))); 
        contourf(taxis,(nC(2:end-1)),csd.data',40,'LineColor','none');hold on;
        box off; set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel');
        colormap jet; try caxis([-cmax cmax]); end
        hold on
        for kk = 1:size(lfpAvg.data,2)
            plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2],'LineWidth',0.75)
        end
        ylim([min(nC) max(nC)]);
        ax = axis;
        plot(ax(1:2), [find(channel_order==hippocampalLayers.pyramidal) find(channel_order==hippocampalLayers.pyramidal)],'color',[.8 .2 1],'LineWidth',1.5);
        text(ax(2), find(channel_order==hippocampalLayers.pyramidal),'Pyr','color',[.8 .2 1]);

        plot(ax(1:2),[find(channel_order==hippocampalLayers.radiatum) find(channel_order==hippocampalLayers.radiatum)],'color',[.5 .5 .1],'LineWidth',1.5);
        text(ax(2),find(channel_order==hippocampalLayers.radiatum), 'Rad','color',[.5 .5 .1]);

        plot(ax(1:2),[find(channel_order==hippocampalLayers.slm) find(channel_order==hippocampalLayers.slm)],'color',[.1 .8 .1],'LineWidth',1.5);
        text(ax(2),find(channel_order==hippocampalLayers.slm),'Slm','color',[.1 .8 .1]);
    catch
    end
    %% Save figures
    saveas(figure(1),strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\AnimalProfiles\',animalID{m},'.png'))
    saveas(figure(1),strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\AnimalProfiles\',animalID{m},'.fig'))
    saveas(figure(1),strcat('Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\AnimalProfiles\',animalID{m},'.eps'),'epsc')       
    
    close all
end
end