function SessBehaviorCSDUpdateAvg(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSessMat',true,@islogical);
addParameter(p,'makePlots',true,@islogical);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'twin',0.2,@isnumeric)
parse(p,varargin{:});

expPath = p.Results.expPath;
makePlots = p.Results.makePlots;
plotfig = p.Results.plotfig;
saveMat = p.Results.saveMat;
saveSessMat = p.Results.saveSessMat;
twin = p.Results.twin;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

zonename = {'Return','Stem','Delay'};

for rr = 1:3
    for cc = 1:2
        for zz = 1:6
            csdData.orCh{rr,cc}{zz} = [];
            csdData.slmCh{rr,cc}{zz} = [];
        end
    end
end

for ii = 1:size(allSess,1)
    fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    file = dir(('*.SessionPulses.Events.mat'));
    load(file.name);
    file = dir(('*.SessionArmChoice.Events.mat'));
    load(file.name);    
    file = dir(('*.session.mat'));
    load(file.name);                 
    file = dir(('*.hippocampalLayers.channelinfo.mat'));
    load(file.name);  
    pyrCh = hippocampalLayers.pyramidal; 
    orCh = hippocampalLayers.oriens;
    slmCh = hippocampalLayers.slm;

    file = dir(('*.region.mat'));
    load(file.name);

    % Only select channels from the shank that includes the pyramidal channel
    for ch = 1:size(sessionInfo.AnatGrps,2)
        if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
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
    efields = fieldnames(sessionPulses);  

%         %% Load LFP
%         lfp = bz_GetLFP(channel_order,'noPrompts', true);
%         % Correct noise and interpolate broken channels
%         lfp = bz_interpolateLFP(lfp,'refChan',refChannel);

    for jj = 1:length(efields)

        load(strcat(efields{jj},'\Analysis\CSD.mat'));

        region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
        target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

        rewardTS = sessionArmChoice.(efields{jj}).timestamps;
        startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
        endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

        for zz = 1:6

            %Populate CSD
            csdProfile.orCh(zz,:) = [0 mean(csd.orCh{zz}.data(find(csd.orCh{zz}.timestamps > -10 & csd.orCh{zz}.timestamps < 10),:)) 0];
            csdData.orCh{region,target}{zz} =  cat(3,csdData.orCh{region,target}{zz},csdProfile.orCh(zz,:));

            csdProfile.slmCh(zz,:) = [0 mean(csd.slmCh{zz}.data(find(csd.slmCh{zz}.timestamps > -10 & csd.slmCh{zz}.timestamps < 10),:)) 0];
            csdData.slmCh{region,target}{zz} =  cat(3,csdData.slmCh{region,target}{zz},csdProfile.slmCh(zz,:));
        end

        if makePlots   
            figure(1)
            set(gcf,'Position',[100 100 600 600])
            set(gcf,'renderer','painters')

            figure(2)
            set(gcf,'Position',[100 100 600 600])
            set(gcf,'renderer','painters')

            for zz = 1:6
                %Plot the CSD                    
                figure(1)
                subplot(3,3,zz);           
                taxis = linspace(-twin,twin,size(csd.orCh{zz}.data,1));
                contourf(taxis,1:size(csd.orCh{zz}.data,2),csd.orCh{zz}.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); 
                colormap jet; caxis([-150 150]); 
                colorbar
                hold on
                for kk = 1:size(lfpAvg.orCh{zz}.data,2)
                    plot(taxis,(lfpAvg.orCh{zz}.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
                end         
                if zz<4
                    subplot(3,3,zz+6)
                    col = [85/243 85/243 85/243];  
                else
                    subplot(3,3,zz+3)
                    col = [8/243 133/243 161/243];
                end
                plot(csdProfile.orCh(zz,:),1:length(channel_order),'color',col,'LineWidth',1.5); 
                set(gca,'YDir','reverse'); xlabel('CSD power'); ylabel('channel'); 
                xlim([-120 120])
                hold on
                if zz<4
                   title(zonename{zz})
                end                     

                %Plot the CSD slm                    
                figure(2)
                subplot(3,3,zz);           
                taxis = linspace(-twin,twin,size(csd.slmCh{zz}.data,1));
                contourf(taxis,1:size(csd.slmCh{zz}.data,2),csd.slmCh{zz}.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); 
                colormap jet; caxis([-150 150]); 
                colorbar
                hold on
                for kk = 1:size(lfpAvg.slmCh{zz}.data,2)
                    plot(taxis,(lfpAvg.slmCh{zz}.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
                end         
                if zz<4
                    subplot(3,3,zz+6)
                    col = [85/243 85/243 85/243];  
                else
                    subplot(3,3,zz+3)
                    col = [8/243 133/243 161/243];
                end
                plot(csdProfile.slmCh(zz,:),1:length(channel_order),'color',col,'LineWidth',1.5); 
                set(gca,'YDir','reverse'); xlabel('CSD power'); ylabel('channel'); 
                xlim([-120 120])
                hold on
                if zz<4
                   title(zonename{zz})
                end  

            end

            if ~isfolder([efields{jj},'\Analysis'])
                mkdir(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}),'Analysis')
            end
            saveas(figure(1),[efields{jj} '\Analysis\CSDor.png'],'png');
            saveas(figure(1),[efields{jj} '\Analysis\CSDor.fig'],'fig');
            saveas(figure(1),[efields{jj} '\Analysis\CSDor.eps'],'epsc');

            saveas(figure(2),[efields{jj} '\Analysis\CSDslm.png'],'png');
            saveas(figure(2),[efields{jj} '\Analysis\CSDslm.fig'],'fig');
            saveas(figure(2),[efields{jj} '\Analysis\CSDslm.eps'],'epsc');                            
        end

        if saveSessMat
            if ~isfolder([efields{jj},'\Analysis'])
                mkdir(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}),'Analysis')
            end
            save([efields{jj} '\Analysis\CSD.mat'], 'csd','lfpAvg');
        end
        close all
        clear rewardTS startDelay events csd lfpAvg
    end

    orCh = hippocampalLayers.oriens; 
    pyrCh = hippocampalLayers.pyramidal; 
    radCh = hippocampalLayers.radiatum; 
    slmCh = hippocampalLayers.slm; 
    channelOrder = hippocampalLayers.channelOrder;
    chanList = [pyrCh radCh slmCh];

    [~,layerInfoIdx] = ismember(chanList,channelOrder);        
    if length(channelOrder)>16
        %correct slm channel
        startCh = layerInfoIdx(3)-8;
        endCh = min(layerInfoIdx(3)+8,length(csdProfile.orCh(2,:)));
        [~,minCSD] = min(csdProfile.orCh(2,startCh:endCh));
        layerInfoIdx(3) = startCh+(minCSD-1);  

        startCh = layerInfoIdx(2)-2;
        endCh = min(layerInfoIdx(2)+2,layerInfoIdx(3)-3);
        [~,maxCSD] = max(csdProfile.orCh(2,startCh:endCh));
        layerInfoIdx(2) = startCh+(maxCSD-1);  

        %add mol layer
        startCh = layerInfoIdx(3)+5;
        endCh = layerInfoIdx(3)+15;
        if length(csdProfile.orCh(2,:))==64
            [~,maxCSD] = max(csdProfile.orCh(2,startCh:endCh));
            layerInfoIdx(4) = startCh+(maxCSD-1);  
        else
            layerInfoIdx(4) = nan;
        end
    else
        layerInfoIdx(4) = nan;
    end
     
    if isnan(layerInfoIdx(4))
        chanList = [channelOrder(layerInfoIdx(1:3)) nan];
    else
        chanList = channelOrder(layerInfoIdx(1:4));
    end
    hippocampalLayers.radiatum = chanList(2); 
    hippocampalLayers.slm = chanList(3); 
    if isnan(layerInfoIdx(4))
        hippocampalLayers.ml = nan; 
    else
        hippocampalLayers.ml = chanList(4); 
    end
    hippocampalLayers.all = [hippocampalLayers.oriens hippocampalLayers.pyramidal hippocampalLayers.radiatum ...
        hippocampalLayers.slm hippocampalLayers.ml];
    
    save(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',allSess(ii).name,'.hippocampalLayersCSD.channelinfo.mat'),'hippocampalLayers');
end

if saveMat 
    save([expPath '\Summ\' 'csdData.mat'], 'csdData');
end
    
animalName = strsplit(expPath,'\');
animalName = animalName{6};

reg = {'CA1','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

if plotfig       
    if strcmp(animalName,'IZ27')==1 || strcmp(animalName,'IZ28')==1 || strcmp(animalName,'IZ29')==1 ...
            || strcmp(animalName,'IZ32')==1 || strcmp(animalName,'IZ33')==1 || strcmp(animalName,'IZ34')==1
        col = [0.3498 0.3498 0.3498;...
            0.1986 0.7214 0.6310;...        
            0.9856 0.7372 0.2537;...
            0.2305 0.2510 0.6173];    

        nC = 1:size(csdcohData.csd{2,1}{1},2);
        figure
        for jj = 1:size(csdcohData.csd,2)

            for kk = 1:3
                for ii = 2:3
                    subplot(3,6,3*(jj-1)+kk)
                    Meanprofile = nanmean(csdcohData.csd{ii,jj}{kk},3);
                    stdProfile = std(csdcohData.csd{ii,jj}{kk},[],3)./sqrt(size(csdcohData.csd,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-3,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col((2*ii)-3,:),'LineWidth',1.5); 
                        xlabel('Avg csd'); ylabel('Channel id');  
                        set(gca,'YDir','reverse');
                        title(strcat(target{jj},' ',zone{kk}(1:(end-1))));
                        hold on
                    end

                    subplot(3,6,3*(jj-1)+kk)
                    Meanprofile = nanmean(csdcohData.csd{ii,jj}{kk+3},3);
                    stdProfile = std(csdcohData.csd{ii,jj}{kk+3},[],3)./sqrt(size(csdcohData.csd,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-2,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col((2*ii)-2,:),'LineWidth',1.5); 
                        xlabel('Avg csd'); ylabel('Channel id');  
                        set(gca,'YDir','reverse');
                        title(strcat(target{jj},' ',zone{kk}(1:(end-1))));
                        hold on
                    end                

                    subplot(3,6,3*(jj-1)+kk+6)
                    Meanprofile = nanmean(csdcohData.cohTheta{ii,jj}{kk},3);
                    stdProfile = std(csdcohData.cohTheta{ii,jj}{kk},[],3)./sqrt(size(csdcohData.cohTheta,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-3,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col((2*ii)-3,:),'LineWidth',1.5); 
                        xlabel('Avg coherence'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end

                    subplot(3,6,3*(jj-1)+kk+6)
                    Meanprofile = nanmean(csdcohData.cohTheta{ii,jj}{kk+3},3);
                    stdProfile = std(csdcohData.cohTheta{ii,jj}{kk+3},[],3)./sqrt(size(csdcohData.cohTheta,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-2,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col((2*ii)-2,:),'LineWidth',1.5); 
                        xlabel('Avg coherence'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end                

                    subplot(3,6,3*(jj-1)+kk+12)
                    Meanprofile = circ_mean(csdcohData.cohPhaseTheta{ii,jj}{kk},[],3);
                    stdProfile = circ_std(csdcohData.cohPhaseTheta{ii,jj}{kk},[],[],3)./sqrt(size(csdcohData.cohPhaseTheta,3));
                    if ~isempty(dev1)
                        dev1 = Meanprofile-stdProfile;
                        dev2 = Meanprofile+stdProfile;
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-3,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'.','color',col((2*ii)-3,:)); 
                        hold on
                        dev1 = Meanprofile-stdProfile+2*pi;
                        dev2 = Meanprofile+stdProfile+2*pi;
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-3,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile+2*pi,nC,'.','color',col((2*ii)-3,:)); 
                        hold on                    
                        xlabel('Avg coherence phase'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end

                    subplot(3,6,3*(jj-1)+kk+12)
                    Meanprofile = circ_mean(csdcohData.cohPhaseTheta{ii,jj}{kk+3},[],3);
                    stdProfile =circ_std(csdcohData.cohPhaseTheta{ii,jj}{kk+3},[],[],3)./sqrt(size(csdcohData.cohPhaseTheta,3));
                    if ~isempty(dev1)
                        dev1 = Meanprofile-stdProfile;
                        dev2 = Meanprofile+stdProfile;
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-2,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'.','color',col((2*ii)-2,:)); 
                        hold on
                        dev1 = Meanprofile-stdProfile+2*pi;
                        dev2 = Meanprofile+stdProfile+2*pi;
                        fill([dev1 flip(dev2)],[nC flip(nC)],col((2*ii)-2,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile+2*pi,nC,'.','color',col((2*ii)-2,:)); 
                        hold on            
                        xlabel('Avg coherence phase'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end                
                end
            end
        end

        saveas(gcf,[expPath '\Summ\' 'csdcohData.png'],'png');
        saveas(gcf,[expPath '\Summ\' 'csdcohData.fig'],'fig');
        saveas(gcf,[expPath '\Summ\' 'csdcohData.eps'],'epsc');
        close all

    else
        col = [0.9856 0.7372 0.2537;...
           0.1986 0.7214 0.6310;...
           0.2305 0.2510 0.6173];
        nC = 1:size(csdcohData.csd{2,1}{1},2);
        figure
        for jj = 1:size(csdcohData.csd,2)
            for kk = 1:3
                for ii = 1:size(csdcohData.csd,1)
                    if ii ==1
                        profile = csdcohData.csd{ii,jj}{kk};
                        coh = csdcohData.cohTheta{ii,jj}{kk};
                        phase = csdcohData.cohPhaseTheta{ii,jj}{kk};
                    else
                        profile = catpad(3,profile,csdcohData.csd{ii,jj}{kk});
                        coh = catpad(3,coh,csdcohData.cohTheta{ii,jj}{kk});
                        phase = catpad(3,phase,csdcohData.cohPhaseTheta{ii,jj}{kk});
                    end

                    subplot(3,6,3*(jj-1)+kk)
                    Meanprofile = nanmean(csdcohData.csd{ii,jj}{kk+3},3);
                    stdProfile = std(csdcohData.csd{ii,jj}{kk+3},[],3)./sqrt(size(csdcohData.csd,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5); 
                        xlabel('Avg csd'); ylabel('Channel id');  
                        set(gca,'YDir','reverse');
                        title(strcat(target{jj},' ',zone{kk}(1:(end-1))));
                        hold on
                    end

                    subplot(3,6,3*(jj-1)+kk+6)
                    Meanprofile = nanmean(csdcohData.cohTheta{ii,jj}{kk+3},3);
                    stdProfile = std(csdcohData.cohTheta{ii,jj}{kk+3},[],3)./sqrt(size(csdcohData.cohTheta,3));
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    if ~isempty(dev1)
                        fill([dev1 flip(dev2)],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5); 
                        xlabel('Avg coherence'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end
                    if ~isempty(dev1)
                        subplot(3,6,3*(jj-1)+kk+12)
                        Meanprofile = nanmean(csdcohData.cohPhaseTheta{ii,jj}{kk+3},3);
                        stdProfile = std(csdcohData.cohPhaseTheta{ii,jj}{kk+3},[],3)./sqrt(size(csdcohData.cohPhaseTheta,3));
                        dev1 = Meanprofile-stdProfile;
                        dev2 = Meanprofile+stdProfile;
                        fill([dev1 flip(dev2)],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
                        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5); 
                        xlabel('Avg coherence phase'); ylabel('Channel id'); 
                        set(gca,'YDir','reverse');
                        hold on
                    end
                end

                subplot(3,6,3*(jj-1)+kk)
                Meanprofile = nanmean(profile,3);
                stdProfile = std(profile,[],3)./sqrt(size(profile,3));
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                if ~isempty(dev1)
                    fill([dev1 flip(dev2)],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5);   
                    set(gca,'YDir','reverse');
                end

                subplot(3,6,3*(jj-1)+kk+6)
                Meanprofile = nanmean(coh,3);
                stdProfile = std(coh,[],3)./sqrt(size(coh,3));
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                if ~isempty(dev1)        
                    fill([dev1 flip(dev2)],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5); 
                    set(gca,'YDir','reverse');
                end

                subplot(3,6,3*(jj-1)+kk+12)
                Meanprofile = circ_mean(phase,[],3);
                stdProfile = circ_std(phase,[],[],3)./sqrt(size(phase,3));
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                if ~isempty(dev1)
                    fill([dev1 flip(dev2)],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5);  
                    set(gca,'YDir','reverse');
                end

            end
        end

        saveas(gcf,[expPath '\Summ\' 'csdcohData.png'],'png');
        saveas(gcf,[expPath '\Summ\' 'csdcohData.fig'],'fig');
        saveas(gcf,[expPath '\Summ\' 'csdcohData.eps'],'epsc');
        close all
    end
end
end