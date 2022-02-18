function SessBehaviorCSDCoherence(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSessMat',true,@islogical);
addParameter(p,'makePlots',true,@islogical);
addParameter(p,'plotfig',false,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'twin',0.2,@isnumeric)
parse(p,varargin{:});

expPath = p.Results.expPath;
makePlots = p.Results.makePlots;
plotfig = p.Results.plotfig;
saveMat = p.Results.saveMat;
saveSessMat = p.Results.saveSessMat;
force = p.Results.force;
twin = p.Results.twin;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

zonename = {'Return','Stem','Delay'};
if exist('Summ\csdcohData.mat','file') && ~force 
    disp('Theta CSD already computed! Loading file.');
    load('Summ\csdcohData.mat');
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                csdcohData.csd{rr,cc}{zz} = [];
                csdcohData.cohTheta{rr,cc}{zz} = [];
                csdcohData.cohPhaseTheta{rr,cc}{zz} = [];
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
        
        file = dir(('*.region.mat'));
        load(file.name);
%         if pyrCh == region.CA1sp
%             disp(strcat('Pyr channel has not changed, skipping session ',num2str(ii),' of ',num2str(size(allSess,1))))
%             skipSaveMat = 1;
%             continue;
%         end
%         
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
        continueAnalysis(1:length(efields)) = 1;
        
        for jj = 1:length(efields)
            if exist([efields{jj} '\Analysis\CSD_Coherence.mat'],'file') && ~force 
                disp('CSD Coh already computed for this session! Loading file.');
                load([efields{jj} '\Analysis\CSD_Coherence.mat']);
                continueAnalysis(jj) = 0;
                region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
                target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
                for zz = 1:6  
                    %Populate CSD
                    csdProfile(zz,:) = [0 mean(csd{zz}.data(find(csd{zz}.timestamps > -10 & csd{zz}.timestamps < 10),:)) 0];
                    csdcohData.csd{region,target}{zz} =  cat(3,csdcohData.csd{region,target}{zz},csdProfile(zz,:));

                    %Populate coherence
                    fRel = find(coherence{zz}.frequency >=6 & coherence{zz}.frequency <=12);
                    cohTheta(zz,:) = mean(coherence{zz}.avgCoherogram(:,fRel),2);
                    cohPhaseTheta(zz,:) = circ_mean(coherence{zz}.avgPhase(:,fRel),[],2);
                    csdcohData.cohTheta{region,target}{zz} =  cat(3,csdcohData.cohTheta{region,target}{zz},cohTheta(zz,:));
                    csdcohData.cohPhaseTheta{region,target}{zz} =  cat(3,csdcohData.cohPhaseTheta{region,target}{zz},cohPhaseTheta(zz,:));
                end  
            clear coherence csd lfpAvg   
            end
        end
                
        %% Load LFP
        if sum(continueAnalysis)>0
            lfp = bz_GetLFP(channel_order,'noPrompts', true);
            % Correct noise and interpolate broken channels
            lfp = bz_interpolateLFP(lfp,'refChan',refChannel);

            coherenceAll = bz_lfpCoherence(orCh,lfp,'forceDetect',true,'showfig',false,'saveMat',true);         

            for jj = 1:length(efields)
                if continueAnalysis(jj) ==1
                    region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
                    target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

                    rewardTS = sessionArmChoice.(efields{jj}).timestamps;
                    startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
                    endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

                    for zz = 1:6
                        %Extract relevant intervals for cross-frequency coupling - 4 cross
                        %modulograms
                        switch zz
                            case 1  %First, no stim trials, return        
                                startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                                endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                                events = [startTS'; endTS'];
                            case 2  %No stim, stem
                                startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                                endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
                                events = [startTS';endTS'];
                            case 3 %No stim, delay
                                startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                                endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
                                events = [startTS';endTS'];  
                            case 4  % Stim, return
                                startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                                endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                                events = [startTS';endTS'];                    
                            case 5   % Stim, stem
                                startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                                endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
                                events = [startTS';endTS'];                      
                            case 6    %stim, delay
                                startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                                endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
                                events = [startTS';endTS'];
                        end

                        [csd{zz}, lfpAvg{zz}] = bz_eventThetaCSD2(lfp,events,'twin',[twin twin],'plotLFP',false,'plotCSD',false,'thetaCh',pyrCh);

                        %Extract relevant timestamps
                        keeptimes = InIntervals(coherenceAll.timebins,events');

                        % Only keep relevant interval data
                        coherence{zz}.frequency = coherenceAll.frequency;
                        coherence{zz}.coherogram(:,:,:) = coherenceAll.coherogram(:,keeptimes,:);
                        coherence{zz}.phase(:,:,:) = coherenceAll.phase(:,keeptimes,:);
                        timestamps = coherenceAll.timebins(keeptimes);
                        coherence{zz}.timebins = timestamps;

                        %Average across time to show along channels
                        avgCoherogram = mean(coherence{zz}.coherogram,2);
                        avgCoherogram = reshape(avgCoherogram,length(coherence{zz}.frequency),length(lfp.channels));
                        coherence{zz}.avgCoherogram = avgCoherogram';

                        avgPhase = circ_mean(coherence{zz}.phase,[],2);
                        avgPhase = reshape(avgPhase,length(coherence{zz}.frequency),length(lfp.channels));
                        coherence{zz}.avgPhase = avgPhase';
                    end

                    for zz = 1:6  
                        %Populate CSD
                        csdProfile(zz,:) = [0 mean(csd{zz}.data(find(csd{zz}.timestamps > -10 & csd{zz}.timestamps < 10),:)) 0];
                        csdcohData.csd{region,target}{zz} =  cat(3,csdcohData.csd{region,target}{zz},csdProfile(zz,:));

                        %Populate coherence
                        fRel = find(coherence{zz}.frequency >=6 & coherence{zz}.frequency <=12);
                        cohTheta(zz,:) = mean(coherence{zz}.avgCoherogram(:,fRel),2);
                        cohPhaseTheta(zz,:) = circ_mean(coherence{zz}.avgPhase(:,fRel),[],2);
                        csdcohData.cohTheta{region,target}{zz} =  cat(3,csdcohData.cohTheta{region,target}{zz},cohTheta(zz,:));
                        csdcohData.cohPhaseTheta{region,target}{zz} =  cat(3,csdcohData.cohPhaseTheta{region,target}{zz},cohPhaseTheta(zz,:));
                    end

                    if makePlots   
                        figure(1)
                        set(gcf,'Position',[100 100 600 600])
                        set(gcf,'renderer','painters')

                        figure(2)
                        set(gcf,'Position',[100 100 600 600])
                        set(gcf,'renderer','painters')

                        figure(3)
                        set(gcf,'Position',[100 100 600 600])
                        set(gcf,'renderer','painters')

                        for zz = 1:6
                            %Plot the CSD                    
                            figure(1)
                            subplot(3,3,zz);           
                            taxis = linspace(-twin,twin,size(csd{zz}.data,1));
                            contourf(taxis,1:size(csd{zz}.data,2),csd{zz}.data',40,'LineColor','none');hold on;
                            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); 
                            colormap jet; caxis([-150 150]); 
                            colorbar
                            hold on
                            for kk = 1:size(lfpAvg{zz}.data,2)
                                plot(taxis,(lfpAvg{zz}.data(:,kk)/1000)+kk-1,'color',[0.2 0.2 0.2])
                            end         
                            if zz<4
                                subplot(3,3,zz+6)
                                col = [85/243 85/243 85/243];  
                            else
                                subplot(3,3,zz+3)
                                col = [8/243 133/243 161/243];
                            end
                            plot(csdProfile(zz,:),1:length(channel_order),'color',col,'LineWidth',1.5); 
                            set(gca,'YDir','reverse'); xlabel('CSD power'); ylabel('channel'); 
                            xlim([-120 120])
                            hold on
                            if zz<4
                               title(zonename{zz})
                            end

                            %Plot the Coherence                    
                            figure(2)
                            subplot(3,3,zz);                         
                            imagesc(coherence{zz}.frequency,1:length(channel_order),coherence{zz}.avgCoherogram)
                            xlabel('Frequency(Hz)')
                            ylabel('Channel')
                            title('Coherogram Amplitude');
                            if zz<4
                                subplot(3,3,zz+6)
                                col = [85/243 85/243 85/243];  
                            else
                                subplot(3,3,zz+3)
                                col = [8/243 133/243 161/243];
                            end
                            plot(cohTheta(zz,:),1:length(channel_order),'color',col,'LineWidth',1.5); 
                            set(gca,'YDir','reverse'); xlabel('Theta coherence'); ylabel('channel'); 
                            xlim([0 1])
                            hold on
                            if zz<4
                               title(zonename{zz})
                            end                    

                            %Plot the Coherence phase    
                            figure(3)
                            subplot(3,3,zz);                         
                            imagesc(coherence{zz}.frequency,1:length(channel_order),coherence{zz}.avgPhase)
                            xlabel('Frequency(Hz)')
                            ylabel('Channel')
                            title('Coherogram Phase');
                            if zz<4
                                subplot(3,3,zz+6)
                                col = [85/243 85/243 85/243];  
                            else
                                subplot(3,3,zz+3)
                                col = [8/243 133/243 161/243];
                            end
                            plot(rad2deg(cohPhaseTheta(zz,:)),1:length(channel_order),'color',col,'LineWidth',1.5); 
                            set(gca,'YDir','reverse'); xlabel('Theta coherence phase'); ylabel('channel'); 
                            xlim([-180 180])
                            hold on
                            if zz<4
                               title(zonename{zz})
                            end                       

                        end

                        if ~isfolder([efields{jj},'\Analysis'])
                            mkdir(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}),'Analysis')
                        end
                        saveas(figure(1),[efields{jj} '\Analysis\CSD.png'],'png');
                        saveas(figure(1),[efields{jj} '\Analysis\CSD.fig'],'fig');
                        saveas(figure(1),[efields{jj} '\Analysis\CSD.eps'],'epsc');

                        saveas(figure(2),[efields{jj} '\Analysis\Coherence.png'],'png');
                        saveas(figure(2),[efields{jj} '\Analysis\Coherence.fig'],'fig');
                        saveas(figure(2),[efields{jj} '\Analysis\Coherence.eps'],'epsc');

                        saveas(figure(3),[efields{jj} '\Analysis\CohPhase.png'],'png');
                        saveas(figure(3),[efields{jj} '\Analysis\CohPhase.fig'],'fig');
                        saveas(figure(3),[efields{jj} '\Analysis\CohPhase.eps'],'epsc');             
                    end

                    if saveSessMat
                        if ~isfolder([efields{jj},'\Analysis'])
                            mkdir(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}),'Analysis')
                        end
                        save([efields{jj} '\Analysis\CSD_Coherence.mat'], 'csd','lfpAvg','coherence');
                    end
                    close all
                    clear rewardTS startDelay events csd lfpAvg coherence
                end
            end
        end
    end
    if saveMat 
        save([expPath '\Summ\' 'csdcohData.mat'], 'csdcohData');
    end
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