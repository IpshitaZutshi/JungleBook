function SessBehaviorCoherenceCSD(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSessMat',true,@islogical);
addParameter(p,'makePlots',true,@islogical);
addParameter(p,'plotfig',false,@islogical);
addParameter(p,'force',true,@islogical);
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
if exist('Summ\cohDataCSD.mat','file') && ~force 
    disp('CSD coherence already computed! Loading file.');
    load('Summ\cohDataCSD.mat');
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
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
        
        file = dir(('*.hippocampalLayersCSD.channelinfo.mat'));
        load(file.name);  
        radCh = hippocampalLayers.oriens;
        channel_order = hippocampalLayers.channelOrder;
        radChIdx = find(channel_order==radCh);
        if radChIdx==1
            radCh = channel_order(2);
        end
        % If there's a reference channel, define
        refChannel = [];
        
        efields = fieldnames(sessionPulses);  
        continueAnalysis(1:length(efields)) = 1;
        
        for jj = 1:length(efields)
            if exist([efields{jj} '\Analysis\Coherence.mat'],'file') && ~force 
                disp('Coh already computed for this session! Loading file.');
                load([efields{jj} '\Analysis\Coherence.mat']);
                continueAnalysis(jj) = 0;
                region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
                target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
                for zz = 1:6  
                    %Populate coherence
                    fRel = find(coherence{zz}.frequency >=6 & coherence{zz}.frequency <=12);
                    cohTheta(zz,:) = mean(coherence{zz}.avgCoherogram(:,fRel),2);
                    cohPhaseTheta(zz,:) = circ_mean(coherence{zz}.avgPhase(:,fRel),[],2);
                    csdcohData.cohTheta{region,target}{zz} =  cat(3,csdcohData.cohTheta{region,target}{zz},cohTheta(zz,:));
                    csdcohData.cohPhaseTheta{region,target}{zz} =  cat(3,csdcohData.cohPhaseTheta{region,target}{zz},cohPhaseTheta(zz,:));
                end  
            clear coherence
            end
        end
                
        %% Load LFP
        if sum(continueAnalysis)>0
            lfp = bz_GetLFP(channel_order,'noPrompts', true);
            % Correct noise and interpolate broken channels
            lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
            csd  = bz_CSDIZ(lfp,'plotCSD',false,'plotLFP',false);
            csd.channels = lfp.channels(2:(end-1));
            
            coherenceAll = bz_lfpCoherence(radCh,csd,'forceDetect',true,'showfig',false,'saveMat',false);         

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
                        avgCoherogram = reshape(avgCoherogram,length(coherence{zz}.frequency),length(csd.channels));
                        coherence{zz}.avgCoherogram = avgCoherogram';

                        avgPhase = circ_mean(coherence{zz}.phase,[],2);
                        avgPhase = reshape(avgPhase,length(coherence{zz}.frequency),length(csd.channels));
                        coherence{zz}.avgPhase = avgPhase';
                    end

                    for zz = 1:6  
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

                        for zz = 1:6
 
                            %Plot the Coherence                    
                            figure(1)
                            subplot(3,3,zz);                         
                            imagesc(coherence{zz}.frequency,1:length(channel_order)-2,coherence{zz}.avgCoherogram)
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
                            plot(cohTheta(zz,:),1:length(channel_order)-2,'color',col,'LineWidth',1.5); 
                            set(gca,'YDir','reverse'); xlabel('Theta coherence'); ylabel('channel'); 
                            xlim([0 1])
                            hold on
                            if zz<4
                               title(zonename{zz})
                            end                    

                            %Plot the Coherence phase    
                            figure(2)
                            subplot(3,3,zz);                         
                            imagesc(coherence{zz}.frequency,1:length(channel_order)-2,coherence{zz}.avgPhase)
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
                            plot(rad2deg(cohPhaseTheta(zz,:)),1:length(channel_order)-2,'color',col,'LineWidth',1.5); 
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

                        saveas(figure(1),[efields{jj} '\Analysis\CSDCoherence.png'],'png');
                        saveas(figure(1),[efields{jj} '\Analysis\CSDCoherence.fig'],'fig');
                        saveas(figure(1),[efields{jj} '\Analysis\RadCoherence.eps'],'epsc');

                        saveas(figure(2),[efields{jj} '\Analysis\CSDCohPhase.png'],'png');
                        saveas(figure(2),[efields{jj} '\Analysis\CSDCohPhase.fig'],'fig');
                        saveas(figure(2),[efields{jj} '\Analysis\CSDCohPhase.eps'],'epsc');             
                    end

                    if saveSessMat
                        if ~isfolder([efields{jj},'\Analysis'])
                            mkdir(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}),'Analysis')
                        end
                        save([efields{jj} '\Analysis\CSDCoherence.mat'], 'coherence');
                    end
                    close all
                    clear rewardTS startDelay events coherence
                end
            end
        end
    end
    if saveMat 
        save([expPath '\Summ\' 'cohDataCSD.mat'], 'csdcohData');
    end
end
end

