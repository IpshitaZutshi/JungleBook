function ModData = SessPeriStimModIndexCSD(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'fixChannels',true,@islogical)
addParameter(p,'refChannel',[],@isnumeric)
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
fixChannels = p.Results.fixChannels;
refChannel = p.Results.refChannel;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
phaserange = 4:0.2:15;
amprange = 25:2:300;

chName =  {'pyramidal','radiatum','slm','ml'};
   
if exist(strcat('Summ\ModIdxCSD.mat'),'file') && ~force 
    disp('ModIdx already computed! Loading file.');
    load(strcat('Summ\ModIdxCSD.mat'));
else
    for ch = 1:4 
        for rr = 1:3
            for cc = 1:2
                for zz = 1:6
                    ModData.(chName{ch}){rr,cc}{zz} = [];
                end
            end
        end
    end
    if exist(strcat('Summ\csdData.mat'),'file')
        load(strcat('Summ\csdData.mat'));
    end  

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.session.mat'));
        load(file.name);  
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        load(file.name);    
        file = dir(('*.hippocampalLayers.channelinfo.mat'));
        load(file.name);
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
            endCh = min(layerInfoIdx(3)+8,size(csdData.orCh{2,1}{2},2));
            [~,minCSD] = min(csdData.orCh{2,1}{2}(1,startCh:endCh,1));
            layerInfoIdx(3) = startCh+(minCSD-1);  

            startCh = layerInfoIdx(2)-2;
            endCh = min(layerInfoIdx(2)+2,layerInfoIdx(3)-3);
            [~,maxCSD] = max(csdData.orCh{2,1}{2}(1,startCh:endCh,1));
            layerInfoIdx(2) = startCh+(maxCSD-1);  

            %add mol layer
            startCh = layerInfoIdx(3)+5;
            endCh = layerInfoIdx(3)+15;
            if size(csdData.orCh{2,1}{2},2)==64
                [~,maxCSD] = max(csdData.orCh{2,1}{2}(1,startCh:endCh,1));
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
            chanList = channelOrder(layerInfoIdx);
        end
        
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)      
                lfp = bz_GetLFP(sessionInfo.AnatGrps(ch).Channels,'noPrompts', true);
                if fixChannels                    
                    refChannel = [];
                    lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
                end
                csd  = bz_CSDIZ(lfp,'plotCSD',false,'plotLFP',false);
                lfpPyrCh = bz_GetLFP(pyrCh,'noPrompts', true);
            end
        end    

        csd.data = csd.data(:,find(sum(lfp.channels == chanList'))-1);
        csd.channels = lfp.channels(find(sum(lfp.channels == chanList')));              
        lfp = csd;
        clear csd
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)

            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

            if exist([efields{jj} '\ModIdxCSD.mat'],'file') && ~force 
                disp(strcat('ModIdx for sess ',num2str(ii),' already computed! Loading file.'));
                load([efields{jj} '\ModIdxCSD.mat']);
            else
                for zz = [1 2 4 5]
                    %Extract relevant intervals for cross-frequency coupling - 6 cross
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
                    
                    for ch = 1:length(chName)

                        if ~isnan(chanList(ch))
                            lfpCh = lfp;
                            lfpCh.data = lfp.data(:,lfp.channels == chanList(ch));
                            lfpCh.channels = lfp.channels(lfp.channels == chanList(ch));

                            if (zz == 3 || zz == 6) && sessionArmChoice.(efields{jj}).delay.dur < 1 
                                comod.(chName{ch}){3} = nan;
                            else
                                comod.(chName{ch}){zz} = bz_ModIndex_IZ(lfpCh,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange,'useCSD',true,'lfpPyrCh',lfpPyrCh);  
                            end
                        else
                            comod.(chName{ch}){zz} = nan;
                        end
                    end
                end                
            end
            for ch = 1:length(chName)
                for zz = [1 2 4 5]            
                   ModData.(chName{ch}){region,target}{zz} =  cat(3,ModData.(chName{ch}){region,target}{zz},comod.(chName{ch}){zz});
                end
            end

            if saveMat
                save([efields{jj} '\ModIdxCSD.mat'], 'comod');
            end
            clear comod rewardTS startDelay events
        end
    end

    if saveMat
        save([expPath '\Summ\' 'ModIdxCSD.mat'], 'ModData');
    end
end

% reg = {'CA1','mEC','Both'};
% zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
% target = {'STEM', 'RETURN'};
% 
% for ch = 1:3
%     for ii = 1:size(ModData,1)
%         figure(ii)
%         set(gcf,'Position',[100 100 1700 600])
%         for jj = 1:size(ModData,2)
%             for kk = 1:length(zone)
%                 subplot(2,length(zone),length(zone)*(jj-1)+kk)  
%                 %meancomod = ModData{ii,jj}{kk}(:,:,2);
%                 meancomod = nanmean(ModData{ii,jj}{kk},3);
%                 imagesc(phaserange(2:end),amprange(2:end),meancomod)
%                 %imagesc(phaserange(2:end),log2(amprange(2:end)),meancomod)
%                 set(gca,'YDir','normal')
%                 colormap jet
%                 if kk <=3
%                     cmax = max(max(nanmean(ModData{ii,jj}{kk},3))); 
%                 else
%                     cmax = max(max(nanmean(ModData{ii,jj}{(kk-3)},3))); 
%                 end
%                 if cmax>0
%                     caxis([0 cmax])
%                 end
%                 %caxis([0 0.0013])
%                 colorbar
%                 %LogScale('y',2)
%                 xlabel('Frequency phase');
%                 ylabel('Frequency amplitude');
%                 title(strcat(target(jj),'.',zone(kk)));
%             end
%         end
%         saveas(figure(ii),strcat(expPath,'\Summ\CSDThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.png'));
%         saveas(figure(ii),strcat(expPath,'\Summ\CSDThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.eps'),'epsc');
%         saveas(figure(ii),strcat(expPath,'\Summ\CSDThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.fig'));
%     end
%     close all
% end
end