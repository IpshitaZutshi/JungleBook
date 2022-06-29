function SessRippleBehaviorCSD(varargin)

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



if exist('Summ\csdRippleBehaviorData.mat','file') && ~force 
    disp('Ripple CSD already computed! Loading file.');
    load('Summ\csdRippleBehaviorData.mat');
else
    for zz = 1:2
        csdData{zz} = [];
        lfpAvg{zz} = [];
    end        
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); 
        file = dir(('*.session.mat'));
        load(file.name);                 
        file = dir(('*.pulses.events.mat'));
        load(file.name);         
        file = dir(('*.hippocampalLayers.channelinfo.mat'));
        load(file.name);  
        file = dir(('*.region.mat'));
        load(file.name);          
        pyrCh = hippocampalLayers.pyramidal;
        %pyrCh = region.CA1so; 

        file = dir(('*.ripples.events.mat'));
        load(file.name);   
  
        % Only select channels from the shank that includes the pyramidal
        % channel
        % Get channels from pyr-11 ch to pyr+38
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
%                 startCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)-11); % changed
%                 endCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)+38); % changed
%                 channel_order = sessionInfo.AnatGrps(ch).Channels(startCh:endCh); % changed
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

        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);               
        
        %% Load LFP
        lfp = bz_GetLFP(channel_order,'noPrompts', true);
        % Correct noise and interpolate broken channels
        lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            
            keepIdx = [];
            
            if (region ~= 2) || (target ~=2)
                continue
            end
            
            for zz = 1:2
                %Extract relevant intervals for ripples
                switch zz
                    case 1  %First, no stim trials, return        
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        events = [startTS'; endTS'];
                    case 2  % Stim, return
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        events = [startTS';endTS'];                    
                end
                
                keepIdx(:,zz) = InIntervals(ripples.peaks,events');
                
                if sum(keepIdx(:,zz))>0
                    st = ripples.peaks(logical(keepIdx(:,zz)))*1250;
                    twin = [0.2*1250 0.2*1250];
                    for e = 1:length(st)
                        lfp_temp(:,:,e) = lfp.data(st(e)-twin(1):st(e)+twin(2),:);
                    end
                    lfpAvg{zz} = cat(3,lfpAvg{zz},lfp_temp);
                end                   
                clear lfp_temp
            end
        end
                

    end
    
    if makePlots
        figure(1)
        set(gcf,'Position',[100 100 600 600])
        set(gcf,'renderer','painters')
    end
    
    for zz = 1:2
       lfp_avg = nanmean(lfpAvg{zz},3)*-1;

       if isempty(lfp_avg)
           continue
       end
        % temporal smoothing
       for ch = 1:size(lfp_avg,2) 
           lfp_avg(:,ch) = smooth(lfp_avg(:,ch),11,'sgolay');
       end

        % spatial smoothing
       for t = 1:size(lfp_avg,1) 
           lfp_avg(t,:) = smooth(lfp_avg(t,:),11,'lowess');
       end
       data = diff(lfp_avg,2,2);
       timestamps = -twin(1):twin(2);
       csdData{zz} = [0 mean(data(find(timestamps > -10 & timestamps < 30),:)) 0];

       if makePlots    
          subplot(2,2,2*(zz-1)+1)
          taxis = linspace(-twin(1),twin(2),size(data,1));
          cmax = max(max(data));      
          contourf(taxis,1:size(data,2),data',40,'LineColor','none');hold on;
          set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel');  
          colormap jet; try caxis([-cmax cmax]); end
          hold on
          for kk = 1:size(lfp_avg,2)
              plot(taxis,(lfp_avg(:,kk)/1000)+kk-1,'k')
          end
           title(num2str(size(lfpAvg{zz},3)))

          subplot(2,2,2)
          if zz == 1
              plot(csdData{zz},1:length(channel_order),'color',[85/243 85/243 85/243],'LineWidth',1.5);     
          else
              plot(csdData{zz},1:length(channel_order),'color',[8/243 133/243 161/243],'LineWidth',1.5);     
              set(gca,'YDir','reverse')
          end
          hold on       
       end
    end

    if saveMat 
        save([expPath '\Summ\' 'csdRippleBehaviorData.mat'], 'csdData');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRippleBehavior.png'),'png');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRippleBehavior.fig'),'fig');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRippleBehavior.eps'),'epsc');
    end
end
close all
end