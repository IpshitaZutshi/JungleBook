function SessRippleCSD(varargin)

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

col = [224/243 163/243 46/243;8/243 133/243 161/243;56/243 61/243 150/243];

if exist('Summ\csdRippleData.mat','file') && ~force 
    disp('Ripple CSD already computed! Loading file.');
    load('Summ\csdRippleData.mat');
else
    for rr = 1:3
        for zz = 1:2
            csdData{rr,zz} = [];
            lfpAvg{rr,zz} = [];
        end        
    end
    
    for ii = 1:6%size(allSess,1)
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
        %pyrCh = hippocampalLayers.pyramidal;
        pyrCh = 27;%region.CA1so; 

        ripples = bz_FindRipples(pwd,pyrCh,'saveMat',false);
%         ripples = removeArtifactsFromEvents(ripples);
        
%         file = dir(('*.ripples.events.mat'));
%         load(file.name);
        pyrCh = hippocampalLayers.pyramidal;
        % Only select channels from the shank that includes the pyramidal
        % channel
        % Get channels from pyr-11 ch to pyr+38
        for ch = 1:size(sessionInfo.AnatGrps,2)
%             if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
%                 startCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)-11); % changed
%                 endCh = (find(sessionInfo.AnatGrps(ch).Channels==pyrCh)+38); % changed
%                 channel_order = sessionInfo.AnatGrps(ch).Channels(startCh:endCh); % changed
%             end
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

        %% Load LFP
        lfp = bz_GetLFP(channel_order,'noPrompts', true);
        % Correct noise and interpolate broken channels
        lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        
        for rr = 1:3
            if exist('pulses')
                if rr <= 2
                    pulTr = (pulses.stimComb==rr);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
                end
            end      
                            
            %Only select the pulses that happened in the home cage,
            %i.e., which were 5 seconds long
            homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
            homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
            pulTr = pulTr & homeCagePulseidx;
            events = pulses.intsPeriods(1,pulTr)';

            %Generate logicals for ripples in pre versus post
            ripple_pre = [];                
            ripple_post = [];               

            ripple_pre(1:length(ripples.peaks)) = 0;                
            ripple_post(1:length(ripples.peaks)) = 0;                

            for pp = 1:length(ripples.peaks)
                tempDiff = ripples.peaks(pp) - events;

                if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                   [~,idxmin] =  min(abs(tempDiff));
                   if tempDiff(idxmin) > 0
                       ripple_post(pp) = 1;
                   elseif tempDiff(idxmin) < 0
                       ripple_pre(pp) = 1;
                   end
                else
                    continue
                end
            end      

            ripple_logical = [logical(ripple_pre)' logical(ripple_post)'];
            for zz = 1:2
                if sum(ripple_logical(:,zz))>0
                    st = ripples.peaks(ripple_logical(:,zz))*1250;
                    twin = [0.2*1250 0.2*1250];
                    for e = 1:length(st)
                        lfp_temp(:,:,e) = lfp.data(st(e)-twin(1):st(e)+twin(2),:);
                    end
                    lfpAvg{rr,zz} = cat(3,lfpAvg{rr,zz},lfp_temp);
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
    
    for rr = 1:3
        for zz = 1:2
           lfp_avg = nanmean(lfpAvg{rr,zz},3)*-1;
           
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
           csdData{rr,zz} = [0 mean(data(find(timestamps > -10 & timestamps < 30),:)) 0];
       
           if makePlots    
              subplot(3,3,3*(rr-1)+zz)
              taxis = linspace(-twin(1),twin(2),size(data,1));
              cmax = max(max(data));      
              contourf(taxis,1:size(data,2),data',40,'LineColor','none');hold on;
              set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel');  
              colormap jet; try caxis([-cmax cmax]); end
              hold on
              for kk = 1:size(lfp_avg,2)
                  plot(taxis,(lfp_avg(:,kk)/1000)+kk-1,'k')
              end
               title(num2str(size(lfpAvg{rr,zz},3)))
              
              subplot(3,3,3*(rr-1)+3)
              if zz == 1
                  plot(csdData{rr,zz},1:length(channel_order),'color',[85/243 85/243 85/243],'LineWidth',1.5);     
              else
                  plot(csdData{rr,zz},1:length(channel_order),'color',col(rr,:),'LineWidth',1.5);     
              end
              hold on       
           end
        end
    end

    if saveMat 
%         save([expPath '\Summ\' 'csdRippleData_shank2.mat'], 'csdData');
%         saveas(figure(1),strcat(expPath,'\Summ\CSDRipple_shank2.png'),'png');
%         saveas(figure(1),strcat(expPath,'\Summ\CSDRipple_shank2.fig'),'fig');
%         saveas(figure(1),strcat(expPath,'\Summ\CSDRipple_shank2.eps'),'epsc');
        save([expPath '\Summ\' 'csdRippleData.mat'], 'csdData');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRipple.png'),'png');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRipple.fig'),'fig');
        saveas(figure(1),strcat(expPath,'\Summ\CSDRipple.eps'),'epsc');        
    end
end
close all
end