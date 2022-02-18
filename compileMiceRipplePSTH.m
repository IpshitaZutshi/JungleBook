function compiledRipples = compileMiceRipplePSTH(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'binSize',0.1,@isnumerical);
addParameter(p,'normalized',true,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
saveMat = p.Results.saveMat;
force = p.Results.force;
binSize = p.Results.binSize;
normalized = p.Results.normalized;

tag = 'mECBilateral'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ18\DG_CA3','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ28\Final','IZ29\Final','IZ30\Final',...
        'IZ31\Final','IZ32\Final','IZ33\Final','IZ34\Final'};  % To add: IZ32, IZ33, IZ34
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline'};%,'IZ33\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

compiledData.mice = [];
compiledData.FR = [];
compiledData.TroughtoPeak = [];
compiledData.putativeClass = [];
compiledData.BurstIndex = [];

compiledData.psthtimes = [];
compiledData.psthstccg = [];
compiledData.analogCh = [];

compiledData.region = [];
compiledData.ripple = [];

if exist([parentDir 'Compiled\compiledRipplePSTH' tag '.mat'],'file') && ~force 
    disp('Compiled ripple data already computed! Loading file.');
    load([parentDir 'Compiled\compiledRipplePSTH' tag '.mat']);
else
    for m = 1:length(mice)

        cd(strcat(parentDir, mice{m}));
        allSess = dir('*_sess*');

        %% Start collecting data
        for ii = 1:size(allSess,1)

           cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
           [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

           %% Load spikes
           if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
                spikes = bz_LoadPhy('noPrompts',true);
           else
               continue
           end

           %% Load pulses
            if exist([sessionInfo.FileName '.pulses.events.mat'],'file') 
                load([sessionInfo.FileName '.pulses.events.mat']);
            else
                continue
            end

            %% Load cell metrics
            if exist([sessionInfo.FileName '.cell_metrics.cellinfo.mat'],'file') 
                load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
            else
                continue
            end
            
            %% Load ripples
            if exist([sessionInfo.FileName '.ripples.events.mat'],'file') 
                load([sessionInfo.FileName '.ripples.events.mat']);
            else
                continue;
            end

            %% calculate PSTH
            for rr = 1:length(reg)
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
                
                %% Now calculate the psth for these specific ripples
                for rippLog = 1:size(ripple_logical,2)
                    if sum(ripple_logical(:,rippLog))~=0
                        st = ripples.peaks(ripple_logical(:,rippLog));
                        for jj = 1:size(spikes.UID,2)
                            % Get psth
                            [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.01,'duration',0.5);
                            compiledData.psthtimes = [compiledData.psthtimes; t'];
                            compiledData.psthstccg = [compiledData.psthstccg; stccg(:,2,1)'];
                            compiledData.analogCh = [compiledData.analogCh; rr];
                            
                            %add cell metrics
                            compiledData.FR = [compiledData.FR; cell_metrics.firingRate(jj)]; 
                            compiledData.TroughtoPeak = [compiledData.TroughtoPeak; cell_metrics.troughToPeak(jj)];
                            compiledData.BurstIndex = [compiledData.BurstIndex; cell_metrics.burstIndex_Royer2012(jj)];

                            % Assign numerical tag to putative class
                            if strcmp(cell_metrics.putativeCellType{jj},'Narrow Interneuron') == 1
                                compiledData.putativeClass = [compiledData.putativeClass; 1];
                            elseif strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell') == 1
                                compiledData.putativeClass = [compiledData.putativeClass; 2];
                            elseif strcmp(cell_metrics.putativeCellType{jj},'Wide Interneuron') == 1
                                compiledData.putativeClass = [compiledData.putativeClass; 3];
                            else 
                                compiledData.putativeClass = [compiledData.putativeClass; 4];                       
                            end                    

                            % Assign numerical tag to region
                            if strcmp(cell_metrics.brainRegion{jj},'CA1') == 1
                                compiledData.region = [compiledData.region; 1];
                            elseif strcmp(cell_metrics.brainRegion{jj},'DG') == 1
                                compiledData.region = [compiledData.region; 2];
                            elseif strcmp(cell_metrics.brainRegion{jj},'CA3') == 1
                                compiledData.region = [compiledData.region; 2];
                            else 
                                compiledData.region = [compiledData.region; 3];                       
                            end

                            %Keep animal id data
                            compiledData.mice = [compiledData.mice; m];    
                            
                            %Ripple pre or post
                            compiledData.ripple = [compiledData.ripple; rippLog]; 
                        end
                    end
                end                

            end
        end
    end
    
    if saveMat
        save([parentDir 'Compiled\compiledRipplePSTH' tag '.mat'], 'compiledData');
    end
end

%% Plot the data
x = compiledData.psthtimes(1,:);

rippleType = {'pre','post'};
classifier = {'IN','Pyr'};
brainReg = {'CA1','CA3/DG'};
YlGnBu=cbrewer('div', 'BrBG', 11);

col = [85/243 85/243 85/243;224/243 163/243 46/243;8/243 133/243 161/243;56/243 61/243 150/243];
figure
set(gcf,'renderer','Painters')    
set(gcf,'Position',[35 50 1900 1000])    
for ii = 1:length(reg)
    for kk = 1:length(classifier)
        for rr = 1:2
             for jj = 1:length(rippleType)
                 
                idxRegion = compiledData.region == rr;                 
                idxclassifier = compiledData.putativeClass ==kk;                    
                idxRipple = compiledData.ripple == jj;
                idxStim = compiledData.analogCh == ii;               
                
                IdxCombined = idxRipple & idxStim & idxclassifier & idxRegion;
                psth = compiledData.psthstccg(IdxCombined,:);

                if normalized
                    psth = zscore(psth,[],2);%(:,:)./nanmean(psth(:,1:20),2);      
                end
                FR_ratio{jj} = nanmean(psth(:,24:28),2);
                [~, idxsort] = sort(FR_ratio{jj});                
            
                subplot(6,8,(4*(rr-1)+16*(ii-1)+2*(kk-1)+jj))
                imagesc(x,1:sum(IdxCombined),psth(idxsort,:))
                ylabel('Neuron ID')
                xlabel('Time (s)')
                title(strcat(reg(ii),', ',rippleType(jj),', ',classifier(kk),',',brainReg(rr)));
                if normalized
                    caxis([-2 5])
                else
                    caxis([0 5])    
                end
                %colormap(YlGnBu)
                
                subplot(6,8,(4*(rr-1)+16*(ii-1)+2*(kk-1)+1+8))
                hold on
                meanpsth = nanmean(psth,1);
                stdpsth = nanstd(psth,1)./sqrt(size(psth,1));               
                hold on
                if jj == 1
                    fill([x'; fliplr(x)'],[meanpsth'-stdpsth';flipud(meanpsth'+stdpsth')],col(1,:),'linestyle','none','FaceAlpha',0.5);                    
                    hi = line(x,meanpsth,'Color',col(1,:));
                else
                    fill([x'; fliplr(x)'],[meanpsth'-stdpsth';flipud(meanpsth'+stdpsth')],col(ii+1,:),'linestyle','none','FaceAlpha',0.5);                       
                    hi = line(x,meanpsth,'Color',col(ii+1,:));                    
                end
                if normalized
                    ylim([-0.5 2])
                end
%                 if kk == 1
%                     ylim([0 1])
%                 else
%                     ylim([0 1])
%                 end
             end
             
             if rr == 1
                data = [FR_ratio{1}' FR_ratio{2}'];
                stats{ii,kk} = groupStats(data,[ones(1,length(FR_ratio{1})) ones(1,length(FR_ratio{2}))*2],'doPlot',false); 
                [stats{ii,kk}.ranksum.p,~,stats{ii,kk}.ranksum.stats] = ranksum(FR_ratio{1},FR_ratio{2});
                clear FR_ratio data
             end
        end 
    end
end

if normalized
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\NormcompiledPSTH',tag,'.png'))
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\NormcompiledPSTH',tag,'.fig'),'fig')    
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\NormcompiledPSTH',tag,'.eps'),'epsc') 
    save(strcat(parentDir,'Compiled\Ripples\Ripple PSTH\NormcompiledPSTH',tag,'.mat'),'stats')
    clear stats
else
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\compiledPSTH',tag,'.png'))
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\compiledPSTH',tag,'.fig'),'fig')    
    saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple PSTH\compiledPSTH',tag,'.eps'),'epsc') 
    save(strcat(parentDir,'Compiled\Ripples\Ripple PSTH\compiledPSTH',tag,'.mat'),'stats')
    clear stats
end

end

