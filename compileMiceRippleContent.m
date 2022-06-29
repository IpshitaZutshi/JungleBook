function compiledRipples = compileMiceRippleContent(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
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
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
         'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
         'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
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

for rr = 1:3
    for zz = 1:3
        compiledData{rr,zz} = [];
    end
end

if exist([parentDir 'Compiled\Ripples\compiledRippleContent' tag '.mat'],'file') && ~force 
    disp('Compiled ripple data already computed! Loading file.');
    load([parentDir 'Compiled\Ripples\compiledRippleContent' tag '.mat']);
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


            %% calculate ripple content
            for rr = 1:length(reg)              
                matPre{rr} = [];
                matPost{rr} = [];
                            
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
                ripple_prestim = [];  
                ripple_pre = [];                
                ripple_post = [];  
                ripple_poststim = [];                

                ripple_pre(1:length(ripples.peaks)) = 0;                
                ripple_post(1:length(ripples.peaks)) = 0;                
                ripple_prestim(1:length(ripples.peaks)) = 0;                
                ripple_poststim(1:length(ripples.peaks)) = 0;                              

                for pp = 1:length(ripples.peaks)
                    tempDiff = ripples.peaks(pp) - events;

                    if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                       [~,idxmin] =  min(abs(tempDiff));
                       if tempDiff(idxmin) > 0
                           ripple_post(pp) = 1;
                       elseif tempDiff(idxmin) < 0
                           ripple_pre(pp) = 1;
                       end
                    elseif min(abs(tempDiff)) >5 && min(abs(tempDiff)) <=10% If a ripple occurs within 5 seconds of a stimulus
                       [~,idxmin] =  min(abs(tempDiff));
                       if tempDiff(idxmin) > 0
                           ripple_poststim(pp) = 1;
                       elseif tempDiff(idxmin) < 0
                           ripple_prestim(pp) = 1;
                       end                   
                    else
                        continue
                    end
                end      
                
                ripple_logical{rr} = [logical(ripple_prestim)' logical(ripple_pre)'...
                    logical(ripple_post)' logical(ripple_poststim)'];     
                
                %% Now calculate the cells active for these specific ripples
                for rippLog = 1:size(ripple_logical{rr},2)
                    if sum(ripple_logical{rr}(:,rippLog))~=0
                        int = ripples.peaks(ripple_logical(:,rippLog));
                        
                        for jj = 1:size(spikes.UID,2)
                            % Get psth
                            if strcmp(cell_metrics.brainRegion{jj},'CA1')~=1 || strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell')~=1
                                continue
                            else
                                [status, interval, index] = InIntervals(spikes.times{jj},int);
                            end
                            a = zeros(1,length(int));
                            a(interval(status>0)) = 1;               
                            if rippLog == 1
                                matPre{rr} = [matPre{rr};a];
                            else
                                matPost{rr} = [matPost{rr};a];
                            end                          
                         end
                    end
                end                
                if ~isempty(matPre{rr})
                    rho_pre = corr(matPre{rr},'Type','Spearman');
                    corrVal_pre = tril(rho_pre,-1);
                    corrVal_pre = corrVal_pre(abs(corrVal_pre)>0);
                    compiledData{rr,1} = [compiledData{rr,1};corrVal_pre];
                end
                if ~isempty(matPost{rr})
                    rho_post = corr(matPost{rr},'Type','Spearman');
                    corrVal_post = tril(rho_post,-1);
                    corrVal_post = corrVal_post(abs(corrVal_post)>0);
                    compiledData{rr,2} = [compiledData{rr,2};corrVal_post];
                end
                
                if ~isempty(matPost{rr}) && ~isempty(matPre{rr})
                    rho_prepost = corr(matPre{rr},matPost{rr},'Type','Spearman');
                    corrVal_prepost = tril(rho_prepost,-1);
                    corrVal_prepost = corrVal_prepost(abs(corrVal_prepost)>0);
                    compiledData{rr,3} = [compiledData{rr,3};corrVal_prepost];                    
                end
            end
        end
    end
    
    if saveMat
        save([parentDir 'Compiled\Ripples\compiledRippleContent' tag '.mat'], 'compiledData');
    end
end

colMat = [85/243 85/243 85/243;...
        56/243 61/243 150/243;...
        224/243 163/243 46/243];  

    %% Plot the data
figure
set(gcf,'Renderer','painters')
set(gcf,'Position',[570 360 1320 622])

for rr = 1:3
    data = [];
    groups = [];
    data = [compiledData{rr,1};compiledData{rr,2};compiledData{rr,3}];
    groups = [ones(length(compiledData{rr,1}),1)*1;ones(length(compiledData{rr,2}),1)*2;ones(length(compiledData{rr,3}),1)*3];
    subplot(1,3,rr)
    stats{rr} = groupStats(data,groups,'inAxis',true,'color',colMat);    
end

saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple corr\RippleCorr',tag,'.png'))
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple corr\RippleCorr',tag,'.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Ripple corr\RippleCorr',tag,'.eps'),'epsc') 
save(strcat(parentDir,'Compiled\Ripples\Ripple corr\RippleCorr',tag,'.mat'),'stats')