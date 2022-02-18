function SessBehaviorThetaCompression(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[7 12],@isnumeric);
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
plotfig = p.Results.plotfig;
passband = p.Results.passband;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
pyrChPlus = 0;

if exist(strcat('Summ\ThetaCompression.mat'),'file') && ~force 
    disp('Theta compression already computed! Loading file.');
    load(strcat('Summ\ThetaCompression.mat'));
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:2
                thetaComp.pairIDs{rr,cc}{zz} = [];
                thetaComp.placefield_difference{rr,cc}{zz} = [];
                thetaComp.placefield_size{rr,cc}{zz} = [];
                thetaComp.placefield_center{rr,cc}{zz} = [];
                thetaComp.ccgs_place_offset{rr,cc}{zz} = [];
                thetaComp.ccgs_time_offset{rr,cc}{zz} = [] ;
                thetaComp.ccgs_phase_offset{rr,cc}{zz} = [];
                thetaComp.ccgs_time{rr,cc}{zz} = [] ;
                thetaComp.ccgs_phase{rr,cc}{zz} = []; 
                thetaComp.compression{rr,cc}{zz} = [] ;
                thetaComp.ccgs_time_peaks{rr,cc}{zz} = [] ;
                thetaComp.ccgs_phase_peaks{rr,cc}{zz} = [];                 
                thetaComp.time{rr,cc}{zz} = [] ;
                thetaComp.phase{rr,cc}{zz} = [];                 
            end
            thetaComp.region{rr,cc} = [];
            thetaComp.putativeCellType{rr,cc} = [];          
            thetaComp.sessCount{rr,cc} = [];
            sessCount{rr,cc} = 1;
        end
    end    

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        file = dir(('*.hippocampalLayers.channelinfo.mat'));
        load(file.name);    
        pyrCh = hippocampalLayers.oriens; 
   
        if ~isfield(session.channelTags,'RippleNoise')
            lfp = bz_GetLFP(pyrCh,'noPrompts', true);
        else
            refChannel = session.channelTags.RippleNoise.channels-1;
            lfp = bz_GetLFP([pyrCh refChannel],'noPrompts', true);
            lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        end
            
        if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
             spikes = bz_LoadPhy('noPrompts',true);
        else
            continue;
        end

        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            arm = sessionArmChoice.(efields{jj}).visitedArm(2:end);
            choice = sessionArmChoice.(efields{jj}).choice(2:end);
            
            cd(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}));
            load(strcat(allSess(ii).name,'.placeFields.cellinfo.mat'))                  
            
            for zz = 1:2%4
                
               fieldMax = []; fieldSize =[];                    

               %Extract relevant intervals for CCGs
               fprintf(' ** Examining zone%3.i\n',zz);
                
               switch zz %% for some reason armchoice arms are FUTURE arm, not current - correct
                    case 1  %First, no stim trials   
%                         startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0 & arm' == 0);
%                         endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0 & arm' == 0)+1);
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
                        events = [startTS'; endTS'];
                    case 2  %stim
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
                        events = [startTS';endTS'];                       
               end
                
               %Collect place field info
               for kk = 1:length(placeFieldStats.mapStats)
                   if ~isnan(placeFieldStats.mapStats{kk,1}{1,zz}.x(1)) || ~isnan(placeFieldStats.mapStats{kk,1}{1,zz+2}.x(1))
                       [~,idxStem1] = max(placeFieldStats.mapStats{kk,1}{1,zz}.x);
                       [~,idxStem2] = max(placeFieldStats.mapStats{kk,1}{1,zz+2}.x);
                       x1 = placeFieldStats.mapStats{kk,1}{1,zz}.x(idxStem1);
                       x2 = placeFieldStats.mapStats{kk,1}{1,zz+2}.x(idxStem2);  
                       
                       sz1 = placeFieldStats.mapStats{kk,1}{1,zz}.size(idxStem1);
                       sz2 = placeFieldStats.mapStats{kk,1}{1,zz+2}.size(idxStem2);

                       fieldMax = [fieldMax mean([x1 x2])];                
                       fieldSize = [fieldSize mean([sz1 sz2])];
                   else                        
                       fieldMax = [fieldMax nan];                     
                       fieldSize = [fieldSize nan];
                   end
               end                
                
               tComp = bz_thetaCompression(spikes,lfp,'passband',passband,'intervals',events','fieldLoc',fieldMax,'fieldSize',fieldSize,'zoneNum',zz);
               
               thetaComp.pairIDs{region,target}{zz} = [thetaComp.pairIDs{region,target}{zz};tComp.pairIDs];
               thetaComp.placefield_difference{region,target}{zz}  = [thetaComp.placefield_difference{region,target}{zz};tComp.placefield_difference];
               thetaComp.placefield_size{region,target}{zz}  = [thetaComp.pairIDs{region,target}{zz};tComp.placefield_size];
               thetaComp.placefield_center{region,target}{zz} = [thetaComp.placefield_center{region,target}{zz};tComp.placefield_center];
               thetaComp.ccgs_place_offset{region,target}{zz} = [thetaComp.ccgs_place_offset{region,target}{zz};tComp.ccgs_place_offset];
               thetaComp.ccgs_time_offset{region,target}{zz} = [thetaComp.ccgs_time_offset{region,target}{zz};tComp.ccgs_time_offset];
               thetaComp.ccgs_phase_offset{region,target}{zz} = [thetaComp.ccgs_phase_offset{region,target}{zz};tComp.ccgs_phase_offset]; 
               thetaComp.ccgs_time{region,target}{zz} = [thetaComp.ccgs_time{region,target}{zz};tComp.ccgs_time];
               thetaComp.ccgs_phase{region,target}{zz} = [thetaComp.ccgs_phase{region,target}{zz};tComp.ccgs_phase]; 
               thetaComp.compression{region,target}{zz} = [thetaComp.compression{region,target}{zz};tComp.compression];
               thetaComp.ccgs_time_peaks{region,target}{zz} = [thetaComp.ccgs_time_peaks{region,target}{zz};tComp.ccgs_time_peaks] ;
               thetaComp.ccgs_phase_peaks{region,target}{zz} = [thetaComp.ccgs_phase_peaks{region,target}{zz};tComp.ccgs_phase_peaks];                                      
               thetaComp.time{rr,cc}{zz} = tComp.t_time;
               thetaComp.phase{rr,cc}{zz} = tComp.t_phase;    
            end
                          
            k = 1;
            % Assign region  & pyr INT to each cell pair
            %RegDef = 1 if 'CA1-CA1 pair', 2 if CA1-DG/CA3 pair, 3 if
            %DG/CA3-DG-CA3 pair
            
            % putCellType - 1 if pyr-pyr, 2 if pyr-int, 3 if int-int
            for i = 1:(length(spikes.times)-1)
                for j = (i+1):length(spikes.times)
                    if strcmp(cell_metrics.brainRegion{i},'CA1')==1 && strcmp(cell_metrics.brainRegion{j},'CA1')==1 
                        regDef(k,1) = 1;
                    elseif strcmp(cell_metrics.brainRegion{i},'CA1')==1 || strcmp(cell_metrics.brainRegion{j},'CA1')==1 
                        regDef(k,1) = 2;
                    else
                        regDef(k,1) = 3;
                    end
                    
                    if strcmp(cell_metrics.putativeCellType{i},'Pyramidal Cell')==1 && strcmp(cell_metrics.putativeCellType{j},'Pyramidal Cell')==1 
                        putCellType(k,1) = 1;
                    elseif strcmp(cell_metrics.putativeCellType{i},'Pyramidal Cell')==1 || strcmp(cell_metrics.putativeCellType{j},'Pyramidal Cell')==1 
                        putCellType(k,1) = 2;
                    else
                        putCellType(k,1) = 3;
                    end                    
                    k = k+1;
                end
            end
            
            thetaComp.region{region,target} = [thetaComp.region{region,target}; regDef];
            thetaComp.putativeCellType{region,target} = [thetaComp.putativeCellType{region,target}; putCellType];
            thetaComp.sessCount{region,target} = [thetaComp.sessCount{region,target}; ones(length(tComp.placefield_center),1)*sessCount{region,target}];
            sessCount{region,target} = sessCount{region,target}+1;
            clear rewardTS startDelay events regDef putCellType
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ThetaCompression.mat'], 'thetaComp','-v7.3');
    end
end

if plotfig
    figure
    for ii = 2:3       
        for jj = 1
            for kk = 1:2
                idxtoKeep =  thetaComp.placefield_center{ii,jj}{kk}==1 & (thetaComp.placefield_difference{ii,jj}{kk}>-15 & thetaComp.placefield_difference{ii,jj}{kk}<15) & ...
                    (thetaComp.ccgs_place_offset{ii,jj}{kk}>-600 & thetaComp.ccgs_place_offset{ii,jj}{kk}<600);%...
                    %& thetaComp.region{ii,jj} == 1 & thetaComp.putativeCellType{ii,jj} == 1;
                
                subplot(4,6,2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (cm)')
                xlabel('Place field dist (ms)')
                    
                subplot(4,6,6+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (cm)')
                xlabel('Theta timescale (ms)')

                subplot(4,6,12+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (ms)')
                xlabel('Theta timescale (ms)')    
                
                subplot(4,6,18+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_phase_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_phase_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (ms)')
                xlabel('Theta phase (rad)')                                     
            end
        end
    end
end
end