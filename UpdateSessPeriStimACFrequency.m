function UpdateSessPeriStimACFrequency(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'passband',[6 14],@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
pyrChPlus = 0;

if exist(strcat('Summ\ACData.mat'),'file') && ~force 
    disp('Autocorrelation already computed! Loading file.');
    load(strcat('Summ\ACData.mat'));
%else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:4
                ACData.placefield{rr,cc}{zz} = [];
                ACData.fieldMax{rr,cc}{zz} = [];                
                ACData.avgRate{rr,cc}{zz} = [];     
%                 ACData.spectrum{rr,cc}{zz} = [];
%                 ACData.freq{rr,cc}{zz} = [];
%                 ACData.LFPfreq{rr,cc}{zz} = [];
%                 ACData.Unitfreq{rr,cc}{zz} = [];
%                 ACData.Unitpower{rr,cc}{zz} = [];
%                 ACData.region{rr,cc}{zz} = [];
%                 ACData.putativeCellType{rr,cc}{zz} = [];
            end
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
        file = dir(('*.region.mat'));
        load(file.name);    
        pyrCh = region.CA1sp; 
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus);
            end
        end    
%         if ~isfield(session.channelTags,'RippleNoise')
%             lfp = bz_GetLFP(pyrChNew,'noPrompts', true);
%         else
%             refChannel = session.channelTags.RippleNoise.channels-1;
%             lfp = bz_GetLFP([pyrChNew refChannel],'noPrompts', true);
%             lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
%         end
%             
%         if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
%              spikes = bz_LoadPhy('noPrompts',true);
%         else
%             continue;
%         end

        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
            load(strcat(efields{jj},'\',allSess(ii).name,'.placeFields.cellinfo.mat'))
            
            for kk = 1:4
                placefield{kk} = [];fieldMax{kk} = [];avgRate{kk} = [];
            end
            for kk = 1:length(placeFieldStats.mapStats)                
                for zz = 1:size(placeFieldStats.mapStats{1},2)
                    if ~isnan(placeFieldStats.mapStats{kk,1}{1,zz}.x(1))
                       placefield{zz} = [placefield{zz}; 1];
                    else
                       placefield{zz} = [placefield{zz}; 0];
                    end
                    fieldMax{zz} = [fieldMax{zz} placeFieldStats.mapStats{kk,1}{1,zz}.x(1)];                
                    avgRate{zz} = [avgRate{zz} placeFieldStats.mapStats{kk,1}{1,zz}.mean(1)]; 
                end
            end
            for zz = 1:4
                ACData.placefield{region,target}{zz} = catpad(3,ACData.placefield{region,target}{zz},placefield{zz});
                ACData.fieldMax{region,target}{zz} = catpad(3,ACData.fieldMax{region,target}{zz},fieldMax{zz});              
                ACData.avgRate{region,target}{zz} = catpad(3,ACData.avgRate{region,target}{zz},avgRate{zz});   
            end
% 
%             rewardTS = sessionArmChoice.(efields{jj}).timestamps;
%             startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
%             endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

%             for zz = 1:6
%                 %Extract relevant intervals for cross-frequency coupling - 4 cross
%                 %modulograms
%                 fprintf(' ** Examining zone%3.i\n',zz);
%                 switch zz
%                     case 1  %First, no stim trials, return        
%                         startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
%                         endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
%                         events = [startTS'; endTS'];
%                     case 2  %No stim, stem
%                         startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
%                         endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
%                         events = [startTS';endTS'];
%                     case 3 %No stim, delay
%                         startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
%                         endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
%                         events = [startTS';endTS'];  
%                     case 4  % Stim, return
%                         startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
%                         endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
%                         events = [startTS';endTS'];                    
%                     case 5   % Stim, stem
%                         startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
%                         endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
%                         events = [startTS';endTS'];                      
%                     case 6    %stim, delay
%                         startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
%                         endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
%                         events = [startTS';endTS'];
%                 end
%                 
%                 if (zz == 3 || zz == 6) && sessionArmChoice.(efields{jj}).delay.dur < 1 
%                     ACDataRes.spectrum = nan;
%                     ACDataRes.freq = nan;
%                     ACDataRes.LFPfreq = nan;
%                     ACDataRes.Unitfreq = nan;
%                     ACDataRes.Unitpower = nan;
%                 else
%                     ACDataRes =  bz_ACFrequency(spikes,'lfp',lfp,'passband',passband,'intervals',events');
%                 end
%                 ACData.spectrum{region,target}{zz} = catpad(3,ACData.spectrum{region,target}{zz},ACDataRes.spectrum);    
%                 ACData.freq{region,target}{zz} = catpad(3,ACData.freq{region,target}{zz},ACDataRes.frequency); 
%                 ACData.LFPfreq{region,target}{zz} = catpad(3,ACData.LFPfreq{region,target}{zz},ACDataRes.lfpfreq); 
%                 ACData.Unitfreq{region,target}{zz} = catpad(3,ACData.Unitfreq{region,target}{zz},ACDataRes.unitfreq); 
%                 ACData.Unitpower{region,target}{zz} = catpad(3,ACData.Unitpower{region,target}{zz},ACDataRes.unitpower); 
%                 ACData.region{region,target}{zz} = catpad(3,ACData.region{region,target}{zz},cell_metrics.brainRegion);
%                 ACData.putativeCellType{region,target}{zz} = catpad(3,ACData.putativeCellType{region,target}{zz},cell_metrics.putativeCellType);
% 
%             end   
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ACData.mat'], 'ACData');
    end
%end
end