function SessUnitISI(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[6 14],@isnumeric);
addParameter(p,'bins',-3:0.05:1.5,@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;
bins = p.Results.bins;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
pyrChPlus = 0;

if exist(strcat('Summ\ISIData.mat'),'file') && ~force 
    disp('ISI already computed! Loading file.');
    load(strcat('Summ\ISIData.mat'));
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                ISIData.ISI{rr,cc}{zz} = [];
                ISIData.LFPfreq{rr,cc}{zz} = [];
                ISIData.region{rr,cc}{zz} = [];
                ISIData.putativeCellType{rr,cc}{zz} = [];
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
        file = dir(('*.hippocampalLayers.channelInfo.mat'));
        load(file.name);    
        orCh = hippocampalLayers.oriens; 
  
        if ~isfield(session.channelTags,'RippleNoise')
            lfp = bz_GetLFP(orCh,'noPrompts', true);
        else
            refChannel = session.channelTags.RippleNoise.channels-1;
            lfp = bz_GetLFP([orCh refChannel],'noPrompts', true);
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
            load(strcat(efields{jj},'\',allSess(ii).name,'.placeFields.cellinfo.mat'))
            
            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

            for zz = [1 2 4 5]
                %Extract relevant intervals for ISI 
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
                
                [ISI,lfpfreq] = getISI(spikes,lfp,events',bins,passband);
                
                ISIData.ISI{region,target}{zz} = [ISIData.ISI{region,target}{zz};ISI];    
                ISIData.LFPfreq{region,target}{zz} = [ISIData.LFPfreq{region,target}{zz};lfpfreq]; 
                ISIData.region{region,target}{zz} = catpad(1,ISIData.region{region,target}{zz},cell_metrics.brainRegion');
                ISIData.putativeCellType{region,target}{zz} = catpad(1,ISIData.putativeCellType{region,target}{zz},cell_metrics.putativeCellType');

            end   
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ISIData.mat'], 'ISIData');
    end
end
end

function [ISI,lfpfreq] = getISI(spikes,lfp,intervals,bins,passband)

    ISI = [];
    lfpfreq = [];
    %% Get peak frequency for every time point in LFP
    if ~isempty(lfp)
        [wave,f_lfp,t_lfp]=getWavelet(double(lfp.data(:,1)),lfp.samplingRate,passband(1),passband(2),40,0);
        [~,mIdx]=max(wave);%get index max power for each timepoint
        freq = f_lfp(mIdx);
        t_lfp = sort(t_lfp);
    end

    for a = 1:length(spikes.times)
        bools = InIntervals(spikes.times{a},intervals);
        s =spikes.times{a}(bools);        
        if length(s)<5
           ISI(a,1:length(bins)-1) = nan;  
           if ~isempty(lfp)
               lfpfreq(a,1) = nan;  
           end
        else   
            ISI(a,:) = histcounts(log10(diff(s)),bins);
            if ~isempty(lfp)           
                %Find lfp power at the spiketimes
                s = sort(s);
                [~, ~, bins1] = histcounts(t_lfp,s);
                [~,y] = unique(bins1);
                closestIndex = y(2:end);
                lfpFreq = freq(closestIndex);
                lfpfreq(a,1) = nanmean(lfpFreq);
            end       
        end
    end
end

