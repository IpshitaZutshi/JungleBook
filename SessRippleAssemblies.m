function SessRippleAssemblies(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
makePlot = p.Results.makePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist(strcat('Summ\RippleAssemblies.mat'),'file') && ~force 
    disp('Ripple Assemblies already computed! Loading file.');
    load(strcat('Summ\RippleAssemblies.mat'));
else
    for rr = 1:3
        RippleAssemblies{rr} = [];
    end
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        
        %Load ripples
        load([sessionInfo.FileName '.ripples.events.mat']);
        load([sessionInfo.FileName '.pulses.events.mat']);
        
        %Load sleep states
        load([sessionInfo.FileName '.SleepState.states.mat']);
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);        
        
        load([sessionInfo.FileName '.SessionPulses.Events.mat']); 
        load([sessionInfo.FileName '.MergePoints.events.mat']); 
        
        %Find home cage intervals
        efields = fields(sessionPulses);
        ts_hc = [];
        for fn = 1:length(MergePoints.foldernames)
            flag = 1;
            for fx = 1:length(efields)
                if strcmp(efields{fx},MergePoints.foldernames{fn}) == 1 %If folder if a maze folder
                    flag = 0;
                end
            end
            if flag
                ts_hc = [ts_hc;MergePoints.timestamps(fn,:)];
            end
        end
        
        %Build spike matrix for sleep periods
        spikeMat = [];
        timeMat = [];   
        events = SleepState.ints.NREMstate;
        % Build matrix, concatenate sleep epochs
        for trials = 1:size(events,1)         
            spkMat = [];
            spkMat_ts = [];
            %check if interval lies within homecage
            status = InIntervals(events(trials,1),ts_hc);
            if status == 0
                continue
            end
            for unit = 1:length(spikes.times)
                if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1
                    continue
                else 
                   spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .025, 'win',events(trials,:));
                   spkMat = [spkMat spkData.data];                                        
                   spkMat_ts = [spkData.timestamps];
                end                               
            end
           spikeMat = [spikeMat;spkMat];
           timeMat = [timeMat;spkMat_ts];            
        end
                
        % Assign logicals for the timeMat 
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

            %Generate logicals for ripples in pre, prestim, post, poststim
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
        end
        
        % Assign logicals for the timeMat         
        timeLog(1,1:length(timeMat)) = 0;
        timeLog(2,1:length(timeMat)) = 0;
        timeLog(3,1:length(timeMat)) = 0;
        
        for rr = 1:3   
            for ripNum  = 1:4
                peakTS = ripples.peaks(ripple_logical{rr}(:,ripNum));
               % peakTS = ripples.timestamps(ripple_logical{rr}(:,ripNum),:);                
                for tt = 1:size(peakTS,1)
                    %Find bins within 500 ms of the ripple peak
                    idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
                    if sum(idx)>0
                        timeLog(rr,idx) = ripNum;
                    end                               
                end
            end
        end

       SpikeCount = zscore(spikeMat)';      
       
       if isempty(SpikeCount)
           continue;
       end              
       % Define assemblies for entire sleep periods
       AssemblyTemplates = assembly_patterns(SpikeCount);
       
       % find cell assemblies
       for aa = 1:size(AssemblyTemplates,2)
            assembliesVector = AssemblyTemplates(:,aa);
            % flip if negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            Vectors(:,aa) = assembliesVector;
%                     assembliesID{length(assembliesID)+1} = find(assembliesVector> 2*std(assembliesVector)) + cellCount;
%                     sessionID(length(sessionID)+1) = ii;
       end
       
      if exist('Vectors','var')==1
           Activities = assembly_activity(Vectors,SpikeCount);
           if makePlot
               figure
               set(gcf,'Position',[10 10 525 936])
               for kk = 1:size(Activities,1)
                   subplot(size(Activities,1),1,kk)
                   plot(Activities(kk,:))
                   hold on
                   plot(timeLog(1,:)*50,'r')
                   plot(timeLog(2,:)*50,'k')
                   %plot(timeLog(3,:)*50)
               end
               saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleAssemblyActivation_onlyRip.fig'))
               saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleAssemblyActivation_onlyRip.png'))
           end
           
           for rr = 1:3
               assemblyExp = [];
               for bin = 1:4
                   assemblyExp(:,bin) = mean(Activities(:,timeLog(rr,:)==bin),2);
               end
               RippleAssemblies{rr} = [RippleAssemblies{rr}; assemblyExp];
           end
           
           if makePlot
              close all
           end         
       end           
       clear Vectors timeLog ripple_logical
    end
    
    if saveMat
        save([expPath '\Summ\' 'RippleOnlyAssemblies.mat'], 'RippleAssemblies','-v7.3');
    end    
end
    
end

