function SessBehaviorAssembliesEntireTrack(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'lowThresh',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'defineSep',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
lowThresh = p.Results.lowThresh;
makePlot = p.Results.makePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;
defineSep = p.Results.defineSep;


if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist(strcat('Summ\AssembliesEntireTrack.mat'),'file') && ~force 
    disp('Assemblies on the entire track already computed! Loading file.');
    load(strcat('Summ\AssembliesEntireTrack.mat'));
else
    for rr = 1:3
        for cc = 1:2    
            Assemblies{rr,cc} = [];
        end
    end
    for ii = 1:size(allSess,1)
        
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        file = dir(('*.Behavior.mat'));
        load(file(1).name);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        if ~lowThresh
            file = dir(('*.ripples.events.mat'));        
            load(file.name);    
            ripplesLowThresh = ripples;
        else
            file = dir(('*.ripplesLowThresh.events.mat'));
            load(file.name);  
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
           
           stim = sessionPulses.(efields{jj}).stim;
           events = [rewardTS(1:length(startDelay)) startDelay(1:end)];
            
           idx = behavior.masks.recording==jj;
           behavTS = behavior.timestamps(idx);

           % Identify ripples within the session
           ripIdxBase = InIntervals(ripplesLowThresh.peaks,events(stim(1:length(startDelay))==0,:));
           ripIdxStim = InIntervals(ripplesLowThresh.peaks,events(stim(1:length(startDelay))==1,:));
            
           % Match ripple numbers between baseline and stim 
           if sum(ripIdxBase)>= sum(ripIdxStim)
              idx = find(ripIdxBase);
              yBase = datasample(idx,sum(ripIdxStim), 'Replace',false);             
              yStim = find(ripIdxStim);
           elseif sum(ripIdxBase)<= sum(ripIdxStim)
              idx = find(ripIdxStim);
              yStim = datasample(idx,sum(ripIdxBase), 'Replace',false);   
              yBase = find(ripIdxBase);
           end
           
           if length(yStim)<5
               continue
           end           
           
           % Take subsets 
           arr = [];
           [yBaseHalf1, idx1]  = datasample(yBase,round(length(yBase)/2), 'Replace',false);       
           arr(1:length(yBase)) = 1;
           arr(idx1) = 0;
           yBaseHalf2  = yBase(logical(arr));
           
           arr = [];
           [yStimHalf1, idx1]  = datasample(yStim,round(length(yStim)/2), 'Replace',false);       
           arr(1:length(yStim)) = 1;
           arr(idx1) = 0;
           yStimHalf2  = yStim(logical(arr));  
           
           spikeMat = [];

           for unit = 1:length(spikes.times)
               if (strcmp(cell_metrics.brainRegion(unit),'CA1')~=1) || (strcmp(cell_metrics.putativeCellType(unit),'Pyramidal Cell')~=1)
                   continue
               else
                   spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .025, 'win',[behavTS(1) behavTS(end)]);
                   spikeMat = [spikeMat spkData.data];
               end
           end
           
           timeMat = spkData.timestamps;
           
           SpikeCount = zscore(spikeMat)';
           AssemblyTemplates = assembly_patterns(SpikeCount);
           
           %Make vectors for ripple times    
           ripLog = [];
           for rr = 1:6
                ripLog(rr,1:length(timeMat)) = 0;
           end
        
           peakTS = ripplesLowThresh.peaks(yBase);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.25) & timeMat<=(peakTS(tt)+0.25);
               if sum(idx)>0
                   ripLog(1,idx) = 1;
               end                               
           end
           
           peakTS = ripplesLowThresh.peaks(yBaseHalf1);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
               if sum(idx)>0
                   ripLog(2,idx) = 1;
               end                               
           end
           
           peakTS = ripplesLowThresh.peaks(yBaseHalf2);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
               if sum(idx)>0
                   ripLog(3,idx) = 1;
               end                               
           end           
           
           peakTS = ripplesLowThresh.peaks(yStim);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
               if sum(idx)>0
                   ripLog(4,idx) = 1;
               end                               
           end           
           
           peakTS = ripplesLowThresh.peaks(yStimHalf1);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
               if sum(idx)>0
                   ripLog(5,idx) = 1;
               end                               
           end
           
           peakTS = ripplesLowThresh.peaks(yStimHalf2);          
           for tt = 1:size(peakTS,1)
               %Find bins within 500 ms of the ripple peak
               idx = timeMat>=(peakTS(tt)-0.5) & timeMat<=(peakTS(tt)+0.5);
               if sum(idx)>0
                   ripLog(6,idx) = 1;
               end                               
           end
           
           % find cell assemblies
           for aa = 1:size(AssemblyTemplates,2)
                assembliesVector = AssemblyTemplates(:,aa);
                % flip is negative weights:
                assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                    /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
                Vectors(:,aa) = assembliesVector;
%                     assembliesID{length(assembliesID)+1} = find(assembliesVector> 2*std(assembliesVector)) + cellCount;
%                     sessionID(length(sessionID)+1) = ii;
           end    
           
           if exist('Vectors','var')==1
               Activities = assembly_activity(Vectors,SpikeCount);

               % Now, expression of each assembly during Baseline ripples
               assembly_base = mean(Activities(:,ripLog(1,:)==1),2);
               assembly_basehalf1 = mean(Activities(:,ripLog(2,:)==1),2);
               assembly_basehalf2 = mean(Activities(:,ripLog(3,:)==1),2);
               
               % Now, expression of each assembly during stim ripples
               assembly_stim = mean(Activities(:,ripLog(4,:)==1),2);
               assembly_stimhalf1 = mean(Activities(:,ripLog(5,:)==1),2);
               assembly_stimhalf2 = mean(Activities(:,ripLog(6,:)==1),2);

               Assemblies{region,target} =  [Assemblies{region,target}; [assembly_base assembly_stim assembly_basehalf1 ...
                   assembly_basehalf2 assembly_stimhalf1 assembly_stimhalf2]];               
           end
           
           clear Vectors
        end

    end

    if saveMat
        if ~lowThresh
            save([expPath '\Summ\' 'AssembliesEntireTrack.mat'], 'Assemblies','-v7.3');
        else
            save([expPath '\Summ\' 'AssembliesEntireTrackLowThresh.mat'], 'Assemblies','-v7.3');
        end
    end
end

