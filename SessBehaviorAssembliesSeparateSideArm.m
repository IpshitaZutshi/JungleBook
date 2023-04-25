function SessBehaviorAssembliesSeparateSideArm(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',false,@islogical);
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

if exist(strcat('Summ\AssembliesSeparateSideArm.mat'),'file') && ~force 
    disp('Assembliesalready computed! Loading file.');
    load(strcat('Summ\AssembliesSeparateSideArm.mat'));
else
    for rr = 1:3
        for cc = 1:2
            Assemblies.expressionBase{rr,cc} = [];
            Assemblies.expressionStim{rr,cc} = [];            
            Assemblies.numberBase{rr,cc} = [];
            Assemblies.numberStim{rr,cc} = [];   
            Assemblies.numberBaseA{rr,cc} = [];
            Assemblies.numberStimA{rr,cc} = [];               
        end
    end
    Assemblies.assembliesIDStim = {};
    Assemblies.assembliesIDBase = {};
    Assemblies.sessionIDStim = [];
    Assemblies.sessionIDBase = [];
    Assemblies.regionIDStim = [];
    Assemblies.regionIDBase = [];
    Assemblies.targetIDStim = [];
    Assemblies.targetIDBase = [];        
    cellCountStim = 0;        
    cellCountBase = 0;     
    sessNum = 0;
    
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
            
            spikeMat = [];
            stimMat = [];   
            % Build matrix, concatenate trials
            for trials = 1:size(events,1)         
               spkMat = [];
               %spikes
               for unit = 1:length(spikes.times)
                   if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1
                       continue
                   else
                       spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .025, 'win',events(trials,:));
                       spkMat = [spkMat spkData.data];
                   end
               end
               spikeMat = [spikeMat;spkMat];
               stimMat = [stimMat;ones(size(spkMat,1),1)*stim(trials)];
            end  
            
           SpikeCount = zscore(spikeMat)';
               
           SpikeCountBase = zscore(spikeMat(stimMat==0,:))';
           AssemblyTemplatesBase = assembly_patterns(SpikeCountBase);
           
           SpikeCountStim = zscore(spikeMat(stimMat==1,:))';
           AssemblyTemplatesStim = assembly_patterns(SpikeCountStim);
           
           % find cell assemblies
           for aa = 1:size(AssemblyTemplatesBase,2)
                assembliesVector = AssemblyTemplatesBase(:,aa);
                % flip is negative weights:
                assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                    /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
                Assemblies.assembliesIDBase{length(Assemblies.assembliesIDBase)+1} = find(assembliesVector> 2*std(assembliesVector)) + cellCountBase;
                Assemblies.sessionIDBase(length(Assemblies.sessionIDBase)+1) = sessNum;
                Assemblies.regionIDBase(length(Assemblies.regionIDBase)+1) = region;                
                Assemblies.targetIDBase(length(Assemblies.targetIDBase)+1) = target; 
                VectorsBase(:,aa) = assembliesVector;
                Assemblies.numberBase{region,target}  = [Assemblies.numberBase{region,target} sum(assembliesVector> 2*std(assembliesVector))];
           end   
           Assemblies.numberBaseA{region,target}  = [Assemblies.numberBaseA{region,target} size(AssemblyTemplatesBase,2)];
               
           % find cell assemblies
           for aa = 1:size(AssemblyTemplatesStim,2)
                assembliesVector = AssemblyTemplatesStim(:,aa);
                % flip is negative weights:
                assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                    /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
                Assemblies.assembliesIDStim{length(Assemblies.assembliesIDStim)+1} = find(assembliesVector> 2*std(assembliesVector)) + cellCountStim;
                Assemblies.sessionIDStim(length(Assemblies.sessionIDStim)+1) = sessNum;
                Assemblies.regionIDStim(length(Assemblies.regionIDStim)+1) = region;                
                Assemblies.targetIDStim(length(Assemblies.targetIDStim)+1) = target; 
                VectorsStim(:,aa) = assembliesVector;
                Assemblies.numberStim{region,target}  = [Assemblies.numberStim{region,target} sum(assembliesVector> 2*std(assembliesVector))];
           end  
           Assemblies.numberStimA{region,target}  = [Assemblies.numberStimA{region,target} size(AssemblyTemplatesStim,2)];
           
           cellCountStim = cellCountStim + size(AssemblyTemplatesStim,1);
           cellCountBase = cellCountBase + size(AssemblyTemplatesBase,1);
           
           if exist('VectorsBase','var')==1 
               ActivitiesBase = assembly_activity(VectorsBase,SpikeCount);
               
               avgExp = [];              
               assembly = ActivitiesBase(:,stimMat==0);
               avgExp(:,1) = mean(assembly,2);                              
               assembly = ActivitiesBase(:,stimMat==1);
               avgExp(:,2) = mean(assembly,2);    
               
               % define stability
               baseidx = find(stimMat==0);
               assembly_first = ActivitiesBase(:,baseidx(1:floor(length(baseidx)/2)));
               assembly_second = ActivitiesBase(:,baseidx(floor(length(baseidx)/2)+1:end)); 
               avgExp(:,3) = mean(assembly_first,2); avgExp(:,4) = mean(assembly_second,2);     
              
               Assemblies.expressionBase{region,target} =  [Assemblies.expressionBase{region,target}; avgExp];
           end
               
           if exist('VectorsStim','var')==1
               
               ActivitiesStim = assembly_activity(VectorsStim,SpikeCount);
               avgExp = [];
               
               assembly = ActivitiesStim(:,stimMat==0);
               avgExp(:,1) = mean(assembly,2);                              
               assembly = ActivitiesStim(:,stimMat==1);
               avgExp(:,2) = mean(assembly,2);               
               % define stability
               baseidx = find(stimMat==1);
               assembly_first = ActivitiesStim(:,baseidx(1:floor(length(baseidx)/2)));
               assembly_second = ActivitiesStim(:,baseidx(floor(length(baseidx)/2)+1:end)); 
               avgExp(:,3) = mean(assembly_first,2); avgExp(:,4) = mean(assembly_second,2);     
              
               Assemblies.expressionStim{region,target} =  [Assemblies.expressionStim{region,target}; avgExp];
           end
           
           clear VectorsBase VectorsStim           
           sessNum = sessNum+1;
        end
    
    end

    if saveMat
        save([expPath '\Summ\' 'AssembliesSeparateSideArm.mat'], 'Assemblies','-v7.3');
    end
end

