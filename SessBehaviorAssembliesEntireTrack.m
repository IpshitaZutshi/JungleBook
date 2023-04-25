function SessBehaviorAssembliesEntireTrack(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'defineSep',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
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

if exist(strcat('Summ\AssembliesSideArm.mat'),'file') && ~force 
    disp('Assembliesalready computed! Loading file.');
    load(strcat('Summ\AssembliesSideArm.mat'));
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
        file = dir(('*.ripples.events.mat'));
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
            
           idx = behavior.masks.recording==jj;
           behavTS = behavior.timestamps(idx);

           % Identify ripples within the session
           ripIdxBase = InIntervals(ripples.peaks,events(stim(1:length(startDelay))==0,:));
           ripIdxStim = InIntervals(ripples.peaks,events(stim(1:length(startDelay))==1,:));
            
           spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .025, 'win',[behavTS(1) behavTS(end)]);
           spikeMat = spkData.data;
           SpikeCount = zscore(spikeMat)';
           AssemblyTemplates = assembly_patterns(SpikeCount);
           
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

               assembly_base = mean(Activities(:,stimMat==0),2);
               assembly_stim = mean(Activities(:,stimMat==1),2);
               baseidx = find(stimMat==0);
               assembly_first_base = mean(Activities(:,baseidx(1:floor(length(baseidx)/2))),2);
               assembly_second_base = mean(Activities(:,baseidx(floor(length(baseidx)/2)+1:end)),2);   
               baseidx = find(stimMat==1);
               assembly_first_stim = mean(Activities(:,baseidx(1:floor(length(baseidx)/2))),2);
               assembly_second_stim = mean(Activities(:,baseidx(floor(length(baseidx)/2)+1:end)),2);                          
               %assemblydiff = (assembly_base-assembly_stim)./(assembly_base+assembly_stim);

               Assemblies{region,target} =  [Assemblies{region,target}; [assembly_base assembly_stim assembly_first_base assembly_second_base ...
                   assembly_first_stim assembly_second_stim]]; 

               if makePlot
                   close all
               end
           end
           clear Vectors
        end

    end

    if saveMat
        save([expPath '\Summ\' 'Assemblies.mat'], 'Assemblies','-v7.3');
    end
end

