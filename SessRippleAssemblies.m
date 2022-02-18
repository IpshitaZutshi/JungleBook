function SessRippleAssemblies(varargin)

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

if exist(strcat('Summ\RippleAssemblies.mat'),'file') && ~force 
    disp('Ripple Assemblies already computed! Loading file.');
    load(strcat('Summ\RippleAssemblies.mat'));
else
    for rr = 1:3
        for cc = 1:2    
            RippleAssemblies{rr,cc} = [];
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
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            stim = sessionPulses.(efields{jj}).stim;
            events = [endDelay(1:length(rewardTS)-1) rewardTS(2:end)];
            
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
               if makePlot
                   figure
                   set(gcf,'Position',[10 10 525 936])
                   for kk = 1:size(Activities,1)
                       subplot(size(Activities,1),1,kk)
                       plot(Activities(kk,:))
                       hold on
                       plot(stimMat*100)
                   end
                   saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\AssemblyActivation.fig'))
                   saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\AssemblyActivation.png'))
               end
               assembly_base = mean(Activities(:,stimMat==0),2);
               assembly_stim = mean(Activities(:,stimMat==1),2);
               baseidx = find(stimMat==0);
               assembly_first = mean(Activities(:,baseidx(1:floor(length(baseidx)/2))),2);
               assembly_second = mean(Activities(:,baseidx(floor(length(baseidx)/2)+1:end)),2);               
               %assemblydiff = (assembly_base-assembly_stim)./(assembly_base+assembly_stim);

               Assemblies{region,target} =  [Assemblies{region,target}; [assembly_base assembly_stim assembly_first assembly_second]]; 

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

