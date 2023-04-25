function SessBehaviorAssembliesSeparateCA3SideArm(varargin)

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

    Assemblies.expressionBase = [];
    Assemblies.expressionmEC = [];  
    Assemblies.expressionCA3 = [];
    Assemblies.expressionBoth = [];  
    
    Assemblies.numberBase = [];
    Assemblies.numbermEC = [];               
    Assemblies.numberCA3 = [];
    Assemblies.numberBoth = [];               
    
    Assemblies.numberBaseA = [];
    Assemblies.numbermECA = [];               
    Assemblies.numberCA3A = [];
    Assemblies.numberBothA = [];   
    
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
        stimMat = []; 
        spikeMat = [];
        
        %Skip stem targeting sessions
        if sessionPulses.(efields{1}).target==1
            continue
        end
        
        for jj = 1:length(efields)
                        
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            
            rewardTS = sessionArmChoice.(efields{jj}).timestamps; 
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';  
            stim = sessionPulses.(efields{jj}).stim;
            events = [rewardTS(1:length(startDelay)) startDelay(1:end)];
            
  
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
               if jj == 1
                   stimMat = [stimMat;ones(size(spkMat,1),1)*stim(trials)];
               elseif jj ==2
                   stimMat = [stimMat;ones(size(spkMat,1),1)*(stim(trials)+2)];
               end
            end  
        end
            
       SpikeCount = zscore(spikeMat)';

       SpikeCountBase = zscore(spikeMat(stimMat==0,:))';
       AssemblyTemplatesBase = assembly_patterns(SpikeCountBase);
       SpikeCountmEC = zscore(spikeMat(stimMat==1,:))';
       AssemblyTemplatesmEC = assembly_patterns(SpikeCountmEC);
       SpikeCountCA3 = zscore(spikeMat(stimMat==2,:))';
       AssemblyTemplatesCA3 = assembly_patterns(SpikeCountCA3);
       SpikeCountBoth = zscore(spikeMat(stimMat==3,:))';
       AssemblyTemplatesBoth = assembly_patterns(SpikeCountBoth);
       
       % find cell assemblies
       for aa = 1:size(AssemblyTemplatesBase,2)
            assembliesVector = AssemblyTemplatesBase(:,aa);
            % flip is negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            VectorsBase(:,aa) = assembliesVector;
            Assemblies.numberBase  = [Assemblies.numberBase sum(assembliesVector> 2*std(assembliesVector))];
       end   
       Assemblies.numberBaseA = [Assemblies.numberBaseA size(AssemblyTemplatesBase,2)];

       % find cell assemblies
       for aa = 1:size(AssemblyTemplatesmEC,2)
            assembliesVector = AssemblyTemplatesmEC(:,aa);
            % flip is negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            VectorsmEC(:,aa) = assembliesVector;
            Assemblies.numbermEC  = [Assemblies.numbermEC sum(assembliesVector> 2*std(assembliesVector))];
       end  
       Assemblies.numbermECA = [Assemblies.numbermECA size(AssemblyTemplatesmEC,2)];

       % find cell assemblies CA3
       for aa = 1:size(AssemblyTemplatesCA3,2)
            assembliesVector = AssemblyTemplatesCA3(:,aa);
            % flip is negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            VectorsCA3(:,aa) = assembliesVector;
            Assemblies.numberCA3  = [Assemblies.numberCA3 sum(assembliesVector> 2*std(assembliesVector))];
       end  
       Assemblies.numberCA3A = [Assemblies.numberCA3A size(AssemblyTemplatesCA3,2)];
       
       % find cell assemblies Both
       for aa = 1:size(AssemblyTemplatesBoth,2)
            assembliesVector = AssemblyTemplatesBoth(:,aa);
            % flip is negative weights:
            assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
                /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
            VectorsBoth(:,aa) = assembliesVector;
            Assemblies.numberBoth = [Assemblies.numberBoth sum(assembliesVector> 2*std(assembliesVector))];
       end  
       Assemblies.numberBothA = [Assemblies.numberBothA size(AssemblyTemplatesBoth,2)];
       
       if exist('VectorsBase','var')==1 
           ActivitiesBase = assembly_activity(VectorsBase,SpikeCount);

           avgExp = [];              
           assembly = ActivitiesBase(:,stimMat==0);
           avgExp(:,1) = mean(assembly,2);                              
           assembly = ActivitiesBase(:,stimMat==1);
           avgExp(:,2) = mean(assembly,2);    
           assembly = ActivitiesBase(:,stimMat==2);
           avgExp(:,3) = mean(assembly,2);    
           assembly = ActivitiesBase(:,stimMat==3);
           avgExp(:,4) = mean(assembly,2);    
           
           % define stability
           baseidx = find(stimMat==0);
           assembly_first = ActivitiesBase(:,baseidx(1:floor(length(baseidx)/2)));
           assembly_second = ActivitiesBase(:,baseidx(floor(length(baseidx)/2)+1:end)); 
           avgExp(:,5) = mean(assembly_first,2); avgExp(:,6) = mean(assembly_second,2);     

           Assemblies.expressionBase =  [Assemblies.expressionBase; avgExp];
       end

       if exist('VectorsmEC','var')==1

           ActivitiesStim = assembly_activity(VectorsmEC,SpikeCount);
           avgExp = [];

           assembly = ActivitiesStim(:,stimMat==0);
           avgExp(:,1) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==1);
           avgExp(:,2) = mean(assembly,2);   
           assembly = ActivitiesStim(:,stimMat==2);
           avgExp(:,3) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==3);
           avgExp(:,4) = mean(assembly,2);   
           
           % define stability
           baseidx = find(stimMat==1);
           assembly_first = ActivitiesStim(:,baseidx(1:floor(length(baseidx)/2)));
           assembly_second = ActivitiesStim(:,baseidx(floor(length(baseidx)/2)+1:end)); 
           avgExp(:,5) = mean(assembly_first,2); avgExp(:,6) = mean(assembly_second,2);     

           Assemblies.expressionmEC =  [Assemblies.expressionmEC; avgExp];
       end

       if exist('VectorsCA3','var')==1

           ActivitiesStim = assembly_activity(VectorsCA3,SpikeCount);
           avgExp = [];

           assembly = ActivitiesStim(:,stimMat==0);
           avgExp(:,1) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==1);
           avgExp(:,2) = mean(assembly,2);    
           assembly = ActivitiesStim(:,stimMat==2);
           avgExp(:,3) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==3);
           avgExp(:,4) = mean(assembly,2);    
           
           % define stability
           baseidx = find(stimMat==2);
           assembly_first = ActivitiesStim(:,baseidx(1:floor(length(baseidx)/2)));
           assembly_second = ActivitiesStim(:,baseidx(floor(length(baseidx)/2)+1:end)); 
           avgExp(:,5) = mean(assembly_first,2); avgExp(:,6) = mean(assembly_second,2);     

           Assemblies.expressionCA3 =  [Assemblies.expressionCA3; avgExp];
       end

       if exist('VectorsBoth','var')==1

           ActivitiesStim = assembly_activity(VectorsBoth,SpikeCount);
           avgExp = [];

           assembly = ActivitiesStim(:,stimMat==0);
           avgExp(:,1) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==1);
           avgExp(:,2) = mean(assembly,2);     
           assembly = ActivitiesStim(:,stimMat==2);
           avgExp(:,3) = mean(assembly,2);                              
           assembly = ActivitiesStim(:,stimMat==3);
           avgExp(:,4) = mean(assembly,2);     
           
           % define stability
           baseidx = find(stimMat==3);
           assembly_first = ActivitiesStim(:,baseidx(1:floor(length(baseidx)/2)));
           assembly_second = ActivitiesStim(:,baseidx(floor(length(baseidx)/2)+1:end)); 
           avgExp(:,5) = mean(assembly_first,2); avgExp(:,6) = mean(assembly_second,2);     

           Assemblies.expressionBoth =  [Assemblies.expressionBoth; avgExp];
       end

       clear VectorsBase VectorsmEC VectorsCA3 VectorsBoth

    end

    if saveMat
        save([expPath '\Summ\' 'AssembliesSeparateSideArm.mat'], 'Assemblies','-v7.3');
    end
end

