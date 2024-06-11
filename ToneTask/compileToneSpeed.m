function Summary = compileToneSpeed

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  
  
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

Summary.r = [];
Summary.p = [];

Summary.r_int = [];
Summary.p_int = [];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.rateMapsTrial.cellinfo.mat']);
    trialMaps = load(file(1).name);    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);    
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);       
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);

    dtime = 0.1;%mean(diff(tracking.timestamps));
    win = [tracking.timestamps(1) tracking.timestamps(end)];
    spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'units','rate','win',win);

    % Interpolate
    timeS = tracking.timestamps;
    vel = interp1(timeS,tracking.position.vy,spkData.timestamps)';
    vel(vel>40) = nan;

    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information         
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end
 

        dataMat = [];
        dataMatTone = [];
        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end  

        spaceMapAvg = nanmean(dataMat,1);
        toneMap = nanmean(dataMatTone,1);        

        corrTone = []; 
        corrSpace = [];
        for pp = 1:6
           for jj = (pp+1):6      
               a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];
           end
        end     
        spatialCorr = nanmean(corrSpace);
        toneCorr = nanmean(corrTone);  

        %% Detect fields
        Field_Info = detectFields(toneMap,'maxFieldSize',40);
        if isempty(Field_Info)
            toneField = 0;
        else 
            toneField = 1;
        end        

        Field_Info = detectFields(spaceMapAvg);
        if isempty(Field_Info)
            spaceField = 0;
        else 
            spaceField = 1;
        end  

        % If its a lick cell, calculate its speed relationship
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) 
            
            spkRate  = spkData.data(:,kk);
            [r,p] = corrcoef(vel,spkRate,'Rows','pairwise');
            Summary.r = [Summary.r; r(1,2)];
            Summary.p = [Summary.p; p(1,2)];

        elseif cellType == 1 && (spaceField==1) && (spatialCorr > 0.1) 
            spkRate  = spkData.data(:,kk);
            [r,p] = corrcoef(vel,spkRate,'Rows','pairwise');
            Summary.r_int = [Summary.r_int; r(1,2)];
            Summary.p_int = [Summary.p_int; p(1,2)];
        end

    end
end

end

