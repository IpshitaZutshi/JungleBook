function UpdateSessPhaseLocking(varargin)

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

if exist(strcat('Summ\PhaseLocking.mat'),'file') && ~force 
    disp('PhaseLocking already computed! Loading file.');
    load(strcat('Summ\PhaseLocking.mat'));
%else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                phaseData.deepSuperficial{rr,cc}{zz} = [];
                phaseData.region{rr,cc}{zz} = [];
                phaseData.putativeCellType{rr,cc}{zz} = [];                
            end
        end
    end    

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        load([sessionInfo.FileName '.deepSuperficialfromRipple.channelinfo.mat']);
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


        efields = fieldnames(sessionPulses);    
        
        deepSup = [];
        
        for unit = 1:length(cell_metrics.maxWaveformCh1)
            if strcmp(deepSuperficialfromRipple.channelClass{cell_metrics.maxWaveformCh1(unit)},'Superficial')==1
                deepSup(1,unit) = 1;
            elseif strcmp(deepSuperficialfromRipple.channelClass{cell_metrics.maxWaveformCh1(unit)},'Deep')==1
                deepSup(1,unit) = 2;
            else
                deepSup(1,unit) = 0;
            end
        end
        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            
            for zz = 1:6
                phaseData.deepSuperficial{region,target}{zz} = catpad(3,phaseData.deepSuperficial{region,target}{zz},...
                    deepSup);
                phaseData.region{region,target}{zz} = catpad(3,phaseData.region{region,target}{zz},cell_metrics.brainRegion);
                phaseData.putativeCellType{region,target}{zz} = catpad(3,phaseData.putativeCellType{region,target}{zz},cell_metrics.putativeCellType);
            end

        end

    end

    if saveMat
        save([expPath '\Summ\' 'PhaseLocking.mat'], 'phaseData');
    end
%end
end