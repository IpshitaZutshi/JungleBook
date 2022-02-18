% Compile data across all sessions

function PF = compilePlaceFields(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');


if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for rr = 1:3
    for cc = 1:2
        for zz = 1:4 %4 fields, left trajectories, right trajectories, stim left and right trajectories
            PF.rateMap{rr,cc}{zz} = [];
            PF.peakRate{rr,cc}{zz} = [];
            PF.avgRate{rr,cc}{zz} = [];
        end
    end
end


%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
    file = dir(('*.SessionPulses.Events.mat'));
    load(file.name);
    
    efields = fieldnames(sessionPulses);    

    for jj = 1:length(efields)
        region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
        target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
        
        load(strcat(efields{jj},'\',allSess(ii).name,'.firingMapsAvg.cellinfo.mat'))
        for kk = 1:length(firingMaps.rateMaps)
            for zz = 1:size(firingMaps.rateMaps{1},2)
                PF.rateMap{region,target}{zz} = [PF.rateMap{region,target}{zz};firingMaps.rateMaps{kk,1}{1,zz}];
                PF.peakRate{region,target}{zz} = [PF.peakRate{region,target}{zz};max(firingMaps.rateMaps{kk,1}{1,zz})];
                PF.avgRate{region,target}{zz} = [PF.avgRate{region,target}{zz};nanmean(firingMaps.rateMaps{kk,1}{1,zz})];
            end
        end
        
    end

end

reg = {'CA3','mEC','Both'};
target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};
for ii = 1:length(reg)
    figure
    for jj = 1:length(target)
        for zz = 1:length(zone)
            subplot(2,4,4*(jj-1)+zz)
            if zz == 1 %|| zz ==3
                [~,idx] = max(PF.rateMap{ii,jj}{1,zz},[],2);
                [~,sortidx] = sort(idx,'ascend');
            end
            imagesc(zscore(PF.rateMap{ii,jj}{1,zz}(sortidx,:),0,2))
            title(strcat(target(jj),'.',zone(zz)));
        end
    end
end

