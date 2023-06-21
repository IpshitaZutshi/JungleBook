% Compile data across all sessions

function compileSplitterCellsSessAvg(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'downsample',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
downsample = p.Results.downsample;

tag = 'mEC'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};%'IZ21\Final'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
        'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'}; % To add: IZ23, IZ24, IZ25, IZ26
    reg = {'N/A','mEC','N/A'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for rr = 1:3
    for cc = 1:2
        for zz = 1:4 %4 fields, left trajectories, right trajectories, stim left and right trajectories
            PFStats.location{rr,cc}{zz} = [];
            PFStats.peakRate{rr,cc}{zz} = [];
            PFStats.avgRate{rr,cc}{zz} = [];
            PFStats.size{rr,cc}{zz} = [];
            PFStats.number{rr,cc}{zz} = [];    
            PFStats.mice{rr,cc}{zz} = [];
            PFStats.sess{rr,cc}{zz} = [];
            PFStats.region{rr,cc}{zz} = [];
            PFStats.putativeClass{rr,cc}{zz} = [];  
            PFStats.ratemap{rr,cc}{zz} = []; 
        end
        PFStats.splitter{rr,cc} = [];
    end
end

ret = round((100/175)*80);
stem = round((100/175)*110);
sessNum = 1;
%% Loop through the mice
for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');

    % Start collecting data
    for ii = 1:size(allSess,1)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);   
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);

        efields = fieldnames(sessionPulses);    
        
        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
            if ~downsample
                load(strcat(efields{jj},'\',allSess(ii).name,'.firingMapsAvg.cellinfo.mat'))
                load(strcat(efields{jj},'\',allSess(ii).name,'.placeFields.cellinfo.mat'))
            else
                load(strcat(efields{jj},'\',allSess(ii).name,'_DS.firingMapsAvg.cellinfo.mat'))
                load(strcat(efields{jj},'\',allSess(ii).name,'_DS.placeFields.cellinfo.mat'))
            end
            for kk = 1:length(placeFieldStats.mapStats)
                
                trialMaps = firingMaps.trialMaps{kk,1};
                for zz = 1:4
                    peakFR{zz} = [];
                    rm{zz} = [];
                end
                for zz = 1:size(placeFieldStats.mapStats{1},2)
                   
                  % Collect other info
                    [~,idxStem] = max(placeFieldStats.mapStats{kk,1}{1,zz}.x);                    
                    PFStats.peakRate{region,target}{zz} = [PFStats.peakRate{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.peak(idxStem)];
                    PFStats.avgRate{region,target}{zz} = [PFStats.avgRate{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.mean(idxStem)];    
                    PFStats.mice{region,target}{zz} = [PFStats.mice{region,target}{zz}; m];
                    PFStats.sess{region,target}{zz} = [PFStats.sess{region,target}{zz}; sessNum];
                    PFStats.size{region,target}{zz} = [PFStats.size{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.size(idxStem)];
                    if ~isnan(placeFieldStats.mapStats{kk,1}{1,(zz)}.x(1))
                        PFStats.number{region,target}{zz} = [PFStats.number{region,target}{zz}; 1];
                    else
                        PFStats.number{region,target}{zz} = [PFStats.number{region,target}{zz}; 0];
                    end
                    
                    if zz == 1 || zz == 3
                        [~,idxStemL] = max(placeFieldStats.mapStats{kk,1}{1,1}.x);  
                        locL = placeFieldStats.mapStats{kk,1}{1,1}.x(idxStemL);
                        [~,idxStemR] = max(placeFieldStats.mapStats{kk,1}{1,3}.x);  
                        locR = placeFieldStats.mapStats{kk,1}{1,3}.x(idxStemR);                       
                        if abs(locL-locR)>30 %If the two fields are really far apart, take the field that is closer to the
                            if locL>=60 && locL < 90
                               loc = placeFieldStats.mapStats{kk,1}{1,1}.fieldX(idxStemL,:); 
                               PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,1}.x(idxStemL)];
                            elseif locR>=60 && locR < 90
                                loc = placeFieldStats.mapStats{kk,1}{1,3}.fieldX(idxStemR,:); 
                                PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,3}.x(idxStemR)];
                            else
                                loc = nan;
                                PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};nan];
                            end
                        else % Else just stay within the field
                            loc = placeFieldStats.mapStats{kk,1}{1,zz}.fieldX(idxStem,:);
                            PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.x(idxStem)];
                        end                            
                    elseif zz == 2 || zz == 4
                        [~,idxStemL] = max(placeFieldStats.mapStats{kk,1}{1,2}.x); 
                        locL = placeFieldStats.mapStats{kk,1}{1,2}.x(idxStemL);
                        [~,idxStemR] = max(placeFieldStats.mapStats{kk,1}{1,4}.x);  
                        locR = placeFieldStats.mapStats{kk,1}{1,4}.x(idxStemR);                       
                        if abs(locL-locR)>30 %If the two fields are really far apart, take the field that is closer to the
                            if locL>=60 && locL < 90
                               loc = placeFieldStats.mapStats{kk,1}{1,2}.fieldX(idxStemL,:); 
                               PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,2}.x(idxStemL)];
                            elseif locR>=60 && locR < 90
                                loc = placeFieldStats.mapStats{kk,1}{1,4}.fieldX(idxStemR,:); 
                                PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,4}.x(idxStemR)];
                            else
                                loc = nan;
                                PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};nan];
                            end
                        else % Else just stay within the field
                            loc = placeFieldStats.mapStats{kk,1}{1,zz}.fieldX(idxStem,:);
                            PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.x(idxStem)];
                        end       
                    end
                    
                    numTrials = min([length(firingMaps.stim) length(firingMaps.arm) length(firingMaps.choice)]);
                    if zz == 1
                        idxFR = firingMaps.stim(1:numTrials)'==0 & firingMaps.arm(1:numTrials) == 1 & firingMaps.choice(1:numTrials) == 1; % ARM is opposite to firing maps - why????
                    elseif zz == 2
                        idxFR = firingMaps.stim(1:numTrials)'==1 & firingMaps.arm(1:numTrials) == 1 & firingMaps.choice(1:numTrials) == 1;
                    elseif zz == 3
                        idxFR = firingMaps.stim(1:numTrials)'==0 & firingMaps.arm(1:numTrials) == 0 & firingMaps.choice(1:numTrials) == 1;
                    elseif zz == 4
                        idxFR = firingMaps.stim(1:numTrials)'==1 & firingMaps.arm(1:numTrials) == 0 & firingMaps.choice(1:numTrials) == 1;
                    end

                    for pp = 1:length(idxFR)
                        if idxFR(pp) == 1
                            if sum(~isnan(loc))==2 && ~isempty(trialMaps{pp})
                            %if ~isnan(loc) && ~isempty(trialMaps{pp})
                                peakFR{zz} = [peakFR{zz}; nanmean(trialMaps{pp}(loc(1):loc(2)))];
                                rm{zz} = [rm{zz}; trialMaps{pp}];
                                %peakFR{zz} = [peakFR{zz}; trialMaps{pp}(loc)];
                            else
                                peakFR{zz} = [peakFR{zz};nan];
                                 rm{zz} = [rm{zz};nan(1,100)];
                            end                           
                        end                        
                    end
                    PFStats.ratemap{region,target}{zz} = [PFStats.ratemap{region,target}{zz};nanmean(rm{zz},1)];
                    % Assign numerical tag to putative class
                    if strcmp(cell_metrics.putativeCellType{kk},'Narrow Interneuron') == 1
                        PFStats.putativeClass{region,target}{zz} = [PFStats.putativeClass{region,target}{zz}; 1];
                    elseif strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1
                        PFStats.putativeClass{region,target}{zz} = [PFStats.putativeClass{region,target}{zz}; 2];
                    elseif strcmp(cell_metrics.putativeCellType{kk},'Wide Interneuron') == 1
                        PFStats.putativeClass{region,target}{zz} = [PFStats.putativeClass{region,target}{zz}; 3];
                    else 
                        PFStats.putativeClass{region,target}{zz} = [PFStats.putativeClass{region,target}{zz}; 4];                       
                    end                    

                    % Assign numerical tag to region
                    if strcmp(cell_metrics.brainRegion{kk},'CA1') == 1
                        PFStats.region{region,target}{zz} = [PFStats.region{region,target}{zz}; 1];
                    elseif strcmp(cell_metrics.brainRegion{kk},'DG') == 1
                        PFStats.region{region,target}{zz} = [PFStats.region{region,target}{zz}; 2];
                    elseif strcmp(cell_metrics.brainRegion{kk},'CA3') == 1
                        PFStats.region{region,target}{zz} = [PFStats.region{region,target}{zz}; 3];
                    else 
                        PFStats.region{region,target}{zz} = [PFStats.region{region,target}{zz}; 4];                       
                    end
                end
                if ~isempty(peakFR{1}) && sum(~isnan(peakFR{1}))>2 && ~isempty(peakFR{3}) && sum(~isnan(peakFR{3}))>2
                    pCorrBase = ranksum(peakFR{1},peakFR{3});
                else 
                    pCorrBase = nan;
                end
                if ~isempty(peakFR{2}) && sum(~isnan(peakFR{2}))>2 && ~isempty(peakFR{4}) && sum(~isnan(peakFR{4}))>2
                    pCorrStim = ranksum(peakFR{2},peakFR{4});
                else
                    pCorrStim = nan;
                end
                PFStats.splitter{region,target} = [PFStats.splitter{region,target};[pCorrBase pCorrStim]];
            end
            sessNum = sessNum+1;
        end
    end
end


%%   20(12th bin)<----____          __________------> 175 (100th bin)
%                     |             |
%                     |             |
%                     |             |
%                     |             |-----> 110 (~63rd bin)
%                     |             |
%                     |_____________|
%                              | 
%                              ^ 80 (~45th bin)         

ret = round((100/175)*80);
stem = round((100/175)*110);

target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};

colMat = [0 0 0;... %Black
    224/243 163/243 46/243;...  %Yellow    
    8/243 133/243 161/243;...   %Teal    
    56/243 61/243 150/243;...
    180/243 180/243 180/243];     %Gray
%Plot differently depending on whether plotting mEC only, or manipulations
%with separate analog channels, or CA3

%% If plotting CA1
if strcmp(tag,'CA1') == 1 || strcmp(tag,'mECBilateral') == 1 || strcmp(tag,'mEC') == 1 
    
    figure
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1200 700])
    
    manipNum = [];
    datPeak =[];
    datMean = [];
    datSplitter = [];
    datLoc = [];
    datManip = [];
    datManipID = [];
    datPeakID = [];
    manipNumID = [];
    datPyr = [];
    
    if strcmp(tag,'mEC') == 1 
        iiRange = 2;
        splitPercB{2} =[];
        splitPercS{2} =[];
    else
        iiRange = [1 2 3];
        for nn = 1:3
            splitPercB{nn} =[];
            splitPercS{nn} =[];
        end
    end

    for ii = iiRange
        for jj = 1%:length(target)            

            % Separate CA1, select principal cells with place fields in any
            % of the 4 conditions
            idxSelect = PFStats.putativeClass{ii,jj}{1,1} == 2 & PFStats.region{ii,jj}{1,1} == 1 &...
                (PFStats.number{ii,jj}{1,1} == 1 | PFStats.number{ii,jj}{1,2} == 1 | PFStats.number{ii,jj}{1,3} == 1 | PFStats.number{ii,jj}{1,4} == 1);                

            for zz = 1:2        

                %PlaceFieldidx_curr = idxSelect;
                PlaceFieldidx_curr = PFStats.putativeClass{ii,jj}{1,zz} == 2 & PFStats.region{ii,jj}{1,zz} == 1 &...
                (PFStats.number{ii,jj}{1,zz} == 1 | PFStats.number{ii,jj}{1,zz+2} == 1);               
                PyrCellidx = PFStats.putativeClass{ii,jj}{1,zz} == 2 & PFStats.region{ii,jj}{1,zz} == 1; 

                locationA = PFStats.location{ii,jj}{zz};
                locationB = PFStats.location{ii,jj}{zz+2};
                locMean = round(nanmean([locationA locationB],2));
                difflocAB = abs(locationA-locationB);
                idxSplitter = locMean >=stem & locMean <=90;% & difflocAB <= 30;

                location1 = PFStats.location{ii,jj}{1};
                location2 = PFStats.location{ii,jj}{2};
                location3 = PFStats.location{ii,jj}{3};
                location4 = PFStats.location{ii,jj}{4};                
                MaxLoc = max([location1 location2 location3 location4],[],2);
                MinLoc = min([location1 location2 location3 location4],[],2);
                diffloc1234 = MaxLoc-MinLoc;
                idxSplitterID = location1>=stem & location2>=stem & location3>=stem & location4>=stem &...
                    location1<=90 & location2<=90 & location3<=90 & location4<=90 & diffloc1234 <=10;

                idxSess = unique(PFStats.sess{ii,jj}{1});
                tempSplitter = PFStats.splitter{ii,jj}(:,zz);
                
                for nn = 1:length(idxSess)
                    if zz ==1
                        numSplitter = PFStats.sess{ii,jj}{1}==idxSess(nn) & idxSplitter & tempSplitter <= 0.05 & PlaceFieldidx_curr;
                        numPlace = PFStats.sess{ii,jj}{1}==idxSess(nn) & idxSplitter & PlaceFieldidx_curr;
                        a = [sum(numSplitter)./sum(numPlace) sum(numPlace)];
                        splitPercB{ii} = [splitPercB{ii}; a];
                    else
                        numSplitter = PFStats.sess{ii,jj}{1}==idxSess(nn) & idxSplitter & tempSplitter <= 0.05 & PlaceFieldidx_curr;
                        numPlace = PFStats.sess{ii,jj}{1}==idxSess(nn) & idxSplitter & PlaceFieldidx_curr;
                        a = [sum(numSplitter)./sum(numPlace) sum(numPlace)];
                        splitPercS{ii} = [splitPercS{ii}; a];
                    end
                end
            end
        end
    end  
    
    if strcmp(tag,'mEC') == 1 
        finalDist = [splitPercB{2}(:,1) splitPercS{2}(:,1)];
        stats = groupStats({finalDist(:,1),finalDist(:,2)},[],'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM');
        ylim([0 0.7])

    else
        splitPercB{3}(2,:) = nan;
        splitPercS{3}(2,:) = nan;
        finalDist{1} = [splitPercB{1}(:,1) splitPercS{1}(:,1)];
        finalDist{2} = [splitPercB{2}(:,1) splitPercS{2}(:,1)];
        finalDist{3} = [splitPercB{3}(:,1) splitPercS{3}(:,1)];
    
        subplot(1,3,1)
        stats{1} = groupStats({finalDist{1}(:,1),finalDist{1}(:,2)},[],'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM');
        ylim([0 0.7])

        subplot(1,3,2)
        stats{2} = groupStats({finalDist{2}(:,1),finalDist{2}(:,2)},[],'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM');
        ylim([0 0.7])
        
        subplot(1,3,3)
        stats{3} = groupStats({finalDist{3}(:,1),finalDist{3}(:,2)},[],'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM');
        ylim([0 0.7])
    end
    if ~downsample
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsSummary',tag,'.png'));
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsSummary',tag,'.eps'),'epsc');
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsSummary',tag,'.fig'));
        save(strcat(parentDir,'Compiled\Place Fields\SplitterCellsSummary',tag,'.mat'),'stats')  
    else
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsDS',tag,'.png'));
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsDS',tag,'.eps'),'epsc');
        saveas(gcf,strcat(parentDir,'Compiled\Place Fields\SplitterCellsDS',tag,'.fig'));
    end
              
end
    
end

