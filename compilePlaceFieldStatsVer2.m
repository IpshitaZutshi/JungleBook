% Compile data across all sessions

function compilePlaceFieldStatsVer2(varargin)

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

tag = 'CA3'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};%'IZ21\Final'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
        'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'}; % To add: IZ23, IZ24, IZ25, IZ26
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};%'IZ27\Saline','IZ28\Saline','IZ29\Saline',
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
%             PFStats.information1{rr,cc}{zz} = [];
%             PFStats.information2{rr,cc}{zz} = [];
%             PFStats.sparsity{rr,cc}{zz} = [];           
%             PFStats.selectivity{rr,cc}{zz} = [];             
            PFStats.central.information1{rr,cc}{zz} = [];
            PFStats.central.information2{rr,cc}{zz} = [];
            PFStats.central.sparsity{rr,cc}{zz} = [];           
            PFStats.central.selectivity{rr,cc}{zz} = [];  
            PFStats.side.information1{rr,cc}{zz} = [];
            PFStats.side.information2{rr,cc}{zz} = [];
            PFStats.side.sparsity{rr,cc}{zz} = [];           
            PFStats.side.selectivity{rr,cc}{zz} = [];  
            PFStats.side.stability{rr,cc}{zz} = [];  
            PFStats.side.rate{rr,cc}{zz} = [];
            PFStats.central.stability{rr,cc}{zz} = [];  
            PFStats.central.rate{rr,cc}{zz} = [];             
            PFStats.mice{rr,cc}{zz} = [];
            PFStats.region{rr,cc}{zz} = [];
            PFStats.putativeClass{rr,cc}{zz} = [];    
        end
    end
end

ret = round((100/175)*80);
stem = round((100/175)*110);

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
                for zz = 1:size(placeFieldStats.mapStats{1},2)
                    
                    Rate_Map = firingMaps.rateMaps{kk,1}{1,zz};
                    snSpikes = firingMaps.countMaps{kk,1}{1,zz};
                    sTimeSpent = firingMaps.occupancy{kk,1}{1,zz};
                    
                    %get place field stats
%                     [Information_1,Information_2,Sparsity,~,Selectivity] = PlaceCellInfo(Rate_Map, snSpikes, sTimeSpent);
%                     PFStats.information1{region,target}{zz} = [PFStats.information1{region,target}{zz};Information_1];
%                     PFStats.information2{region,target}{zz} = [PFStats.information2{region,target}{zz};Information_2];
%                     PFStats.sparsity{region,target}{zz} = [PFStats.sparsity{region,target}{zz};Sparsity];           
%                     PFStats.selectivity{region,target}{zz} = [PFStats.selectivity{region,target}{zz};Selectivity];                      
                    
                    [Information_1,Information_2,Sparsity,~,Selectivity] = PlaceCellInfo(Rate_Map(1:ret), snSpikes(1:ret), sTimeSpent(1:ret));
                    PFStats.side.information1{region,target}{zz} = [PFStats.side.information1{region,target}{zz};Information_1];
                    PFStats.side.information2{region,target}{zz} = [PFStats.side.information2{region,target}{zz};Information_2];
                    PFStats.side.sparsity{region,target}{zz} = [PFStats.side.sparsity{region,target}{zz};Sparsity];           
                    PFStats.side.selectivity{region,target}{zz} = [PFStats.side.selectivity{region,target}{zz};Selectivity];  
                     
                    [Information_1,Information_2,Sparsity,~,Selectivity] = PlaceCellInfo(Rate_Map(stem:95), snSpikes(stem:95), sTimeSpent(stem:95));
                    PFStats.central.information1{region,target}{zz} = [PFStats.central.information1{region,target}{zz};Information_1];
                    PFStats.central.information2{region,target}{zz} = [PFStats.central.information2{region,target}{zz};Information_2];
                    PFStats.central.sparsity{region,target}{zz} = [PFStats.central.sparsity{region,target}{zz};Sparsity];           
                    PFStats.central.selectivity{region,target}{zz} = [PFStats.central.selectivity{region,target}{zz};Selectivity];  
                    
                    PFStats.side.rate{region,target}{zz} = [PFStats.side.rate{region,target}{zz};sum(snSpikes(1:ret))./sum(sTimeSpent(1:ret))];
                    PFStats.central.rate{region,target}{zz} = [PFStats.central.rate{region,target}{zz};sum(snSpikes(stem:95))./sum(sTimeSpent(stem:95))];
                   
                    % get place field stability
                    if ~isempty(firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+1})&& ~isempty(firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+2})
                        corrMap = corr(firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+1}(1:ret)',firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+2}(1:ret)','Rows','pairwise','Type','Spearman');
                        PFStats.side.stability{region,target}{zz} = [PFStats.side.stability{region,target}{zz};corrMap];
                        corrMap = corr(firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+1}(stem:95)',firingMaps.sessCorrMaps{kk,1}{1,2*(zz-1)+2}(stem:95)','Rows','pairwise','Type','Spearman');
                        PFStats.central.stability{region,target}{zz} = [PFStats.central.stability{region,target}{zz};corrMap];                    
                    else
                        PFStats.side.stability{region,target}{zz} = [PFStats.side.stability{region,target}{zz};nan];
                        PFStats.central.stability{region,target}{zz} = [PFStats.central.stability{region,target}{zz};nan];  
                    end
                    % Collect other info
                    PFStats.location{region,target}{zz} = [PFStats.location{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.x(1)];
                    PFStats.peakRate{region,target}{zz} = [PFStats.peakRate{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.peak(1)];
                    PFStats.avgRate{region,target}{zz} = [PFStats.avgRate{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.mean(1)];    
                    PFStats.mice{region,target}{zz} = [PFStats.mice{region,target}{zz}; m];
                    PFStats.size{region,target}{zz} = [PFStats.size{region,target}{zz};placeFieldStats.mapStats{kk,1}{1,zz}.size(1)];
                    if ~isnan(placeFieldStats.mapStats{kk,1}{1,(zz)}.x(1))
                        PFStats.number{region,target}{zz} = [PFStats.number{region,target}{zz}; 1];
                    else
                        PFStats.number{region,target}{zz} = [PFStats.number{region,target}{zz}; 0];
                    end
                    
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
            end
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

colMat = [85/243 85/243 85/243;... %Gray
    224/243 163/243 46/243;...  %Yellow    
    8/243 133/243 161/243;...   %Teal    
    56/243 61/243 150/243];     %Navy           

%Plot differently depending on whether plotting mEC only, or manipulations
%with separate analog channels, or CA3

%% If plotting CA1
if strcmp(tag,'CA1') == 1 || strcmp(tag,'mECBilateral') == 1 
    
    figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1500 700])
    figure(2)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1800 500]) 
    
    manipNum = [];
    
    locName = {'side','central'};
    for loc = 1:length(locName)
        datPeak.(locName{loc}) =[];
        datMean.(locName{loc}) = [];
        datSize.(locName{loc}) = [];
        datNum.(locName{loc}) = [];
        manip.(locName{loc}) = [];
        
        datInfo1.(locName{loc}) = [];
        datInfo2.(locName{loc}) = [];
        datSparsity.(locName{loc}) = [];
        datSelectivity.(locName{loc}) = [];
        datStability.(locName{loc}) = [];
        datManip.(locName{loc}) = [];
        datManip2.(locName{loc}) = [];
    end
    rateSide = [];
    rateCenter = [];
    manipRate = [];   
    for ii = 1:length(reg)
        for jj = 1%:length(target)
            
            % Separate CA1, select principal cells with place fields in any
            % of the 4 conditions
            idxSelect = PFStats.putativeClass{ii,jj}{1,1} == 2 & PFStats.region{ii,jj}{1,1} == 1;% &...
 %               (PFStats.number{ii,jj}{1,1} == 1 | PFStats.number{ii,jj}{1,2} == 1 | PFStats.number{ii,jj}{1,3} == 1 | PFStats.number{ii,jj}{1,4} == 1);                
            
            for zz = 1:length(zone)             
%                 if zz == 1 || zz == 3
%                     location1 = PFStats.location{ii,jj}{zz};
%                     location2 = PFStats.location{ii,jj}{zz+1};
%                     %Make sure only overlapping fields are considered
%                     idxOverlap = abs(location1-location2) <= 10;
%                     if jj == 1 % If target is stem
%                         idxstim = location1 >= stem  & idxOverlap & idxSelect;
%                         idxnostim = location1 <=stem & idxOverlap & idxSelect;
%                     elseif jj == 2 % If target is return
%                         idxstim = location1<= ret & idxOverlap & idxSelect;
%                         idxnostim = location1 >= ret & idxOverlap & idxSelect;
%                     end                            
%                     
%                     %% Plot peak rate
%                     figure(1)
%                     subplot(3,5,5*(ii-1)+1)
%                     s = scatter(PFStats.peakRate{ii,jj}{zz}(idxnostim),PFStats.peakRate{ii,jj}{zz+1}(idxnostim),4,col(1,:),'filled');
%                     %s.MarkerFaceAlpha = 0.7;
%                     xlim([0 40])
%                     ylim([0 40])
%                     xlabel('Peak Rate Baseline')
%                     ylabel('Peak Rate Stim')
% 
%                     hold on
%                     s = scatter(PFStats.peakRate{ii,jj}{zz}(idxstim),PFStats.peakRate{ii,jj}{zz+1}(idxstim),4,col(ii+1,:),'filled');
%                     %s.MarkerFaceAlpha = 0.7;
%                     xlim([0 40])
%                     ylim([0 40])
%                     xlabel('Peak Rate Baseline')
%                     ylabel('Peak Rate Stim')
%                     title(strcat('Region: ',reg{ii}))
%                     refline(1,0)
%                 end
%                 

                if zz == 1 || zz == 3
                    manipNum = [manipNum;1];
                else
                    manipNum = [manipNum;(ii+1)];
                end
                
                figure(2)                
                if zz == 1 || zz == 3
                    subplot(2,9,ii)
                    scatter(PFStats.side.rate{ii,jj}{zz}(idxSelect),PFStats.side.rate{ii,jj}{zz+1}(idxSelect),4,colMat(ii+1,:),'filled')
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Side arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])
                    rateSide = [rateSide [PFStats.side.rate{ii,jj}{zz}(idxSelect)';PFStats.side.rate{ii,jj}{zz+1}(idxSelect)']];
                    rateCenter = [rateCenter [PFStats.central.rate{ii,jj}{zz}(idxSelect)';PFStats.central.rate{ii,jj}{zz+1}(idxSelect)']];
                    manipRate = [manipRate ones(1,sum(idxSelect))*ii];
                        
                    subplot(2,9,9+ii)
                    scatter(PFStats.central.rate{ii,jj}{zz}(idxSelect),PFStats.central.rate{ii,jj}{zz+1}(idxSelect),4,colMat(ii+1,:),'filled')
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Central arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])                    
                end
 
           idxSelect = PFStats.putativeClass{ii,jj}{1,1} == 2 & PFStats.region{ii,jj}{1,1} == 1 &...
                (PFStats.number{ii,jj}{1,1} == 1 | PFStats.number{ii,jj}{1,2} == 1 | PFStats.number{ii,jj}{1,3} == 1 | PFStats.number{ii,jj}{1,4} == 1);                
                 
                PlaceFieldidx_curr = idxSelect;%PFStats.putativeClass{ii,jj}{1,zz} == 2 & PFStats.region{ii,jj}{1,zz} == 1 & PFStats.number{ii,jj}{1,zz} == 1; 
                PyrCellidx = PFStats.putativeClass{ii,jj}{1,zz} == 2 & PFStats.region{ii,jj}{1,zz} == 1; 
                
                location = PFStats.location{ii,jj}{zz};
                fieldLocation = zeros(length(location),1);
                fieldLocation(location >= stem) = 2;
                fieldLocation(location < stem) = 1;
                                                
                for loc = 1:length(locName)
                     
                    datPeak.(locName{loc}) = [datPeak.(locName{loc});PFStats.peakRate{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datMean.(locName{loc}) = [datMean.(locName{loc});PFStats.avgRate{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datSize.(locName{loc}) = [datSize.(locName{loc});PFStats.size{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datInfo1.(locName{loc}) = [datInfo1.(locName{loc});PFStats.information1{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datInfo2.(locName{loc}) = [datInfo2.(locName{loc});PFStats.information2{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datSparsity.(locName{loc}) = [datSparsity.(locName{loc});PFStats.sparsity{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datSelectivity.(locName{loc}) = [datSelectivity.(locName{loc});PFStats.selectivity{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datInfo1.(locName{loc}) = [datInfo1.(locName{loc});PFStats.(locName{loc}).information1{ii,jj}{zz}(PyrCellidx)];
                    datInfo2.(locName{loc}) = [datInfo2.(locName{loc});PFStats.(locName{loc}).information2{ii,jj}{zz}(PyrCellidx)];
                    datSparsity.(locName{loc}) = [datSparsity.(locName{loc});PFStats.(locName{loc}).sparsity{ii,jj}{zz}(PyrCellidx)];
                    datSelectivity.(locName{loc}) = [datSelectivity.(locName{loc});PFStats.(locName{loc}).selectivity{ii,jj}{zz}(PyrCellidx)];

                    datNum.(locName{loc}) = [datNum.(locName{loc});sum(PFStats.number{ii,jj}{zz}(fieldLocation==loc))./numel(PFStats.number{ii,jj}{zz}(PyrCellidx))];                                  
                    datStability.(locName{loc}) = [datStability.(locName{loc});PFStats.(locName{loc}).stability{ii,jj}{zz}(PyrCellidx)];
                    if zz == 1 || zz == 3
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*1];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*1];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*1];
                    else
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*(ii+1)];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*(ii+1)];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*(ii+1)];
                    end

                end                
            end            
        end
    end  

    %% Now deal with other data
    %Fig 1 peak Rate, mean Rate, bits/cm
    
    %Rate, infor1 sparsity selectivity
    %stability
    for loc = 1:length(locName)
            
        figure(1)        
        subplot(2,4,4*(loc-1)+1)
        stats.(locName{loc}).peakRate = groupStats(datPeak.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Peak Rate (Hz)')
        
        subplot(2,4,4*(loc-1)+2)
        stats.(locName{loc}).meanRate = groupStats(datMean.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Mean Rate (Hz)')        

        subplot(2,4,4*(loc-1)+3)
        stats.(locName{loc}).info2 = groupStats(datInfo2.(locName{loc}),datManip2.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('information (bits/sec)')    
        
        subplot(2,4,4*(loc-1)+4)
        stats.(locName{loc}).Selectivity = groupStats(datSelectivity.(locName{loc}),datManip2.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Selectivity')       
        
        figure(2)  
        subplot(2,9,9*(loc-1)+5)
        data = [nanmean(datNum.(locName{loc})(manipNum==1)) nanmean(datNum.(locName{loc})(manipNum==2)) nanmean(datNum.(locName{loc})(manipNum==3)) nanmean(datNum.(locName{loc})(manipNum==4))];
        bar(data,'k')
        ylabel('% of place fields')
        title(num2str(data))
        
        subplot(2,9,9*(loc-1)+6)
        stats.(locName{loc}).size = groupStats(datSize.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Field Size')

        subplot(2,9,9*(loc-1)+7)
        stats.(locName{loc}).info1 = groupStats(datInfo1.(locName{loc}),datManip2.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('information (bits/spike)')
                
        subplot(2,9,9*(loc-1)+8)
        stats.(locName{loc}).Sparsity = groupStats(datSparsity.(locName{loc}),datManip2.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Sparsity')
               
        subplot(2,9,9*(loc-1)+9)
        stats.(locName{loc}).Stability = groupStats(datStability.(locName{loc}),datManip2.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Stability')         
    end
    
    figure(2)
    for ii = 1:3
%         side.(reg{ii}) = (point_to_line_distance(rateSide(:,manipRate == ii)',[0,0],[7,7]));
%         idxminus = rateSide(1,manipRate == ii) > rateSide(2,manipRate == ii); 
%        side.(reg{ii})(idxminus) = rateSide(2,manipRate == ii)./rateSide(1,manipRate == ii);%-1*side.(reg{ii})(idxminus);
%         center.(reg{ii}) = (point_to_line_distance(rateCenter(:,manipRate == ii)',[0,0],[7,7]));
%         idxminus = rateCenter(1,manipRate == ii) > rateCenter(2,manipRate == ii); 
%        center.(reg{ii})(idxminus) = -1*center.(reg{ii})(idxminus);   
        side.(reg{ii}) = rateSide(2,manipRate == ii)./rateSide(1,manipRate == ii);
       %side.(reg{ii}) = side.(reg{ii})(side.(reg{ii})<100);
        center.(reg{ii}) = rateCenter(2,manipRate == ii)./rateCenter(1,manipRate == ii);
        %center.(reg{ii}) = center.(reg{ii})(center.(reg{ii})<100);
        subplot(2,9,4)  
        a = histcounts(side.(reg{ii}),[0:0.15:3],'Normalization','probability');
        plot(0:0.15:2.85,smooth(a),'color',colMat(ii+1,:),'LineWidth',1.5)
        xlim([0 3])
        hold on
        
        subplot(2,9,9+4)  
        a = histcounts(center.(reg{ii}),[0:0.15:3],'Normalization','probability');
        plot(0:0.15:2.85,smooth(a),'color',colMat(ii+1,:),'LineWidth',1.5)
        xlim([0 3])
        hold on      
    end
%     subplot(2,9,4)    
%     nhist(side,'samebins','smooth','proportion','binfactor',3)
%     xlim([0 4])
%     subplot(2,9,9+4)    
%     nhist(center,'samebins','smooth','proportion','binfactor',3)
%     xlim([0 4])
    
    if ~downsample
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.fig'));

        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.png'));
        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.eps'),'epsc');
        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.fig'));
        save(strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'.mat'),'stats');
    else
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1_DS.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1_DS.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1_DS.fig'));

        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2_DS.png'));
        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2_DS.eps'),'epsc');
        saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2_DS.fig'));
        save(strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'_DS.mat'),'stats');
    end

%% If plotting mEC    
elseif strcmp(tag,'mEC') == 1
    
    
%% If plotting CA3    
elseif strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    
colMat = [85/243 85/243 85/243;... %Gray  
    8/243 133/243 161/243;...   %Teal   
    224/243 163/243 46/243;...  %Yellow      
    56/243 61/243 150/243];     %Navy      

     figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1500 700])
    figure(2)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1800 500]) 
    
    manipNum = [];
    
    locName = {'side','central'};
    for loc = 1:length(locName)
        datPeak.(locName{loc}) =[];
        datMean.(locName{loc}) = [];
        datSize.(locName{loc}) = [];
        datNum.(locName{loc}) = [];
        manip.(locName{loc}) = [];
        
        datInfo1.(locName{loc}) = [];
        datInfo2.(locName{loc}) = [];
        datSparsity.(locName{loc}) = [];
        datSelectivity.(locName{loc}) = [];
        datStability.(locName{loc}) = [];
        datManip.(locName{loc}) = [];
        datManip2.(locName{loc}) = [];
    end
    rateSide = [];
    rateCenter = [];
    manipRate = []; 
    for ii = 2:length(reg)
        for jj = 1%:length(target)
            
            % Separate CA1, select principal cells with place fields in any
            % of the 4 conditions
            idxSelect = PFStats.putativeClass{ii,jj}{1,1} == 1 & PFStats.region{ii,jj}{1,1} == 1;% &...
%                 (PFStats.number{2,jj}{1,1} == 1 | PFStats.number{2,jj}{1,2} == 1 | PFStats.number{2,jj}{1,3} == 1 | PFStats.number{2,jj}{1,4} == 1 |...
%                 PFStats.number{3,jj}{1,1} == 1 | PFStats.number{3,jj}{1,2} == 1 | PFStats.number{3,jj}{1,3} == 1 | PFStats.number{3,jj}{1,4} == 1);                
%             
            for zz = 1:length(zone)             
%                 if zz == 1 || zz == 3
%                     location1 = PFStats.location{ii,jj}{zz};
%                     location2 = PFStats.location{ii,jj}{zz+1};
%                     %Make sure only overlapping fields are considered
%                     idxOverlap = abs(location1-location2) <= 10;
%                     if jj == 1 % If target is stem
%                         idxstim = location1 >= stem  & idxOverlap & idxSelect;
%                         idxnostim = location1 <=stem & idxOverlap & idxSelect;
%                     elseif jj == 2 % If target is return
%                         idxstim = location1<= ret & idxOverlap & idxSelect;
%                         idxnostim = location1 >= ret & idxOverlap & idxSelect;
%                     end                            
%                     
%                     %% Plot peak rate
%                     figure(1)
%                     subplot(3,5,5*(ii-1)+1)
%                     s = scatter(PFStats.peakRate{ii,jj}{zz}(idxnostim),PFStats.peakRate{ii,jj}{zz+1}(idxnostim),4,col(1,:),'filled');
%                     %s.MarkerFaceAlpha = 0.7;
%                     xlim([0 40])
%                     ylim([0 40])
%                     xlabel('Peak Rate Baseline')
%                     ylabel('Peak Rate Stim')
% 
%                     hold on
%                     s = scatter(PFStats.peakRate{ii,jj}{zz}(idxstim),PFStats.peakRate{ii,jj}{zz+1}(idxstim),4,col(ii+1,:),'filled');
%                     %s.MarkerFaceAlpha = 0.7;
%                     xlim([0 40])
%                     ylim([0 40])
%                     xlabel('Peak Rate Baseline')
%                     ylabel('Peak Rate Stim')
%                     title(strcat('Region: ',reg{ii}))
%                     refline(1,0)
%                 end
%                 

                if ii == 2 && (zz == 1 || zz == 3)
                    manipNum = [manipNum;1];
                elseif ii == 2 && (zz == 2 || zz == 4)
                    manipNum = [manipNum;2];
                elseif ii == 3 && (zz == 1 || zz == 3)
                    manipNum = [manipNum;3];
                elseif ii == 3 && (zz == 2 || zz == 4)
                    manipNum = [manipNum;4];                    
                end
                
                figure(2)                
                if ii == 2 && (zz == 1 || zz == 3)
                    
                    rateSide = [rateSide [PFStats.side.rate{ii,jj}{zz}(idxSelect)';PFStats.side.rate{ii,jj}{zz+1}(idxSelect)']];
                    rateCenter = [rateCenter [PFStats.central.rate{ii,jj}{zz}(idxSelect)';PFStats.central.rate{ii,jj}{zz+1}(idxSelect)']];
                    manipRate = [manipRate ones(1,sum(idxSelect))*1];
                    
                    subplot(2,9,2)
                    s = scatter(PFStats.side.rate{2,jj}{zz}(idxSelect),PFStats.side.rate{ii,jj}{zz+1}(idxSelect),4,colMat(3,:),'filled');
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Side arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])
                    
                    subplot(2,9,9+2)
                    s = scatter(PFStats.central.rate{2,jj}{zz}(idxSelect),PFStats.central.rate{ii,jj}{zz+1}(idxSelect),4,colMat(3,:),'filled');
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Central arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])    
                    
                elseif ii ==3
                    
                    if zz == 1 || zz == 3
                       subplot(2,9,1)
                       scatCol = colMat(2,:);
                       comp = zz;
                       manipIdx = 2;
                    else
                       subplot(2,9,3)
                       scatCol = colMat(4,:);
                       comp = zz-1;
                       manipIdx = 3;
                    end
                    s = scatter(PFStats.side.rate{2,jj}{comp}(idxSelect),PFStats.side.rate{ii,jj}{zz}(idxSelect),4,scatCol,'filled');
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Side arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])
                    
                    if zz == 1 || zz == 3
                       subplot(2,9,10)
                       scatCol = colMat(2,:);
                    else
                       subplot(2,9,12)
                       scatCol = colMat(4,:);
                    end
                    s = scatter(PFStats.central.rate{2,jj}{comp}(idxSelect),PFStats.central.rate{ii,jj}{zz}(idxSelect),4,scatCol,'filled');
                    xlabel('Rate Baseline')
                    ylabel('Rate Stim')
                    title('Central arm')
                    refline(1)
                    hold on
                    xlim([0 7])
                    ylim([0 7])   
                    
                    rateSide = [rateSide [PFStats.side.rate{2,jj}{comp}(idxSelect)';PFStats.side.rate{ii,jj}{zz}(idxSelect)']];
                    rateCenter = [rateCenter [PFStats.central.rate{2,jj}{comp}(idxSelect)';PFStats.central.rate{ii,jj}{zz}(idxSelect)']];
                    manipRate = [manipRate ones(1,sum(idxSelect))*manipIdx];
        
                end
                

            idxSelect = PFStats.putativeClass{ii,jj}{1,1} == 1 & PFStats.region{ii,jj}{1,1} == 1;% &...
%                  (PFStats.number{2,jj}{1,1} == 1 | PFStats.number{2,jj}{1,2} == 1 | PFStats.number{2,jj}{1,3} == 1 | PFStats.number{2,jj}{1,4} == 1 |...
%                  PFStats.number{3,jj}{1,1} == 1 | PFStats.number{3,jj}{1,2} == 1 | PFStats.number{3,jj}{1,3} == 1 | PFStats.number{3,jj}{1,4} == 1);                
                PlaceFieldidx_curr = idxSelect;                  
%
                PyrCellidx = PFStats.putativeClass{ii,jj}{1,zz} == 1 & PFStats.region{ii,jj}{1,zz} == 1; 
                
                location = PFStats.location{ii,jj}{zz};
                fieldLocation = zeros(length(location),1);
                fieldLocation(location >= stem) = 2;
                fieldLocation(location < stem) = 1;
                                                
                for loc = 1:length(locName)
                     
                    datPeak.(locName{loc}) = [datPeak.(locName{loc});PFStats.peakRate{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datMean.(locName{loc}) = [datMean.(locName{loc});PFStats.avgRate{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datSize.(locName{loc}) = [datSize.(locName{loc});PFStats.size{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    
%                     datInfo1.(locName{loc}) = [datInfo1.(locName{loc});PFStats.information1{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datInfo2.(locName{loc}) = [datInfo2.(locName{loc});PFStats.information2{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datSparsity.(locName{loc}) = [datSparsity.(locName{loc});PFStats.sparsity{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
%                     datSelectivity.(locName{loc}) = [datSelectivity.(locName{loc});PFStats.selectivity{ii,jj}{zz}(PlaceFieldidx_curr & (fieldLocation==loc))];
                    datInfo1.(locName{loc}) = [datInfo1.(locName{loc});PFStats.(locName{loc}).information1{ii,jj}{zz}(PyrCellidx)];
                    datInfo2.(locName{loc}) = [datInfo2.(locName{loc});PFStats.(locName{loc}).information2{ii,jj}{zz}(PyrCellidx)];
                    datSparsity.(locName{loc}) = [datSparsity.(locName{loc});PFStats.(locName{loc}).sparsity{ii,jj}{zz}(PyrCellidx)];
                    datSelectivity.(locName{loc}) = [datSelectivity.(locName{loc});PFStats.(locName{loc}).selectivity{ii,jj}{zz}(PyrCellidx)];

                    datNum.(locName{loc}) = [datNum.(locName{loc});sum(PFStats.number{ii,jj}{zz}(fieldLocation==loc))./numel(PFStats.number{ii,jj}{zz}(PyrCellidx))];                                  
                    datStability.(locName{loc}) = [datStability.(locName{loc});PFStats.(locName{loc}).stability{ii,jj}{zz}(PyrCellidx)];
                    if ii == 2 && (zz == 1 || zz == 3)
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*1];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*1];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*1];
                    elseif ii == 2 && (zz == 2 || zz == 4)
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*2];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*2];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*2];
                    elseif ii == 3 && (zz == 1 || zz == 3)
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*3];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*3];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*3];
                    elseif ii == 3 && (zz == 2 || zz == 4)
                        datManip.(locName{loc}) = [datManip.(locName{loc});ones(sum(PlaceFieldidx_curr),1)*4];
                        datManip2.(locName{loc}) = [datManip2.(locName{loc});ones(sum(PyrCellidx),1)*4];
                        manip.(locName{loc}) = [manip.(locName{loc});ones(sum(PlaceFieldidx_curr & (fieldLocation==loc)),1)*4];
                        
                    end

                end                
            end            
        end
    end  

    %% Now deal with other data
    %Fig 1 peak Rate, mean Rate, bits/cm
    
    %Rate, infor1 sparsity selectivity
    %stability
    for loc = 1:length(locName)
            
        figure(1)        
        subplot(2,4,4*(loc-1)+1)
        stats.(locName{loc}).peakRate = groupStats(datPeak.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Peak Rate (Hz)')
        
        subplot(2,4,4*(loc-1)+2)
        stats.(locName{loc}).meanRate = groupStats(datMean.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Mean Rate (Hz)')        

        subplot(2,4,4*(loc-1)+3)
        stats.(locName{loc}).info2 = groupStats(datInfo2.(locName{loc}),datManip2.(locName{loc}),'repeatedMeasures',true,'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('information (bits/sec)')    
        
        subplot(2,4,4*(loc-1)+4)
        stats.(locName{loc}).Selectivity = groupStats(datSelectivity.(locName{loc}),datManip2.(locName{loc}),'repeatedMeasures',true,'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Selectivity')       
        
        figure(2)  
        subplot(2,9,9*(loc-1)+5)
        data = [nanmean(datNum.(locName{loc})(manipNum==1)) nanmean(datNum.(locName{loc})(manipNum==2)) nanmean(datNum.(locName{loc})(manipNum==3)) nanmean(datNum.(locName{loc})(manipNum==4))];
        bar(data,'k')
        ylabel('% of place fields')
        title(num2str(data))
        
        subplot(2,9,9*(loc-1)+6)
        stats.(locName{loc}).size = groupStats(datSize.(locName{loc}),manip.(locName{loc}),'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Field Size')

        subplot(2,9,9*(loc-1)+7)
        stats.(locName{loc}).info1 = groupStats(datInfo1.(locName{loc}),datManip2.(locName{loc}),'repeatedMeasures',true,'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('information (bits/spike)')
                
        subplot(2,9,9*(loc-1)+8)
        stats.(locName{loc}).Sparsity = groupStats(datSparsity.(locName{loc}),datManip2.(locName{loc}),'repeatedMeasures',true,'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Sparsity')
               
        subplot(2,9,9*(loc-1)+9)
        stats.(locName{loc}).Stability = groupStats(datStability.(locName{loc}),datManip2.(locName{loc}),'repeatedMeasures',true,'inAxis',true,'color',colMat);
        title(locName{loc})
        ylabel('Stability')         
    end
    
    figure(2)
    for ii = 1:3
%         side.(reg{ii}) = (point_to_line_distance(rateSide(:,manipRate == ii)',[0,0],[7,7]));
%         idxminus = rateSide(1,manipRate == ii) > rateSide(2,manipRate == ii); 
%         side.(reg{ii})(idxminus) = -1*side.(reg{ii})(idxminus);
%         center.(reg{ii}) = (point_to_line_distance(rateCenter(:,manipRate == ii)',[0,0],[7,7]));
%         idxminus = rateCenter(1,manipRate == ii) > rateCenter(2,manipRate == ii); 
%         center.(reg{ii})(idxminus) = -1*center.(reg{ii})(idxminus);        
        side.(reg{ii}) = rateSide(2,manipRate == ii)./rateSide(1,manipRate == ii);
        %side.(reg{ii}) = side.(reg{ii})(side.(reg{ii})<100);
        center.(reg{ii}) = rateCenter(2,manipRate == ii)./rateCenter(1,manipRate == ii);
        %center.(reg{ii}) = center.(reg{ii})(center.(reg{ii})<100);
        subplot(2,9,4)  
        a = histcounts(side.(reg{ii}),[0:0.15:3],'Normalization','probability');
        plot(0:0.15:2.85,smooth(a),'color',colMat(ii+1,:),'LineWidth',1.5)
        xlim([0 3])
        hold on
        
        subplot(2,9,9+4)  
        a = histcounts(center.(reg{ii}),[0:0.15:3],'Normalization','probability');
        plot(0:0.15:2.85,smooth(a),'color',colMat(ii+1,:),'LineWidth',1.5)
        xlim([0 3])
        hold on              
    end
%     subplot(2,9,4)    
%     nhist(side,'samebins','smooth','proportion','binfactor',3)
%     xlim([0 4])
%     subplot(2,9,9+4)    
%     nhist(center,'samebins','smooth','proportion','binfactor',3)
%     xlim([0 4])
    
%     if ~downsample
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.png'));
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.eps'),'epsc');
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1.fig'));
% 
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.png'));
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.eps'),'epsc');
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2.fig'));
%         save(strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'.mat'),'stats');
%     else
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1DS.png'));
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1DS.eps'),'epsc');
%         saveas(figure(1),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'1DS.fig'));
% 
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2DS.png'));
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2DS.eps'),'epsc');
%         saveas(figure(2),strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'2DS.fig'));
%         save(strcat(parentDir,'Compiled\Place Fields\2_PlaceFieldMeasures',tag,'DS.mat'),'stats');
%     end
end
    
end

