function Summary = sessionRemappingTrajectorynoToneVer2

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
    'IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
plotfig = 0;

YlGnBu=cbrewer('seq', 'YlGnBu', 11);

AllCellsCorr = [];
AllCellsSimPre = [];
AllCellsSim = [];

AllCellsRateMapLin1 = [];
AllCellsRateMapLin2 = [];
AllCellsRateMapTone1 = [];
AllCellsRateMapTone2 = [];
AllCellsRateMapLinEnd1 = [];
AllCellsRateMapLinEnd2 = [];

AllCellsRateMapLinFirstHalf = [];
AllCellsRateMapLinSecondHalf = [];      
AllCellsRateMapToneFirstHalf = [];
AllCellsRateMapToneSecondHalf = [];

for ii = 1:length(sess)
    
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.rateMapsTrial.cellinfo.mat']);
    firingMapsTrial = load(file(1).name);    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);    
    
    % Cluster trajectories and extract similar trials
    [trialNum, trialLabel, trajectLabel] = clusterTrajectories('plotfig',false);
    
    % Check if the post-tone trials have more similarity to the no tone or
    % tone    
    typePre = mode(trajectLabel(trialLabel==1));
    typeTone = mode(trajectLabel(trialLabel==2));
    typePost = mode(trajectLabel(trialLabel==3));

    if typePost == typeTone
        sim = 1;
    else
        sim = 0 ;
    end
    
    if typePre == typeTone
        simPre  = 1;
    else
        simPre  = 0;
    end
    
    for kk=1:length(cell_metrics.UID)
        for zz = 1:10
            Maps{zz} = [];
        end
        avgMaps = [];
        
        %% First, only extract cells with a field in either of the three
        %trial types
        linMapInit = firingMaps.forward.rateMaps{kk}{1};
        linMapEnd = firingMaps.forward.rateMaps{kk}{28};
        toneMap = firingMaps.forward.rateMaps{kk}{7};
        
        %% Detect fields
        Field1 = detectFields(linMapInit);
        Field2 = detectFields(toneMap);
        Field3 = detectFields(linMapEnd);

        simil  = [];
        similPre = [];
        
        if isempty(Field1) && isempty(Field2)% && isempty(Field3)
            continue
        end

        %% Now extract firing map of each category according to trajectories        
        for zz = 1:6
            switch zz
                case 1
                    idx = (trialLabel == 1) & trajectLabel == 1;
                case 2
                    idx = (trialLabel == 1) & trajectLabel == 2;                     
                case 3
                    idx = trialLabel == 2 & trajectLabel == 1;
                case 4
                    idx = trialLabel == 2 & trajectLabel == 2;
                case 5
                    idx = trialLabel == 3 & trajectLabel == 1;                     
                case 6
                    idx = trialLabel == 3 & trajectLabel == 2;                    
                                 
            end
            
            trialIdx  = trialNum(idx);
            for tt = 1:length(trialIdx)
                Maps{zz} = [Maps{zz};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx(tt)}];
            end
            
            if isempty(Maps{zz})
                avgMaps(zz,1:50) = nan;
            else
                avgMaps(zz,:) = nanmean(Maps{zz},1);
            end
            
        end
        
        % Also take subsets of within stability
        idx = find(trialLabel == 1 & trajectLabel == typePre);
        trialIdx1 = trialNum(idx(1:floor(length(idx)/2)));
        trialIdx2 = trialNum(idx(floor(length(idx)/2)+1:end));
        
        for tt = 1:length(trialIdx1)
            Maps{7} = [Maps{7};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx1(tt)}];
        end
        for tt = 1:length(trialIdx2)
            Maps{8} = [Maps{8};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx2(tt)}];
        end
            
        if isempty(Field1)
            avgMaps(7,1:50) = nan;
            avgMaps(8,1:50) = nan;
        else
            avgMaps(7,:) = nanmean(Maps{7},1);
            avgMaps(8,:) = nanmean(Maps{8},1);
        end
        
        % Also take subsets of within tone stability
        idx = find(trialLabel == 2 & trajectLabel == typeTone);
        trialIdx1 = trialNum(idx(1:floor(length(idx)/2)));
        trialIdx2 = trialNum(idx(floor(length(idx)/2)+1:end));
        
        for tt = 1:length(trialIdx1)
            Maps{9} = [Maps{9};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx1(tt)}];
        end
        for tt = 1:length(trialIdx2)
            Maps{10} = [Maps{10};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx2(tt)}];
        end
            
        if isempty(Field2)
            avgMaps(9,1:50) = nan;
            avgMaps(10,1:50) = nan;
        else
            avgMaps(9,:) = nanmean(Maps{9},1);
            avgMaps(10,:) = nanmean(Maps{10},1);
        end
        
        AllCellsRateMapLin1 = [AllCellsRateMapLin1;avgMaps(1,:)];
        AllCellsRateMapLin2 = [AllCellsRateMapLin2;avgMaps(2,:)];
        AllCellsRateMapTone1 = [AllCellsRateMapTone1;avgMaps(3,:)];
        AllCellsRateMapTone2 = [AllCellsRateMapTone2;avgMaps(4,:)];
        AllCellsRateMapLinFirstHalf = [AllCellsRateMapLinFirstHalf;avgMaps(7,:)];
        AllCellsRateMapLinSecondHalf = [AllCellsRateMapLinSecondHalf;avgMaps(8,:)];
        AllCellsRateMapToneFirstHalf = [AllCellsRateMapToneFirstHalf;avgMaps(9,:)];
        AllCellsRateMapToneSecondHalf = [AllCellsRateMapToneSecondHalf;avgMaps(10,:)];
        
        corrMatrix = corrcoef(avgMaps','Rows','pairwise');
        C = tril(corrMatrix,-1);
        C = C(isnan(C) | C~=0);

        simil = [simil sim];
        similPre = [similPre simPre];
        
        if plotfig
            saveLoc = strcat(expPath,sess{ii},'\TrajectMaps');
            if ~isfolder('TrajectMaps')
                mkdir('TrajectMaps')
            end  
            
            figure
            subplot(6,2,1)
            imagesc(linMapInit)
            colormap(YlGnBu)
            title(strcat('Max:',num2str(max(linMapInit)),' Hz'))
            
            subplot(6,2,5)
            imagesc(toneMap)
            title(strcat('Max:',num2str(max(toneMap)),' Hz'))
            
            subplot(6,2,9)
            imagesc(linMapEnd)
            title(strcat('Max:',num2str(max(linMapEnd)),' Hz'))
            
            for zz = 1:6
                subplot(6,2,2*(zz-1)+2)
                imagesc(Maps{zz})
                title(strcat('Max:',num2str(max(max(Maps{zz}))),' Hz'))
            end
            
            subplot(6,2,11)
            imagesc(avgMaps)
            title('Average trajectory maps')
            
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.png'],'png');
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.fig'],'fig');
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.eps'],'epsc');
            close all
        end
            
        AllCellsCorr = [AllCellsCorr;C'];
        AllCellsSim = [AllCellsSim;simil];
        AllCellsSimPre = [AllCellsSimPre;similPre];
    end
end

a = corr(AllCellsRateMapLinFirstHalf',AllCellsRateMapLinSecondHalf','rows','pairwise');
b = diag(a);

a1 = corr(AllCellsRateMapLin1',AllCellsRateMapLin2','rows','pairwise');
b1 = diag(a1);

a2 = corr(AllCellsRateMapLin1',AllCellsRateMapTone1','rows','pairwise');
b2 = diag(a2);

a3 = corr(AllCellsRateMapLin1',AllCellsRateMapTone2','rows','pairwise');
b3 = diag(a3);

a4 = corr(AllCellsRateMapLin2',AllCellsRateMapTone2','rows','pairwise');
b4 = diag(a4);

a5 = corr(AllCellsRateMapLin2',AllCellsRateMapTone1','rows','pairwise');
b5 = diag(a5);

a6 = corr(AllCellsRateMapTone1',AllCellsRateMapTone2','rows','pairwise');
b6 = diag(a6);

a7 = corr(AllCellsRateMapToneFirstHalf',AllCellsRateMapToneSecondHalf','rows','pairwise');
b7 = diag(a7);

Summary.comp1 = [b;b7];
Summary.comp2 = [b1; b6];
Summary.comp3 = [b2; b4];
Summary.comp4 = [b3; b5];

end

