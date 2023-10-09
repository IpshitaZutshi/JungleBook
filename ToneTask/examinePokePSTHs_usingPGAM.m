function Summary = examinePokePSTHs_usingPGAM(varargin)

p = inputParser;
addParameter(p,'plotfig',false,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;

%% Look at the following types of pokes - 
%   1) Forward run - no tone trials, end, 
%   2) Forward run - tone trials, 
%   3) Spontaneous middle ports,facing forward
%   4) Spontaneous middle ports, running back
%   5) Home port

% Next, separate correct and incorrect for cells tuned to type 2, and type 5

% Finally, for the ones tuned to 2, separate the 6 ports

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24','IZ47\Final\IZ47_230710_sess25',...
    'IZ47\Final\IZ47_230712_sess27'};  

PGAMpath = 'Z:\Buzsakilabspace\LabShare\AthinaApostolelli\PGAM\postprocess';
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

for tt = 1:15
    Summary.psthReward.lickTypes{tt} = [];
end

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);    
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);       
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*.DigitalIn.Events.mat']);
    load(file(1).name);
    
    fileloc = strcat(PGAMpath,'\',sess{ii}(1:4),'\',sess{ii}(12:end),'\',sess{ii}(12:end),'.postprocessPGAM.mat');
    load(fileloc)
     
    lickCellIDs = resultsPGAM.tuned_licksP;
    
    licktimes = [];
    lickport1 = [];
    lickport = [2 3 4 5 6 7];
    for ll = 1:5
        licktimes = [licktimes; digitalIn.timestampsOn{lickport(ll)+2}]; 
        lickport1 = [lickport1; ones(length(digitalIn.timestampsOn{lickport(ll)+2}),1)*lickport(ll)];
    end    
    
    %Sort licktimes
    [licktimesSorted,lickIdx] = sort(licktimes,'ascend');
    lickportSorted = lickport1(lickIdx);
    idxNew = [1;find(diff(lickportSorted)~=0)+1];
    
    tsLick2 = licktimesSorted(idxNew);
    
    %% Calculate PSTHS for each type of lick

    for tt = 1:15
        
        st = [];
        if tt ==1 %Forward run - no tone trials, end, 
            st = (behavTrials.timestamps(behavTrials.linTrial==1,2)-0.03);
        elseif tt == 2 %Forward run - tone trials, end, 
            st = (behavTrials.timestamps(behavTrials.linTrial==0,2)-0.03);                    
        elseif tt == 3 %Spontaneous licks,facing forward
            %No tone trials, spontaneous forward
            ts1 = behavTrials.timestamps(behavTrials.linTrial==1,1);
            ts2 = behavTrials.timestamps(behavTrials.linTrial==1,2)-0.03;
            intPer = InIntervals(tsLick2,[ts1 ts2]);
            st1 = tsLick2(intPer,1);   

            %Tone trials, spontaneous
            idxLin = behavTrials.linTrial==0;
            tsFwd = behavTrials.timestamps(idxLin,1);
            tsRet = behavTrials.timestamps(idxLin,2);
            ints = [tsRet(1:(end-1)) tsFwd(2:end)];
            intPer = InIntervals(tsLick2(:,1),ints);  
            st2 = tsLick2(intPer,1);   

            st = [st1; st2];

            %only take the ones where licks are in the forward
            %direction
            idx = [];
            v = [];
            for ss = 1:length(st)
                [~,idx(ss)] = min(abs(tracking.timestamps-st(ss)));
                v(ss) = tracking.position.vy(idx(ss));
            end                       

            if isempty(st)
                continue
            else
                st = st(v>1);  
            end           

        elseif tt == 4 %Spontaneous licks,facing backward
            %No tone trials, spontaneous return
            idxLin = behavTrials.linTrial==1;
            tsFwd = behavTrials.timestamps(idxLin,1);
            tsRet = behavTrials.timestamps(idxLin,2);
            ints = [tsRet(1:(end-1)) tsFwd(2:end)];
            intPer = InIntervals(tsLick2(:,1),ints);  
            st1 = tsLick2(intPer,1);   

            %Tone trials, spontaneous
            idxLin = behavTrials.linTrial==0;
            tsFwd = behavTrials.timestamps(idxLin,1);
            tsRet = behavTrials.timestamps(idxLin,2);
            ints = [tsRet(1:(end-1)) tsFwd(2:end)];
            intPer = InIntervals(tsLick2(:,1),ints);  
            st2 = tsLick2(intPer,1);   

            st = [st1; st2];

            %only take the ones where licks are in the reverse
            %direction
            idx = [];
            v = [];
            for ss = 1:length(st)
                [~,idx(ss)] = min(abs(tracking.timestamps-st(ss)));
                v(ss) = tracking.position.vy(idx(ss));
            end                       

            if ~isempty(st)                            
                st = st(v<-1);  
            end           

        elseif tt == 5 % home port
            st = (behavTrials.timestamps(behavTrials.linTrial==0,1)-0.03);
        elseif tt == 6 %Forward correct
            st = behavTrials.timestamps(behavTrials.linTrial ==0 & behavTrials.correct ==1,2)-0.03;
        elseif tt ==7 %Forward incorrect
            st = behavTrials.timestamps(behavTrials.linTrial ==0 & behavTrials.correct ==0,2)-0.03;
        elseif tt ==8 % return correct
            idx = find(behavTrials.linTrial(1:(end-1)) ==0 & behavTrials.correct(1:(end-1)) ==1);
            homestartidx = idx+1;
            st = behavTrials.timestamps(homestartidx,1)-0.03;
        elseif tt ==9 % return incorrect     
            idx = find(behavTrials.linTrial(1:(end-1)) ==0 & behavTrials.correct(1:(end-1)) ==0);
            homestartidx = idx+1;
            st = behavTrials.timestamps(homestartidx,1)-0.03;
        elseif tt >=10
            idx = behavTrials.lickLoc==(tt-10) & behavTrials.linTrial ==0;
            st = behavTrials.timestamps(idx,2)-0.03;
        end
            
        for kk=1:length(lickCellIDs)        
            if ~isempty(st)
                [stccg, tPSTH] = CCG({spikes.times{lickCellIDs(kk)} st},[],'binSize',0.1,'duration',4,'norm','rate');                
                Summary.psthReward.lickTypes{tt} = [Summary.psthReward.lickTypes{tt}; stccg(:,2,1)'];               
            else
                fillArr(1,1:41) = nan;
                Summary.psthReward.lickTypes{tt} = [Summary.psthReward.lickTypes{tt}; fillArr];               
            end
        end
    end
end

save('C:\Users\ipshi\Desktop\Summary.mat', 'Summary'); 

if plotfig

    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')

    YlGnBu=cbrewer('seq', 'YlGnBu', 11);
    colormap(YlGnBu)

    idxT = tPSTH<0 & tPSTH>=-0.3;
    avgRate = nanmean(psthReward{1}(:,idxT),2);
    for tt = 1:3
        newpsth{tt} = psthReward{tt}(avgRate>1,:);
    end
    idxT = tPSTH<0 & tPSTH>=-0.5;
    [~,idxMax2] = max(newpsth{1}(:,idxT),[],2);
    [~,idxMax] = sort(idxMax2);
    for tt = 1:3
        subplot(2,3,tt)
        %temp = zscore(newpsth{tt},[],2);
        %h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
        %set(h, 'AlphaData', ~isnan(temp))
        imagesc(tPSTH, 1:size(newpsth{tt},1),newpsth{tt});
        caxis([0 10])
        %ylim([125 350])
        colorbar

        subplot(2,3,tt+3)
        col = [0.5 0.5 0.5];
        meanpsth = nanmean(newpsth{tt},1);
        stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
        hold on
        fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
        hold on
        hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);
        line([0 0],[1 11],'Color','r')
        ylim([1 11])    
    end

    saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.png'));
    saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.eps'),'epsc');
    saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.fig'));
end
end

