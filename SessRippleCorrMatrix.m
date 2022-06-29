function SessRippleCorrMatrix(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',true,@islogical);
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

col1 = cbrewer('seq','Blues',11);

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist(strcat('Summ\RippleCorrMatrix.mat'),'file') && ~force 
    disp('Ripple corr matrix already computed! Loading file.');
    load(strcat('Summ\RippleCorrMatrix.mat'));
else
    for rr = 1:3
        corrMatrix{rr} = [];
    end
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        
        %Load ripples
        load([sessionInfo.FileName '.ripples.events.mat']);
        load([sessionInfo.FileName '.pulses.events.mat']);
        
        %Load sleep states
        load([sessionInfo.FileName '.SleepState.states.mat']);
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);        
        
        load([sessionInfo.FileName '.SessionPulses.Events.mat']); 
        load([sessionInfo.FileName '.MergePoints.events.mat']); 
        
        %Find home cage intervals
        efields = fields(sessionPulses);
        ts_hc = [];
        for fn = 1:length(MergePoints.foldernames)
            flag = 1;
            for fx = 1:length(efields)
                if strcmp(efields{fx},MergePoints.foldernames{fn}) == 1 %If folder if a maze folder
                    flag = 0;
                end
            end
            if flag
                ts_hc = [ts_hc;MergePoints.timestamps(fn,:)];
            end
        end
        
        %Build spike matrix for sleep periods
        spikeMat = [];
        timeMat = [];   
        events = SleepState.ints.NREMstate;
        % Build matrix, concatenate sleep epochs
        for trials = 1:size(events,1)         
            spkMat = [];
            spkMat_ts = [];
            %check if interval lies within homecage
            status = InIntervals(events(trials,1),ts_hc);
            if status == 0
                continue
            end
            for unit = 1:length(spikes.times)
                if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1
                    continue
                else 
                   spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .05, 'win',events(trials,:));
                   spkMat = [spkMat spkData.data];                                        
                   spkMat_ts = [spkData.timestamps];
                end                               
            end
           spikeMat = [spikeMat;spkMat];
           timeMat = [timeMat;spkMat_ts];            
        end
                
        % Assign logicals for the timeMat 
        for rr = 1:3              
            if exist('pulses')
                if rr <= 2
                    pulTr = (pulses.stimComb==rr);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
                end
            end

            %Only select the pulses that happened in the home cage,
            %i.e., which were 5 seconds long
            homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
            homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
            pulTr = pulTr & homeCagePulseidx;
            events = pulses.intsPeriods(1,pulTr)';

            %Generate logicals for ripples in pre, prestim, post, poststim
            ripple_prestim = [];  
            ripple_pre = [];                
            ripple_post = [];  
            ripple_poststim = [];               

            ripple_pre(1:length(ripples.peaks)) = 0;                
            ripple_post(1:length(ripples.peaks)) = 0;                
            ripple_prestim(1:length(ripples.peaks)) = 0;                
            ripple_poststim(1:length(ripples.peaks)) = 0;   
            
            for pp = 1:length(ripples.peaks)
                tempDiff = ripples.peaks(pp) - events;

                if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                   [~,idxmin] =  min(abs(tempDiff));
                   if tempDiff(idxmin) > 0
                       ripple_post(pp) = 1;
                   elseif tempDiff(idxmin) < 0
                       ripple_pre(pp) = 1;
                   end
                elseif min(abs(tempDiff)) >5 && min(abs(tempDiff)) <=10% If a ripple occurs within 5 seconds of a stimulus
                   [~,idxmin] =  min(abs(tempDiff));
                   if tempDiff(idxmin) > 0
                       ripple_poststim(pp) = 1;
                   elseif tempDiff(idxmin) < 0
                       ripple_prestim(pp) = 1;
                   end                   
                else
                    continue
                end
            end      

            ripple_logical{rr} = [logical(ripple_prestim)' logical(ripple_pre)'...
                logical(ripple_post)' logical(ripple_poststim)'];        
        end
        
        % Assign logicals for the timeMat    
        timeLog(1,1:length(timeMat)) = 0;
        timeLog(2,1:length(timeMat)) = 0;
        timeLog(3,1:length(timeMat)) = 0;
        
%         for rr = 1:3              
%             if exist('pulses')
%                 if rr <= 2
%                     pulTr = (pulses.stimComb==rr);
%                 else
%                     pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
%                 end
%             end
% 
%             %Only select the pulses that happened in the home cage,
%             %i.e., which were 5 seconds long
%             homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
%             homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
%             pulTr = pulTr & homeCagePulseidx;
%             events = pulses.intsPeriods(1,pulTr)';
%             
%             for trials = 1:length(events)
%                 idx = find(timeMat>= events(trials) & timeMat< (events(trials)+5));
%                 if ~isempty(idx)
%                     timeLog(rr,idx) = 3;
%                 end
%                 idx = find(timeMat>= (events(trials)-5) & timeMat< events(trials));
%                 if ~isempty(idx)
%                     timeLog(rr,idx) = 2;
%                 end
%                 idx = find(timeMat>= (events(trials)-10) & timeMat< (events(trials)-5));
%                 if ~isempty(idx)
%                     timeLog(rr,idx) = 1;
%                 end
%                 idx = find(timeMat>= (events(trials)+5) & timeMat< (events(trials)+10));
%                 if ~isempty(idx)
%                     timeLog(rr,idx) = 4;
%                 end                
%             end
%         end        

        for rr = 1:3   
            for ripNum  = 1:4
               % peakTS = ripples.peaks(ripple_logical{rr}(:,ripNum));
                peakTS = ripples.timestamps(ripple_logical{rr}(:,ripNum),:);                
                for tt = 1:size(peakTS,1)
                    %Find bins within 500 ms of the ripple peak
                    idx = timeMat>=(peakTS(tt,1)-0.5) & timeMat<=(peakTS(tt,2)+0.5);
                    if sum(idx)>0
                        timeLog(rr,idx) = ripNum;
                    end                               
                end
            end
        end
       
       if isempty(spikeMat)
           continue;
       end    
       
       for rr = 1:3
           spkMat_prestim = spikeMat(timeLog(rr,:)==1,:);
           spkMat_pre = spikeMat(timeLog(rr,:)==2,:);
           spkMat_post = spikeMat(timeLog(rr,:)==3,:);
           spkMat_poststim = spikeMat(timeLog(rr,:)==4,:);
           
           celltokeep = (sum(spkMat_prestim,1)>0) & (sum(spkMat_pre,1)>0) & (sum(spkMat_post,1)>0) & (sum(spkMat_poststim,1)>0);
            
           if sum(celltokeep) < 5 
               continue;
           end
           
           corr_prestim = corr(spkMat_prestim(:,celltokeep),'Type','Spearman','Rows','pairwise');
           corr_pre = corr(spkMat_pre(:,celltokeep),'Type','Spearman','Rows','pairwise');
           corr_post = corr(spkMat_post(:,celltokeep),'Type','Spearman','Rows','pairwise');
           corr_poststim = corr(spkMat_poststim (:,celltokeep),'Type','Spearman','Rows','pairwise');
              
           Z = linkage(corr_pre,'weighted');
           T = cluster(Z,'Maxclust',5);
           [~,idxsort] = sort(T);
           corr_sorted_pre = corr_pre(idxsort,idxsort);
           corr_sorted_pre(corr_sorted_pre==0) = 0.0000000001; % change this value because it gives errors later
           corr_prestim_pre = corr_prestim(idxsort,idxsort);
           corr_prestim_pre(corr_prestim_pre==0) = 0.0000000001; % change this value because it gives errors later
           corr_post_pre = corr_post(idxsort,idxsort);
           corr_post_pre(corr_post_pre==0) = 0.0000000001; % change this value because it gives errors later
           corr_poststim_pre = corr_poststim(idxsort,idxsort);
           corr_poststim_pre(corr_poststim_pre==0) = 0.0000000001; % change this value because it gives errors later
           
           
           Z = linkage(corr_post,'weighted');
           T = cluster(Z,'Maxclust',5);
           [~,idxsortstim] = sort(T);
           corr_sorted_post = corr_post(idxsortstim,idxsortstim);
           corr_sorted_post(corr_sorted_post==0) = 0.0000000001; % change this value because it gives errors later
           corr_prestim_post = corr_prestim(idxsortstim,idxsortstim);
           corr_prestim_post(corr_prestim_post==0) = 0.0000000001; % change this value because it gives errors later
           corr_pre_post = corr_pre(idxsortstim,idxsortstim);
           corr_pre_post(corr_pre_post==0) = 0.0000000001; % change this value because it gives errors later
           corr_poststim_post = corr_poststim(idxsortstim,idxsortstim);
           corr_poststim_post(corr_poststim_post==0) = 0.0000000001; % change this value because it gives errors later
                 
           basecorr = corrcoef(corr_sorted_pre(tril(corr_sorted_pre,-1)~=0),corr_prestim_pre(tril(corr_prestim_pre,-1)~=0),'rows','complete');             
           poststimcorr = corrcoef(corr_sorted_pre(tril(corr_sorted_pre,-1)~=0),corr_poststim_pre(tril(corr_poststim_pre,-1)~=0),'rows','complete'); 
           bscorr = corrcoef(corr_sorted_pre(tril(corr_sorted_pre,-1)~=0),corr_post_pre(tril(corr_post_pre,-1)~=0),'rows','complete'); 
            
           S_basecorr = corrcoef(corr_sorted_post(tril(corr_sorted_post,-1)~=0),corr_prestim_post(tril(corr_prestim_post,-1)~=0),'rows','complete');             
           S_poststimcorr = corrcoef(corr_sorted_post(tril(corr_sorted_post,-1)~=0),corr_poststim_post(tril(corr_poststim_post,-1)~=0),'rows','complete'); 
           S_bscorr = corrcoef(corr_sorted_post(tril(corr_sorted_post,-1)~=0),corr_pre_post(tril(corr_pre_post,-1)~=0),'rows','complete'); 
           
           
           corrMatrix{rr} = [corrMatrix{rr}; basecorr(1,2) bscorr(1,2) poststimcorr(1,2) S_basecorr(1,2) S_bscorr(1,2) S_poststimcorr(1,2)];    
            
            if makePlot
                figure
                subplot(2,4,1)
                imagesc(corr_prestim_pre)
                caxis([0 0.3])

                subplot(2,4,2)
                imagesc(corr_sorted_pre)
                caxis([0 0.3])
                
                subplot(2,4,3)
                imagesc(corr_post_pre)
                caxis([0 0.3])
                
                subplot(2,4,4)
                imagesc(corr_poststim_pre)
                caxis([0 0.3])           
                
                subplot(2,4,5)
                imagesc(corr_prestim_post)
                caxis([0 0.3])
                
                subplot(2,4,6)
                imagesc(corr_pre_post)
                caxis([0 0.3])   

                subplot(2,4,7)
                imagesc(corr_sorted_post)
                caxis([0 0.3])
                
                subplot(2,4,8)
                imagesc(corr_poststim_post)
                caxis([0 0.3])   
                
                colormap(col1)
                
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCorr',num2str(rr),'.fig'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCorr',num2str(rr),'.png'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCorr',num2str(rr),'.eps'),'epsc')      
                close all
            end            
       end    
       clear timeLog ripple_logical
    end
    
    if saveMat
        save([expPath '\Summ\' 'RippleCorrMatrix.mat'], 'corrMatrix','-v7.3');
    end    
end
    
end

