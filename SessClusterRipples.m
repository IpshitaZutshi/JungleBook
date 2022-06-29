function SessClusterRipples(varargin)

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

if exist(strcat('Summ\RippleClusters.mat'),'file') && ~force 
    disp('Ripple clusters already computed! Loading file.');
    load(strcat('Summ\RippleClusters.mat'));
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

        %% calculate ripple content
        for rr = 1:3              
            rippleMat{rr} = [];
            rippleID{rr} = [];

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

            %Generate logicals for ripples in pre versus post
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

            ripple_logical = [logical(ripple_prestim)' logical(ripple_pre)'...
                logical(ripple_post)' logical(ripple_poststim)'];     

             %% Now calculate the cells active for these specific ripples
            for rippLog = 1:size(ripple_logical,2)
                ripX =[];
                ripLog = [];
                if sum(ripple_logical(:,rippLog))~=0
                    int = ripples.peaks(ripple_logical(:,rippLog));
                    for jj = 1:size(spikes.UID,2)
                        % Get psth
                        if strcmp(cell_metrics.brainRegion{jj},'CA1')~=1 || strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell')~=1
                            continue
                        else
                            [status, interval, index] = InIntervals(spikes.times{jj},[int-0.5 int+0.5]);
                        end
                        a = zeros(1,length(int));
                        a(interval(status>0)) = 1;   
                        ripX = [ripX; a];
                     end
                end
                ripLog(1:size(ripX,2)) = rippLog;
                rippleMat{rr} = [rippleMat{rr} ripX]; 
                rippleID{rr} = [rippleID{rr} ripLog];
            end   
            
           spkMat_prestim = rippleMat{rr}(:,rippleID{rr}==1)';
           spkMat_pre = rippleMat{rr}(:,rippleID{rr}==2)';
           spkMat_post = rippleMat{rr}(:,rippleID{rr}==3)';
           spkMat_poststim = rippleMat{rr}(:,rippleID{rr}==4)';
           
           if isempty(spkMat_prestim) || isempty(spkMat_pre) ||...
                   isempty(spkMat_post) || isempty(spkMat_poststim)
               continue;
           end 
           
           celltokeep = (sum(spkMat_prestim,1)>0) & (sum(spkMat_pre,1)>0) & (sum(spkMat_post,1)>0) & (sum(spkMat_poststim,1)>0);
            
           if sum(celltokeep) < 5 
               continue;
           end
           
           corr_prestim = corr(spkMat_prestim,'Type','Spearman','Rows','pairwise');
           corr_pre = corr(spkMat_pre,'Type','Spearman','Rows','pairwise');
           corr_post = corr(spkMat_post,'Type','Spearman','Rows','pairwise');
           corr_poststim = corr(spkMat_poststim,'Type','Spearman','Rows','pairwise');
              
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
                
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCluster',num2str(rr),'.fig'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCluster',num2str(rr),'.png'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\SummaryFigures\RippleCluster',num2str(rr),'.eps'),'epsc')      
                close all
            end               
            clear ripple_logical
        end
    end
    
    if saveMat
        save([expPath '\Summ\' 'RippleClusters.mat'], 'corrMatrix','-v7.3');
    end    
end
    
end

