function SessBehaviorCorrMatrix(varargin)

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

if exist(strcat('Summ\CorrMatrix.mat'),'file') && ~force 
    disp('Corr matrix already computed! Loading file.');
    load(strcat('Summ\CorrMatrix.mat'));
else
    for rr = 1:3
        for cc = 1:2
            CorrMatrix{rr,cc} = [];    
        end
    end
    
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
        
        for jj = 1:length(efields)
                        
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps; 
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            stim = sessionPulses.(efields{jj}).stim;
            events = [endDelay(1:length(rewardTS)-1) rewardTS(2:end)];
            
            spikeMat = [];
            stimMat = [];   
            % Build matrix, concatenate trials
            for trials = 1:size(events,1)         
               spkMat = [];
               %spikes
               for unit = 1:length(spikes.times)
                   if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1 || strcmp(cell_metrics.putativeCellType(unit),'Pyramidal Cell')~=1
                       continue
                   else
                       spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .125, 'win',events(trials,:));
                       spkMat = [spkMat spkData.data];
                   end
               end
               spikeMat = [spikeMat;spkMat];
               stimMat = [stimMat;ones(size(spkMat,1),1)*stim(trials)];
            end  
            
            baseidx = find(stimMat==0);
            SpikeCountFirstBase = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
            SpikeCountSecondBase = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);        
            
            baseidx = find(stimMat==1);
            SpikeCountFirstStim = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
            SpikeCountSecondStim = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);                     
            
            celltokeep = (sum(SpikeCountFirstBase,1)>0) & (sum(SpikeCountSecondBase,1)>0) & (sum(SpikeCountFirstStim,1)>0) & (sum(SpikeCountSecondStim,1)>0);
            
            if sum(celltokeep) < 5 
                continue;
            end            
            corrfirstbase = corr(SpikeCountFirstBase(:,celltokeep),'Type','Spearman','Rows','pairwise');
            corrsecondbase = corr(SpikeCountSecondBase(:,celltokeep),'Type','Spearman','Rows','pairwise');
            corrfirststim = corr(SpikeCountFirstStim(:,celltokeep),'Type','Spearman','Rows','pairwise');
            corrsecondstim = corr(SpikeCountSecondStim(:,celltokeep),'Type','Spearman','Rows','pairwise');
              
            Z = linkage(corrfirstbase,'weighted');
            T = cluster(Z,'Maxclust',5);
            [~,idxsort] = sort(T);
            corr_sorted_first = corrfirstbase(idxsort,idxsort);
            corr_sorted_first(corr_sorted_first==0) = 0.0000000001; % change this value because it gives errors later
            corr_sorted_second = corrsecondbase(idxsort,idxsort);
            corr_sorted_second(corr_sorted_second==0) = 0.0000000001; % change this value because it gives errors later

            Z = linkage(corrfirststim,'weighted');
            T = cluster(Z,'Maxclust',5);
            [~,idxsortstim] = sort(T);
            corr_sorted_stimfirst = corrfirststim(idxsortstim,idxsortstim);      
            corr_sorted_stimfirst(corr_sorted_stimfirst==0) = 0.0000000001; % change this value because it gives errors later
            corr_sorted_stimsecond = corrsecondstim(idxsortstim,idxsortstim);
            corr_sorted_stimsecond(corr_sorted_stimsecond==0) = 0.0000000001; % change this value because it gives errors later            
            corr_stimbase = corrfirststim(idxsort,idxsort);
            corr_stimbase(corr_stimbase==0) = 0.0000000001; % change this value because it gives errors later  
            
            basecorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_sorted_second(tril(corr_sorted_second,-1)~=0),'rows','complete');             
            stimcorr = corrcoef(corr_sorted_stimfirst(tril(corr_sorted_stimfirst,-1)~=0),corr_sorted_stimsecond(tril(corr_sorted_stimsecond,-1)~=0),'rows','complete'); 
            bscorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_stimbase(tril(corr_stimbase,-1)~=0),'rows','complete'); 
            
            CorrMatrix{region,target} = [CorrMatrix{region,target}; basecorr(1,2) bscorr(1,2) stimcorr(1,2)];    
            
            if makePlot
                figure
                subplot(3,2,1)
                imagesc(corr_sorted_first)
                caxis([0 0.3])

                subplot(3,2,2)
                imagesc(corr_sorted_second)
                caxis([0 0.3])
                
                subplot(3,2,3)
                imagesc(corr_sorted_stimfirst)
                caxis([0 0.3])
                
                subplot(3,2,4)
                imagesc(corr_sorted_stimsecond)
                caxis([0 0.3])           
                
                subplot(3,2,5)
                imagesc(corrfirstbase(idxsortstim,idxsortstim))
                caxis([0 0.3])
                
                subplot(3,2,6)
                imagesc(corrfirststim(idxsort,idxsort))
                caxis([0 0.3])   
                
                colormap(col1)
                
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.fig'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.png'))
                saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.eps'),'epsc')      
                close all
            end            
        end    
    end

    if saveMat
        save([expPath '\Summ\' 'CorrMatrix.mat'], 'CorrMatrix','-v7.3');
    end
end

