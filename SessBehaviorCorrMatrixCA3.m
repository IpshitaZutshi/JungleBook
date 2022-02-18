function SessBehaviorCorrMatrixCA3(varargin)

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

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

col1 = cbrewer('seq','Blues',11);

if exist(strcat('Summ\CorrMatrix.mat'),'file') && ~force 
    disp('Corr matrix already computed! Loading file.');
    load(strcat('Summ\CorrMatrix.mat'));
else
    
    CorrMatrix = [];

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
        
        efields = fieldnames(sessionPulses);            
        stimMat = []; 
        spikeMat = [];
        
        for jj = 1:length(efields)
                        
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps; 
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            stim = sessionPulses.(efields{jj}).stim;
            events = [endDelay(1:length(rewardTS)-1) rewardTS(2:end)];
            
  
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
               if jj == 1
                   stimMat = [stimMat;ones(size(spkMat,1),1)*stim(trials)];
               elseif jj ==2
                   stimMat = [stimMat;ones(size(spkMat,1),1)*(stim(trials)+2)];
               end
            end  
        end

            
        baseidx = find(stimMat==0);
        SpikeCountFirstBase = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
        SpikeCountSecondBase = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);        

        baseidx = find(stimMat==1);
        SpikeCountFirstmEC = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
        SpikeCountSecondmEC = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);             
        
        baseidx = find(stimMat==2);
        SpikeCountFirstCA3 = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
        SpikeCountSecondCA3 = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);                     
        
        baseidx = find(stimMat==3);
        SpikeCountFirstBoth = spikeMat(baseidx(1:floor(length(baseidx)/2)),:);
        SpikeCountSecondBoth = spikeMat(baseidx(floor(length(baseidx)/2)+1:end),:);               

        celltokeep = (sum(SpikeCountFirstBase,1)>0) & (sum(SpikeCountSecondBase,1)>0) & (sum(SpikeCountFirstmEC,1)>0) & (sum(SpikeCountSecondmEC,1)>0) &...
            (sum(SpikeCountFirstCA3,1)>0) & (sum(SpikeCountSecondCA3,1)>0) & (sum(SpikeCountFirstBoth,1)>0) & (sum(SpikeCountSecondBoth,1)>0);

        if sum(celltokeep) < 5 
            continue;
        end
        
        corrfirstbase = corr(SpikeCountFirstBase(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrsecondbase = corr(SpikeCountSecondBase(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrfirstmec = corr(SpikeCountFirstmEC(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrsecondmec = corr(SpikeCountSecondmEC(:,celltokeep),'Type','Spearman','Rows','pairwise');        
        corrfirstca3 = corr(SpikeCountFirstCA3(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrsecondca3 = corr(SpikeCountSecondCA3(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrfirstboth = corr(SpikeCountFirstBoth(:,celltokeep),'Type','Spearman','Rows','pairwise');
        corrsecondboth = corr(SpikeCountSecondBoth(:,celltokeep),'Type','Spearman','Rows','pairwise');
        
        Z = linkage(corrfirstbase,'weighted');
        T = cluster(Z,'Maxclust',5);
        [~,idxsort] = sort(T);
        corr_sorted_first = corrfirstbase(idxsort,idxsort);
        corr_sorted_first(corr_sorted_first==0) = 0.0000000001; % change this value because it gives errors later
        corr_sorted_second = corrsecondbase(idxsort,idxsort);
        corr_sorted_second(corr_sorted_second==0) = 0.0000000001; % change this value because it gives errors later
        
        Z = linkage(corrfirstmec,'weighted');
        T = cluster(Z,'Maxclust',5);
        [~,idxsortmec] = sort(T);
        corr_sorted_mecfirst = corrfirstmec(idxsortmec,idxsortmec);
        corr_sorted_mecfirst(corr_sorted_mecfirst==0) = 0.0000000001; % change this value because it gives errors later
        corr_sorted_mecsecond = corrsecondmec(idxsortmec,idxsortmec);
        corr_sorted_mecsecond(corr_sorted_mecsecond==0) = 0.0000000001; % change this value because it gives errors later        

        Z = linkage(corrfirstca3,'weighted');
        T = cluster(Z,'Maxclust',5);
        [~,idxsortca3] = sort(T);
        corr_sorted_ca3first = corrfirstca3(idxsortca3,idxsortca3);      
        corr_sorted_ca3first(corr_sorted_ca3first==0) = 0.0000000001; % change this value because it gives errors later
        corr_sorted_ca3second = corrsecondca3(idxsortca3,idxsortca3);
        corr_sorted_ca3second(corr_sorted_ca3second==0) = 0.0000000001; % change this value because it gives errors later       
        
        Z = linkage(corrfirstboth,'weighted');
        T = cluster(Z,'Maxclust',5);
        [~,idxsortboth] = sort(T);
        corr_sorted_bothfirst = corrfirstboth(idxsortboth,idxsortboth);      
        corr_sorted_bothfirst(corr_sorted_bothfirst==0) = 0.0000000001; % change this value because it gives errors later
        corr_sorted_bothsecond = corrsecondboth(idxsortboth,idxsortboth);
        corr_sorted_bothsecond(corr_sorted_bothsecond==0) = 0.0000000001; % change this value because it gives errors later            

        corr_mecbase = corrfirstmec(idxsort,idxsort);
        corr_mecbase(corr_mecbase==0) = 0.0000000001; % change this value because it gives errors later  
        
        corr_ca3base = corrfirstca3(idxsort,idxsort);
        corr_ca3base(corr_ca3base==0) = 0.0000000001; % change this value because it gives errors later  

        corr_bothbase = corrfirstboth(idxsort,idxsort);
        corr_bothbase(corr_bothbase==0) = 0.0000000001; % change this value because it gives errors later  

        corr_ca3both = corrfirstboth(idxsortca3,idxsortca3);
        corr_ca3both(corr_ca3both==0) = 0.0000000001; % change this value because it gives errors later  
        
        basecorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_sorted_second(tril(corr_sorted_second,-1)~=0),'rows','complete');             
        meccorr = corrcoef(corr_sorted_mecfirst(tril(corr_sorted_mecfirst,-1)~=0),corr_sorted_mecsecond(tril(corr_sorted_mecsecond,-1)~=0),'rows','complete');             
        ca3corr = corrcoef(corr_sorted_ca3first(tril(corr_sorted_ca3first,-1)~=0),corr_sorted_ca3second(tril(corr_sorted_ca3second,-1)~=0),'rows','complete'); 
        bothcorr = corrcoef(corr_sorted_bothfirst(tril(corr_sorted_bothfirst,-1)~=0),corr_sorted_bothsecond(tril(corr_sorted_bothsecond,-1)~=0),'rows','complete'); 
        mecbasecorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_mecbase(tril(corr_mecbase,-1)~=0),'rows','complete'); 
        ca3basecorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_ca3base(tril(corr_ca3base,-1)~=0),'rows','complete'); 
        bothbasecorr = corrcoef(corr_sorted_first(tril(corr_sorted_first,-1)~=0),corr_bothbase(tril(corr_bothbase,-1)~=0),'rows','complete'); 
        ca3bothcorr = corrcoef(corr_sorted_ca3first(tril(corr_sorted_ca3first,-1)~=0),corr_ca3both(tril(corr_ca3both,-1)~=0),'rows','complete'); 

        CorrMatrix = [CorrMatrix; basecorr(1,2) mecbasecorr(1,2) meccorr(1,2) ca3basecorr(1,2) ca3corr(1,2) bothbasecorr(1,2) bothcorr(1,2) ca3bothcorr(1,2)];    

        if makePlot           
            figure
            subplot(6,2,1)
            imagesc(corr_sorted_first)
            caxis([0 0.3])

            subplot(6,2,2)
            imagesc(corr_sorted_second)
            caxis([0 0.3])
            
            subplot(6,2,3)
            imagesc(corr_sorted_mecfirst)
            caxis([0 0.3])

            subplot(6,2,4)
            imagesc(corr_sorted_mecsecond)
            caxis([0 0.3])                 

            subplot(6,2,5)
            imagesc(corr_sorted_ca3first)
            caxis([0 0.3])

            subplot(6,2,6)
            imagesc(corr_sorted_ca3second)
            caxis([0 0.3])           

            subplot(6,2,7)            
            imagesc(corr_sorted_bothfirst)
            caxis([0 0.3])
            
            subplot(6,2,8)            
            imagesc(corr_sorted_bothsecond)
            caxis([0 0.3])            
            
            subplot(6,2,9)
            imagesc(corr_mecbase)
            caxis([0 0.3])   
            
            subplot(6,2,10)
            imagesc(corr_ca3base)
            caxis([0 0.3])   

            subplot(6,2,11)
            imagesc(corr_bothbase)
            caxis([0 0.3])   
            
            subplot(6,2,12)
            imagesc(corr_ca3both)
            caxis([0 0.3])   
            colorbar
            
            colormap(col1)
            saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.fig'))
            saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.png'))
            saveas(gcf,strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\Analysis\SessCorr.eps'),'epsc')      
            close all
        end            
    end    

    if saveMat
        save([expPath '\Summ\' 'CorrMatrix.mat'], 'CorrMatrix','-v7.3');
    end
end

