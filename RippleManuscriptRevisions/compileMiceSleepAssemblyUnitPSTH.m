function compileMiceSleepAssemblyUnitPSTH(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'binSize',0.1,@isnumerical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
saveMat = p.Results.saveMat;
force = p.Results.force;
binSize = p.Results.binSize;

mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
 'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ30\Final','IZ31\Final',...
 'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline',...
 'IZ33\Saline','IZ34\Saline'};

for ii =1:2
    compiledData{ii}.psthtimes = [];
    compiledData{ii}.psthstccg = [];
    compiledData{ii}.cellID = [];
    compiledData{ii}.mouseID = [];
    compiledData{ii}.sessID = [];
end

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');

    %% Start collecting data
    for ii = 1:size(allSess,1)

        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

        %% Check if assembly units have been detected
        if exist('SleepAssemblyUnits.mat')
            load('SleepAssemblyUnits.mat')
        else 
            continue
        end
        %% Load spikes
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);         

       %% Load pulses
        load([sessionInfo.FileName '.pulses.events.mat']);
        
       %% Extract IDs of CA1 cells
       uid = [];
       for unit = 1:length(spikes.times)
            if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1 || strcmp(cell_metrics.putativeCellType(unit),'Pyramidal Cell')~=1
                continue
            else 
                uid = [uid;cell_metrics.UID(unit)];
            end
       end      
       spikeID = [];
       %% Extract cells that were significantly participating in an assembly
       for assemb = 1:length(assembliesID)
           currUnits = assembliesID{assemb};
           spikeID = [spikeID; uid(currUnits)];
       end
       
       spkID = unique(spikeID);
       
       idxMatch = ismember(uid,spkID);
       spkIDNA  = uid(~idxMatch);
       
       %% Load event timestamps - pulses in the home cage
       pulTr = (pulses.stimComb==2);
       homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
       homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
       pulTr = pulTr & homeCagePulseidx;
       st = pulses.intsPeriods(1,pulTr)';       
       
       %% Calculate PSTH
       for jj = 1:length(spkID)
            % Skip interneurons
            if strcmp(cell_metrics.putativeCellType(spkID(jj)),'Pyramidal Cell')~=1
                continue
            end
            % Get psth
            [stccg, t] = CCG({spikes.times{spkID(jj)} st},[],'binSize',binSize,'duration',21,'normtype','rate');
            compiledData{1}.psthtimes = [compiledData{1}.psthtimes; t'];
            compiledData{1}.psthstccg = [compiledData{1}.psthstccg; stccg(:,2,1)'];
            compiledData{1}.cellID = [compiledData{1}.cellID;spkID(jj)];
            compiledData{1}.mouseID = [compiledData{1}.mouseID; m];
            compiledData{1}.sessID = [compiledData{1}.sessID; ii];            
       end
       
       for jj = 1:length(spkIDNA)
            % Skip interneurons
            if strcmp(cell_metrics.putativeCellType(spkIDNA(jj)),'Pyramidal Cell')~=1
                continue
            end           
            % Get psth
            [stccg, t] = CCG({spikes.times{spkIDNA(jj)} st},[],'binSize',binSize,'duration',21,'normtype','rate');
            compiledData{2}.psthtimes = [compiledData{2}.psthtimes; t'];
            compiledData{2}.psthstccg = [compiledData{2}.psthstccg; stccg(:,2,1)'];
            compiledData{2}.cellID = [compiledData{2}.cellID;spkIDNA(jj)];
            compiledData{2}.mouseID = [compiledData{2}.mouseID; m];
            compiledData{2}.sessID = [compiledData{2}.sessID; ii];            
       end       
    end
end


col = [0.9856 0.7372 0.2537;...
       0.1986 0.7214 0.6310;...
       0.2305 0.2510 0.6173];

% YlGnBu=cbrewer('seq', 'YlGnBu', 11);
% YlGnBu = YlGnBu(11:-1:1,:);
% colormap(YlGnBu)

colormap default

figure
set(gcf,'renderer','painters')
subplot(2,2,1)
imagesc(compiledData{2}.psthtimes(1,:),1:1:size(compiledData{1}.psthstccg,1),zscore(compiledData{1}.psthstccg,[],2))
caxis([-1.5 1.5])

subplot(2,2,2)
imagesc(compiledData{2}.psthtimes(1,:),1:1:size(compiledData{2}.psthstccg,1),zscore(compiledData{2}.psthstccg,[],2))
caxis([-1.5 1.5])

subplot(2,2,[3 4])
plot(compiledData{2}.psthtimes(1,:),nanmean(compiledData{1}.psthstccg,1),'r')
meanpsth = nanmean(compiledData{1}.psthstccg,1);
stdpsth = nanstd(compiledData{1}.psthstccg,[],1)./sqrt(size(compiledData{1}.psthstccg,1));
x = compiledData{2}.psthtimes(1,:);
hold on
fill([x'; fliplr(x)'],[meanpsth'-stdpsth';flipud(meanpsth'+stdpsth')],'r','linestyle','none','FaceAlpha',0.5);

meanpsth = nanmean(compiledData{2}.psthstccg,1);
stdpsth = nanstd(compiledData{2}.psthstccg,[],1)./sqrt(size(compiledData{2}.psthstccg,1));
x = compiledData{2}.psthtimes(1,:);
hold on
fill([x'; fliplr(x)'],[meanpsth'-stdpsth';flipud(meanpsth'+stdpsth')],'k','linestyle','none','FaceAlpha',0.5);
plot(compiledData{2}.psthtimes(1,:),nanmean(compiledData{2}.psthstccg,1),'k')


set(gcf,'renderer','Painters')    
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Revisions\compiledSleepAssembliesPSTH.png'))
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Revisions\compiledSleepAssembliesPSTH.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Revisions\compiledSleepAssembliesPSTH.eps'),'epsc') 

end        