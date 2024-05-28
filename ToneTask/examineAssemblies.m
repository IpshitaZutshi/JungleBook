function Activities = examineAssemblies

sessloc = pwd;%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)

file = dir(['*TrialBehavior.Behavior.mat']);
load(file(1).name);
file = dir(['*Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*spikes.cellinfo.mat']);
load(file(1).name);
file = dir(['*cell_metrics.cellinfo.mat']);
load(file.name);
file = dir(['*rateMapsAvg.cellinfo.mat']);
load(file(1).name);


%% First define assemblies in behavior
dtime = mean(diff(tracking.timestamps));
spkData = bz_SpktToSpkmat(spikes.times,'dt',dtime, 'win',[behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
spikeMat = spkData.data';

%Only select pyramidal cells
logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
%Only select cells with a rate> 0.1 Hz
rate = sum(spikeMat,2)./(length(spkData.timestamps)*(1/30));
logVector2 = rate>0.1;
keepCells = logicalVector'& logVector2;
spikeMat = spikeMat(keepCells,:);

SpikeCount = zscore(spikeMat,[],1);
AssemblyTemplates = assembly_patterns(SpikeCount);

for aa = 1:size(AssemblyTemplates,2)
    assembliesVector = AssemblyTemplates(:,aa);
    % flip if negative weights:
    assembliesVector = assembliesVector * assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector)))...
        /abs(assembliesVector(find(max(abs(assembliesVector)) == abs(assembliesVector))));
    Vectors(:,aa) = assembliesVector;

    %For each assembly, define significant cells
    assembliesID{aa} = find(assembliesVector> 2*std(assembliesVector));
   % [~,maxCell(aa)] = max(assembliesVector);

  %  ordCell(aa) = find(sortIdx_new==maxCell(aa));
end   

% Interesting assemblies in IZ43_sess4, 21 and 17

%% For each assembly, calculate average forward and return ratemaps for
% NT1, TONE, and NT2, Return (averaged across significant cells for the assembly)
Activities = assembly_activity(Vectors,SpikeCount);

%For each assembly, extract times when it was "active"
for aa = 1:size(Activities,1)
    idx = find(Activities(aa,:)>3*std(Activities(aa,:)));
    assemblyTimes{aa} = spkData.timestamps(idx);
    for ii = 1:8
        %calculate ccg of each assembly
        if ii == 7
            [stccg, t] = CCG({assemblyTimes{aa} behavTrials.timestamps(behavTrials.linTrial==1,2)},[],'binSize',0.05,'duration',6,'norm','rate');
        else
            [stccg, t] = CCG({assemblyTimes{aa} behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.lickLoc==ii-1,2)},[],'binSize',0.05,'duration',10,'norm','rate');
        end

        psth{ii}(aa,:) = stccg(:,2,1);
    end
    
end


% for aa = 1:length(assembliesID) 
%     NT1 = [];
%     NT2 = [];
%     Return = [];
%     ToneAll = [];   
%     for bb = 1:length(assembliesID{aa})
%         locid = assembliesID{aa}(bb);
%         NT1 = [NT1;firingMaps.forward.rateMaps{locid}{1}];
%         if length(firingMaps.forward.rateMaps{locid})==28
%             NT2 = [NT2;firingMaps.forward.rateMaps{locid}{28}];
%         end
%         % Average Tone
%         dataMat = [];
%         for jj = 2:7           
%             dataMat = [dataMat;firingMaps.forward.rateMaps{locid}{jj}];
%         end 
%         ToneAll = [ToneAll;nanmean(dataMat,1)];
%         Return = [Return;firingMaps.reverse.rateMaps{locid}{2}];
%     end
%     Maps.NT1(aa,:) = nanmean(NT1,1);
%     Maps.NT2(aa,:) = nanmean(NT2,1);
%     Maps.Return(aa,:) = nanmean(Return,1);
%     Maps.Tone(aa,:) = nanmean(ToneAll,1);
% end

[~,maxPSTH]  = max(psth{3},[],2);
[~,sortidx] = sort(maxPSTH,'descend');
figure
%[~,sortidx] = sort(ordCell,'ascend');
for ii = 1:7
    subplot(1,7,ii)    
    imagesc(t, 1:size(psth{ii},1),psth{ii})%(sortidx,:))
    caxis([0 7])
    set(gca,'YDir','normal')
    hold on
    line([0 0],[1 22],'LineWidth',1.5,'Color','w')
end

startT = behavTrials.timestamps(45,1);
endT = behavTrials.timestamps(60,2);
[~,startIdx] = min(abs(spkData.timestamps-startT));
[~,endIdx] = min(abs(spkData.timestamps-endT));

figure
for kk = 1:22
    plot(spkData.timestamps(startIdx:endIdx),Activities(kk,startIdx:endIdx))
    hold on
end
for ii = 45:60
    line([behavTrials.timestamps(ii,2) behavTrials.timestamps(ii,2)],[-10 60])
end
plot(tracking.timestamps,tracking.position.y)
plot(tracking.timestamps,tracking.position.vy)
scatter(tracking.timestamps(idxFrame),tracking.position.y(idxFrame),10,'r')
xlim([startT endT])

%Include assembly times only within the forward run, happening 0.5 seconds BEFORE the lick
idxTone = behavTrials.linTrial==0;
idx1 = InIntervals(spikes.times{28},[behavTrials.timestamps(idxTone,1) behavTrials.timestamps(idxTone,2)-1]);
ids = find(idx1);
figure
plot(tracking.timestamps,tracking.position.y)
hold on
plot(tracking.timestamps,tracking.position.v)
vMat = []; 
idxFrame = [];
for ii = 1:length(ids)
    [~,posid] = min(abs(tracking.timestamps-spikes.times{28}(ids(ii))));
    if tracking.position.y(posid)>12 && tracking.position.y(posid)<98
        idxFrame = [idxFrame posid];
        scatter(tracking.timestamps(posid),tracking.position.y(posid),10,'r','filled')    
        vMat = [vMat; tracking.position.v(posid-66:posid+66)'];   
    end
end
figure
plot(nanmean(vMat,1))
for ii = 25
    figure
    for jj = 1:3
        subplot(1,3,jj)
        frame = read(videoObj1,idxFrame(ii)+10*(jj-2));
        imagesc(frame(50:end,200:450,:))  
    end
end

end