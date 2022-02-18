function getPhotometryMaps(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;


 %% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    for ii = 1:size(MergePoints.foldernames,2)
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'test*']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            fprintf('Computing [photometry signal in %s folder \n',MergePoints.foldernames{ii});           
            photometry = getPhotometry;
            
            file = dir(['*.Tracking.Behavior.mat']);
            load(file.name);
            
            file = dir(['*.TrialBehavior.Events.mat']);
            load(file.name);
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end

%% Now separate trial types
for pf = 1:(size(behavTrials.timestamps,1)-1)   
    %forward run
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    % Find closest photometry signal
    photo = [];
    idxTS = find(idx);
    for ii  = 1:length(idxTS)
        currTS = tracking.timestamps(idxTS(ii));
        [~,idxP] = min(abs(photometry.timestamps-currTS));
        photo(ii) = photometry.corrSignal(idxP);
    end
    positions.forward{pf} = [tracking.position.y(idx) tracking.position.v(idx) photo'];

    %return run
    [idx] = InIntervals(tracking.timestamps,[behavTrials.timestamps(pf,2) behavTrials.timestamps(pf+1,1)]);
    % Find closest photometry signal
    photo = [];
    idxTS = find(idx);
    for ii  = 1:length(idxTS)
        currTS = tracking.timestamps(idxTS(ii));
        [~,idxP] = min(abs(photometry.timestamps-currTS));
        photo(ii) = photometry.corrSignal(idxP);
    end
    positions.return{pf} = [tracking.position.y(idx) tracking.position.v(idx) photo'];  
end

linMap = linspace(1,112,56); 
% First lin trials
idx = find(behavTrials.linTrial);
%figure
for kk = 1:7
    photoMap{kk} = [];
    photoSpeed{kk} = [];
end
for ii = 1:length(idx)
    % First photometry and spatial position
    trialMap = [];
    for bb = 1:(length(linMap)-1)
        idxLoc = positions.forward{idx(ii)}(:,1)>=linMap(bb) & positions.forward{idx(ii)}(:,1)<linMap(bb+1);
        trialMap(bb)  = nanmean(positions.forward{idx(ii)}(idxLoc,3));
    end
    photoMap{1} = [photoMap{1}; trialMap];
    
    % Next photometry and speed
    data = [positions.forward{idx(ii)}(:,2) positions.forward{idx(ii)}(:,3)]; 
    photoSpeed{1} = [photoSpeed{1}; data];
%     subplot(4,4,1)
%     plot(positions.forward{idx(ii)}(:,1),positions.forward{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%     hold on
%     subplot(4,4,5)
%     plot(positions.return{idx(ii)}(:,1),positions.return{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%     hold on
end

for kk = 1:3
    idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & behavTrials.toneGain(1:(end-1)) ==(kk-1));
    for ii = 1:length(idx)        
        % First photometry and spatial position
        trialMap = [];
        for bb = 1:(length(linMap)-1)
            idxLoc = positions.forward{idx(ii)}(:,1)>=linMap(bb) & positions.forward{idx(ii)}(:,1)<linMap(bb+1);
            trialMap(bb)  = nanmean(positions.forward{idx(ii)}(idxLoc,3));
        end
        photoMap{1+kk} = [photoMap{1+kk}; trialMap];

        % Next photometry and speed
        data = [positions.forward{idx(ii)}(:,2) positions.forward{idx(ii)}(:,3)]; 
        photoSpeed{1+kk} = [photoSpeed{1+kk}; data];
%         subplot(4,4,1+kk)        
%         plot(positions.forward{idx(ii)}(:,1),positions.forward{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%         hold on
%         
%         subplot(4,4,5+kk)        
%         plot(positions.return{idx(ii)}(:,1),positions.return{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%         hold on        
    end
end

for kk = 1:3
    idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==0 & behavTrials.toneGain(1:(end-1)) ==(kk-1));
    for ii = 1:length(idx)        
        % First photometry and spatial position
        trialMap = [];
        for bb = 1:(length(linMap)-1)
            idxLoc = positions.forward{idx(ii)}(:,1)>=linMap(bb) & positions.forward{idx(ii)}(:,1)<linMap(bb+1);
            trialMap(bb)  = nanmean(positions.forward{idx(ii)}(idxLoc,3));
        end
        photoMap{4+kk} = [photoMap{4+kk}; trialMap];

        % Next photometry and speed
        data = [positions.forward{idx(ii)}(:,2) positions.forward{idx(ii)}(:,3)]; 
        photoSpeed{4+kk} = [photoSpeed{4+kk}; data];
%         subplot(4,4,1+kk)        
%         plot(positions.forward{idx(ii)}(:,1),positions.forward{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%         hold on
%         
%         subplot(4,4,5+kk)        
%         plot(positions.return{idx(ii)}(:,1),positions.return{idx(ii)}(:,3),'Color',[0.5 0.5 0.5])
%         hold on        
    end
end
save([basepath filesep sessionInfo.session.name '.photoMapsTrial.mat'],'photoMap','photoSpeed');

figure
set(gcf,'Color','w')
linMap = linspace(1,112,56); 
for kk = 1:7
    if kk == 1
        subplot(2,4,1)
        col = [0.5 0.5 0.5];
    elseif kk >1 && kk <5
        subplot(2,4,kk)
        col = [8/243 133/243 161/243];
    else
        subplot(2,4,5+(kk-4))
        col = [175/243 54/243 60/243];
    end
    imagesc(linMap(1:(end-1)),1:size(photoMap{kk}),photoMap{kk})
   % caxis([0 0.008])
    
    subplot(2,4,5)   
    plot(linMap(1:(end-1)), nanmean(photoMap{kk}),'Color',col,'LineWidth',1.5)
    hold on
    xlabel('Position on track')
end
end