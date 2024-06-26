function getToneMapsAvg7PortsError(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'toneMap',true,@islogical);
addParameter(p,'var','stim',@ischar);

parse(p,varargin{:});
basepath = p.Results.basepath;
toneMap = p.Results.toneMap;
plotfig = p.Results.plotfig;
var = p.Results.var;

%% Deal with inputse
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

fprintf('Computing place fields\n');

% Compute three sets of maps. Forward mapped to space, forward, mapped to
% tone, and reverse. 
%Get the index for different conditions
linIdx = find(behavTrials.linTrial==1);
jumpLin = find(diff(linIdx)>1);

if isempty(jumpLin)
    idx{1}= behavTrials.linTrial==1; % Idx of initial runs  
    idx{30}(1:length(behavTrials.linTrial)) = 0;
    idx{30}(linIdx(1:round(length(linIdx)/2))) = 1; % Idx of first half of initial runs
    
    idx{31}(1:length(behavTrials.linTrial)) = 0;
    idx{31}(linIdx(round(length(linIdx)/2)+1:end)) = 1; % Idx of second half of initial runs
    
    idx{32} = []; % Idx of later runs         
else
    idx{1}(1:length(behavTrials.linTrial)) = 0;
    idx{1}(linIdx(1:jumpLin)) = 1; % Idx of initial runs
    
    idx{30}(1:length(behavTrials.linTrial)) = 0;
    idx{30}(linIdx(1:round(jumpLin/2))) = 1; % Idx of first half of initial runs
    
    idx{31}(1:length(behavTrials.linTrial)) = 0;
    idx{31}(linIdx(round(jumpLin/2)+1:jumpLin)) = 1; % Idx of second half of initial runs
    
    idx{32}(1:length(behavTrials.linTrial)) = 0;
    idx{32}(linIdx(jumpLin+1:end)) = 1; % Idx of later runs      
end
for zz = [1 30 31 32]
    idx{zz} = logical(idx{zz});
end

%No (var) strials
idx{2} = behavTrials.toneGain ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{3} = behavTrials.toneGain ==1 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{4} = behavTrials.toneGain ==2 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{5} = behavTrials.toneGain ==3 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{6} = behavTrials.toneGain ==4 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{7} = behavTrials.toneGain ==5 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{8} = behavTrials.lickLoc ==0 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{9} = behavTrials.lickLoc ==1 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{10} = behavTrials.lickLoc ==2 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{11} = behavTrials.lickLoc ==3 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{12} = behavTrials.lickLoc ==4 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{13} = behavTrials.lickLoc ==5 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;    

% (var) trials
idx{14} = behavTrials.toneGain ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{15} = behavTrials.toneGain ==1 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{16} = behavTrials.toneGain ==2 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{17} = behavTrials.toneGain ==3 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{18} = behavTrials.toneGain ==4 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{19} = behavTrials.toneGain ==5 & behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{20} = behavTrials.lickLoc ==0 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{21} = behavTrials.lickLoc ==1 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{22} = behavTrials.lickLoc ==2 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{23} = behavTrials.lickLoc ==3 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{24} = behavTrials.lickLoc ==4 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{25} = behavTrials.lickLoc ==5 & behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;  

%correct vs incorrect
idx{26} = behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;
idx{27} = behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==0;  
idx{28} = behavTrials.correct==0 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;
idx{29} = behavTrials.correct==1 & behavTrials.linTrial ==0 & behavTrials.(var) ==1;  

gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
freqExp = log10(22000/1000);

for ii = 1:length(idx)    
    if ii<26
        [idxPos] = InIntervals(tracking.timestamps,behavTrials.timestamps(idx{ii}(1:(end-1)),:));   
        positions.forward{ii} = [tracking.timestamps(idxPos) tracking.position.x(idxPos) tracking.position.y(idxPos)];
        if ii==1
            positions.tone{ii} = [tracking.timestamps(idxPos) tracking.position.x(idxPos) tracking.position.y(idxPos)*nan];
            kk = 0;
        elseif ii<=7
            kk = ii-1;
        elseif ii<=13
            kk = ii-7;
        elseif ii<=19
            kk = ii-13;
        elseif ii<=25
            kk = ii-19;            
        end

        if kk >0
            y = tracking.position.y(idxPos); 
            tonepos = [];
            for jj = 1:length(y)
                tonepos(jj) = (y(jj)*gain(kk))/120;
                %tonepos(jj) = (1000*(10.^(freqExp*freq)));
            end
            tonepos(tonepos>1.1) = nan;        
            if isempty(tracking.timestamps(idxPos))
                positions.tone{ii} = [tracking.timestamps(idxPos) tracking.position.x(idxPos) tracking.position.y(idxPos)];    
            else
                positions.tone{ii} = [tracking.timestamps(idxPos) tracking.position.x(idxPos) tonepos'];    
            end
                
        end  
        if ii == 1 || ii == 7
            idPos = find(idx{ii}(1:(end-1)));
            % Only take moments when vy is negative
            [~,~,~,vx,vy,~,~] = KalmanVel(tracking.position.x,tracking.position.y,tracking.timestamps,2);
            x_temp = tracking.position.x;
            x_temp(vy>=-2) = nan;
            
            y_temp = tracking.position.y;
            y_temp(vy>=-2) = nan;
            
            [idxPos] = InIntervals(tracking.timestamps,[behavTrials.timestamps(idPos,2) behavTrials.timestamps(idPos+1,1)]);
            if ii== 1
                positions.reverse{ii}= [tracking.timestamps(idxPos) x_temp(idxPos) y_temp(idxPos)];             
            elseif ii == 7
                positions.reverse{2}= [tracking.timestamps(idxPos) x_temp(idxPos) y_temp(idxPos)];    
            end
        end
    elseif ii<30
        idPos = find(idx{ii}(1:(end-1)));
        [idxPos] = InIntervals(tracking.timestamps,[behavTrials.timestamps(idPos,2) behavTrials.timestamps(idPos+1,1)]);  
        positions.reverse{ii-23}= [tracking.timestamps(idxPos) x_temp(idxPos) y_temp(idxPos)]; 
    else
        [idxPos] = InIntervals(tracking.timestamps,behavTrials.timestamps(idx{ii}(1:(end-1)),:));
        if ~isempty(idxPos)
            positions.forward{ii-4} = [tracking.timestamps(idxPos) tracking.position.x(idxPos) tracking.position.y(idxPos)];
        end
        
        idPos = find(idx{ii}(1:(end-1)));
        [idxPos] = InIntervals(tracking.timestamps,[behavTrials.timestamps(idPos,2) behavTrials.timestamps(idPos+1,1)]);  
        positions.reverse{ii-23}= [tracking.timestamps(idxPos) x_temp(idxPos) y_temp(idxPos)];        
    end
    
end

firingMaps.forward = bz_getRateMaps(positions.forward,spikes,'xRange',[0 6],'yRange',[0 125], 'binSize',2.5,'saveMat',false);
if toneMap
   firingMaps.tone = bz_getRateMaps(positions.tone,spikes,'xRange',[0 6],'yRange',[0 1.1], 'binSize',0.022,'minOccupancy',0,'saveMat',false);    
end
firingMaps.reverse = bz_getRateMaps(positions.reverse,spikes,'xRange',[0 6],'yRange',[0 125], 'binSize',2.5,'saveMat',false);

firingMaps.linTrial = behavTrials.linTrial(1:(end-1));
firingMaps.toneTrial = behavTrials.toneTrial(1:(end-1));
firingMaps.toneGain = behavTrials.toneGain(1:(end-1));
firingMaps.correct = behavTrials.correct(1:(end-1));
firingMaps.numLicks = behavTrials.numLicks(1:(end-1),:);

save([sessionInfo.FileName '.rateMapsAvgError.cellinfo.mat'],'firingMaps'); 

labels = {'forward','tone','reverse'};
col = [176/243 223/243 229/243; 149/243 200/243 216/243; 137/243 207/243 240/243;
        70/243 130/243 180/243; 16/243 52/243 166/243;0/243 0/243 128/243];
colRev = [0.5 0.5 0.5; 175/243 54/243 60/243; 70/243 148/243 73/243; 175/243 54/243 60/243; 70/243 148/243 73/243];
a = linspace(2000,22000,50);
b = linspace(0,125,50);

if plotfig
    
    if ~isfolder('FiringMapAvg')
        mkdir('FiringMapAvg')
    end
    for pf = 1:length(firingMaps.forward.rateMaps)
       figure
       set(gcf,'Renderer','painters')
       set(gcf,'Color','w')
       set(gcf,'Position',[1154 358 1682 394])
       %plot linear track   
       subplot(4,6,1)
       dataMat = firingMaps.forward.rateMaps{pf}{1};
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('Linear Track')
       subplot(4,6,7)
       plot(b, dataMat,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
       
       dataMat = [];
       dataMatTone = [];
       for ii = 2:7           
           dataMat = [dataMat;firingMaps.forward.rateMaps{pf}{ii}];
           subplot(4,6,8)
           plot(b, firingMaps.forward.rateMaps{pf}{ii},'Color',col(ii-1,:),'LineWidth',1.5)
           hold on
           
           dataMatTone = [dataMatTone;firingMaps.tone.rateMaps{pf}{ii}];
           subplot(4,6,9)
           plot(a, firingMaps.tone.rateMaps{pf}{ii},'Color',col(ii-1,:),'LineWidth',1.5)           
           hold on
       end
       subplot(4,6,2)
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('Space, correct')
       
       subplot(4,6,3)
       h = imagesc(a,size(dataMatTone,2),dataMatTone);
       set(h, 'AlphaData', ~isnan(dataMatTone))
       title('Tone, correct')
       
       dataMat = [];
       dataMatTone = [];
       for ii = 8:13       
           dataMat = [dataMat;firingMaps.forward.rateMaps{pf}{ii}];
           subplot(4,6,10)
           plot(b, firingMaps.forward.rateMaps{pf}{ii},'Color',col(ii-7,:),'LineWidth',1.5)
           hold on
           
           dataMatTone = [dataMatTone;firingMaps.tone.rateMaps{pf}{ii}];
           subplot(4,6,11)
           plot(a, firingMaps.tone.rateMaps{pf}{ii},'Color',col(ii-7,:),'LineWidth',1.5)     
           hold on
       end
       subplot(4,6,4)
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('Space, incorrect')
       
       subplot(4,6,5)
       h = imagesc(a,size(dataMatTone,2),dataMatTone);
       set(h, 'AlphaData', ~isnan(dataMatTone))
       title('Tone, incorrect')       
       
       dataMat = [];
       for ii = 1:5
           dataMat = [dataMat;firingMaps.reverse.rateMaps{pf}{ii}];
           subplot(4,6,12)
           plot(a, firingMaps.reverse.rateMaps{pf}{ii},'Color',colRev(ii,:),'LineWidth',1.5) 
           hold on
       end
       subplot(4,6,6)
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('Return')     

       dataMat = [];
       dataMatTone = [];
       for ii = 14:19           
           dataMat = [dataMat;firingMaps.forward.rateMaps{pf}{ii}];
           subplot(4,6,20)
           plot(b, firingMaps.forward.rateMaps{pf}{ii},'Color',col(ii-13,:),'LineWidth',1.5)
           hold on
           
           dataMatTone = [dataMatTone;firingMaps.tone.rateMaps{pf}{ii}];
           subplot(4,6,21)
           plot(a, firingMaps.tone.rateMaps{pf}{ii},'Color',col(ii-13,:),'LineWidth',1.5)           
           hold on
       end
       subplot(4,6,14)
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('(var) Space, correct')
       
       subplot(4,6,15)
       h = imagesc(a,size(dataMatTone,2),dataMatTone);
       set(h, 'AlphaData', ~isnan(dataMatTone))
       title('(var) Tone, correct')
       
       dataMat = [];
       dataMatTone = [];
       for ii = 20:25       
           dataMat = [dataMat;firingMaps.forward.rateMaps{pf}{ii}];
           subplot(4,6,22)
           plot(b, firingMaps.forward.rateMaps{pf}{ii},'Color',col(ii-19,:),'LineWidth',1.5)
           hold on
           
           dataMatTone = [dataMatTone;firingMaps.tone.rateMaps{pf}{ii}];
           subplot(4,6,23)
           plot(a, firingMaps.tone.rateMaps{pf}{ii},'Color',col(ii-19,:),'LineWidth',1.5)     
           hold on
       end
       subplot(4,6,16)
       h = imagesc(b,size(dataMat,2),dataMat);
       set(h, 'AlphaData', ~isnan(dataMat))
       title('(var) Space, incorrect')
       
       subplot(4,6,17)
       h = imagesc(a,size(dataMatTone,2),dataMatTone);
       set(h, 'AlphaData', ~isnan(dataMatTone))
       title('(var) Tone, incorrect')        
       
       saveas(gcf,['FiringMapAvg',filesep ,'cell_' num2str(pf) '.png'],'png');
       saveas(gcf,['FiringMapAvg',filesep ,'cell_' num2str(pf) '.eps'],'epsc');
       saveas(gcf,['FiringMapAvg',filesep ,'cell_' num2str(pf) '.fig'],'fig');
       close all;

    end   
end

% Calculate correlation
% for pf = 1:length(firingMaps.forward.rateMaps)
%     corrSpace = [];
%     corrTone = [];
%     dataMat = [];
%     dataMatTone = [];
%     for ii = 3:7           
%         dataMat = [dataMat;firingMaps.forward.rateMaps{pf}{ii}];   
%         a = fillmissing(firingMaps.tone.rateMaps{pf}{ii},'linear');
%         dataMatTone = [dataMatTone;a];
%     end    
%     for ii = 1:5
%        for jj = (ii+1):5
%            a = corrcoef(dataMat(ii,:),dataMat(jj,:),'rows','complete');
%            corrSpace = [corrSpace a(1,2)];
%            a = corrcoef(dataMatTone(ii,:),dataMatTone(jj,:),'rows','complete');         
%            corrTone = [corrTone a(1,2)];
%        end
%     end
%     avgRate  = nanmean(nanmean(dataMat,1));
%     maxRate = max(max(dataMat));
%     if avgRate<15 & maxRate>4
%         corrSpatial(pf) = nanmean(corrSpace);
%         corrFrequency(pf) = nanmean(corrTone);
%     else
%         corrSpatial(pf) = nan;
%         corrFrequency(pf) = nan;
%     end
%         
% end
% 
% figure
% scatter(corrSpatial,corrFrequency)
% hold on
% ylim([-0.4 1])
% refline(1)
% 
% find(corrFrequency>corrSpatial)
end