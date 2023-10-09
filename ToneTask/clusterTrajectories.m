function [posXtrial,posXLabel,T,positions,lin_begin,lin_end,tone6,posX,posY] = clusterTrajectories(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'plotfig',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
plotfig = p.Results.plotfig;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
end

% Port locations

%portLoc = [1 5 12 20 32 40 55 80 103 120];
portLoc = [1 11.6, 32.27 55.53 79.62 102.79 120];

% Initial lin trials
linIdx = find(behavTrials.linTrial==1);
jumpLin = find(diff(linIdx)>1);  
if isempty(jumpLin)
    lin_begin = linIdx;
else
    lin_begin = linIdx(1:jumpLin);
end

if isempty(jumpLin)
    lin_end = [];
else
    lin_end = linIdx(jumpLin+1:end);
end

tone6 = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & ...
            behavTrials.toneGain(1:(end-1)) == 5);
if plotfig
    figure
end
for ii = 1:length(lin_begin)    
   if plotfig
    subplot(3,3,1)
    hold on
    plot(positions.forward{lin_begin(ii)}(:,2),positions.forward{lin_begin(ii)}(:,1),'Color',[0.6 0.6 0.6])
    title(basepath)
   end
   % Downsample to seven points, corresponding to each port
   y = positions.forward{lin_begin(ii)}(:,2);
   x = positions.forward{lin_begin(ii)}(:,1);
   for jj = 1:length(portLoc)
       [~,idClose] = min(abs(y-portLoc(jj)));
       posLinX(ii,jj) = x(idClose);
       posLinY(ii,jj) = y(idClose);      
   end
   if plotfig
       subplot(3,3,2)
       hold on   
       plot(posLinY(ii,:),posLinX(ii,:),'Color',[0.6 0.6 0.6])
   end
end

for ii = 1:(length(tone6))
   if plotfig
    subplot(3,3,4)
    hold on    
    plot(positions.forward{tone6(ii)}(:,2),positions.forward{tone6(ii)}(:,1),'Color',[56/243 61/243 150/243])
   end
   y = positions.forward{tone6(ii)}(:,2);
   x = positions.forward{tone6(ii)}(:,1);
   for jj = 1:length(portLoc)
       [~,idClose] = min(abs(y-portLoc(jj)));
       posTone6X(ii,jj) = x(idClose);
       posTone6Y(ii,jj) = y(idClose);      
   end
   if plotfig
    subplot(3,3,5)
    hold on   
    plot(posTone6Y(ii,:),posTone6X(ii,:),'Color',[56/243 61/243 150/243])   
   end
end

posLinEndX = [];
posLinEndY = [];
for ii = 1:(length(lin_end)-1)
   if plotfig
    subplot(3,3,7)
    hold on    
    plot(positions.forward{lin_end(ii)}(:,2),positions.forward{lin_end(ii)}(:,1),'Color',[0 0 0])
   end
   
   y = positions.forward{lin_end(ii)}(:,2);
   x = positions.forward{lin_end(ii)}(:,1);
   for jj = 1:length(portLoc)
       [~,idClose] = min(abs(y-portLoc(jj)));
       posLinEndX(ii,jj) = x(idClose);
       posLinEndY(ii,jj) = y(idClose);      
   end   

   if plotfig
    subplot(3,3,8)
    hold on   
    plot(posLinEndY(ii,:),posLinEndX(ii,:),'Color',[0 0 0])      
   end
end

%Group trajectories and label them
posX = [posLinX;posTone6X;posLinEndX];
posY = [posLinY;posTone6Y;posLinEndY];
posXLabel = [ones(size(posLinX,1),1)*1;ones(size(posTone6X,1),1)*2;ones(size(posLinEndX,1),1)*3];
posXtrial = [lin_begin; tone6; lin_end(1:(end-1))];

Z =linkage(posX,'complete');
T = cluster(Z,'Maxclust',2);

if plotfig
    subplot(3,3,3)
    hold on
    for ii = 1:(length(unique(T)))
        idx = find(T==ii & posXLabel==1);
        if ii == 1
            col = [0.6 0.6 0.6];
        elseif ii == 2
            col = [56/243 61/243 150/243];
        else
            col = [0 0 0];
        end

        for jj = 1:length(idx)
            plot(posY(idx(jj),:), posX(idx(jj),:),'Color',col)
        end
    end

    subplot(3,3,6)
    hold on
    for ii = 1:(length(unique(T)))
        idx = find(T==ii & posXLabel==2);
        if ii == 1
            col = [0.6 0.6 0.6];
        elseif ii == 2
            col = [56/243 61/243 150/243];
        else
            col = [0 0 0];
        end

        for jj = 1:length(idx)
            plot(posY(idx(jj),:), posX(idx(jj),:),'Color',col)
        end
    end

    subplot(3,3,9)
    hold on
    for ii = 1:(length(unique(T)))
        idx = find(T==ii & posXLabel==3);
        if ii == 1
            col = [0.6 0.6 0.6];
        elseif ii == 2
            col = [56/243 61/243 150/243];
        else
            col = [0 0 0];
        end

        for jj = 1:length(idx)
            plot(posY(idx(jj),:), posX(idx(jj),:),'Color',col)
        end
    end
end

if plotfig
    saveLoc = strcat(basepath,'\SummaryFigures');
    if ~isfolder('SummaryFigures')
        mkdir('SummaryFigures')
    end  
   saveas(gcf,[saveLoc,filesep ,'ClusteredTrajectory.png'],'png');
   saveas(gcf,[saveLoc,filesep ,'ClusteredTrajectory.fig'],'fig');
end    

end