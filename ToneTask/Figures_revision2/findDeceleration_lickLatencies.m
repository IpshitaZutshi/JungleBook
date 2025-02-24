function Deceleration = findDeceleration_lickLatencies(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'plotIndfig',false,@islogical);
addParameter(p,'error',1,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
plotfig = p.Results.plotfig;
plotIndfig = p.Results.plotIndfig;
error = p.Results.error;

cd(basepath)

if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

%% Within each trial, only extract timestamps after the mouse has crossed the first port
acc = gradient(tracking.position.v)./0.033;

binSize = 0.033; % 33 ms in seconds
kernelWidth = 0.2; % 200 ms in seconds
numBinsInKernel = round(kernelWidth / binSize); % Number of bins in the kernel

% Create the Gaussian kernel
gaussianKernel = gausswin(numBinsInKernel);

% Normalize the Gaussian kernel to have unit area
gaussianKernel = gaussianKernel / sum(gaussianKernel);

Deceleration.ts = [];

for tt = 1:length(behavTrials.linTrial)
    if behavTrials.linTrial(tt)==0 && behavTrials.lickLoc(tt) >=2 && behavTrials.correct(tt) == error
        idx = tracking.timestamps>behavTrials.timestamps(tt,1) & tracking.timestamps<(behavTrials.timestamps(tt,2)+0.2);
        curVel = tracking.position.v(idx);
        curAcc = acc(idx);
        curPos = tracking.position.y(idx);
        curTS = tracking.timestamps(idx);

        curAcc1 = curAcc(curPos>5);      
        curTS1 = curTS(curPos>5);
        curPos1 = curPos(curPos>5);

        if length(curAcc1)<15
                continue
        end
        smoothAcc = conv(curAcc1, gaussianKernel, 'same');
        [pks,locs1] = findpeaks(-smoothAcc,'MinPeakProminence',10,'MinPeakHeight',10);

        %% Make sure that the decelerations fall between -5 s and 0.5 s from the end of the trial
        idxTS= zeros(1,length(locs1));

        for ll = 1:length(locs1)
            [~,idx] = min(abs(behavTrials.timestamps(:,2)-curTS1(locs1(ll))));
            ttLick = curTS1(locs1(ll))-behavTrials.timestamps(idx,2);
            if ttLick>=-5 && ttLick<=0.5
                idxTS(ll) = 1;
            end
        end
        locs = locs1(logical(idxTS));
        Deceleration.ts = [Deceleration.ts;curTS1(locs)];

        if plotIndfig
            figure
            subplot(2,1,1)
            plot(curTS1, curAcc1)
            hold on
            plot(curTS1, smoothAcc)
            scatter(curTS1(locs),smoothAcc(locs),'o')

            subplot(2,1,2)
            plot(curPos1, curAcc1)
            hold on
            plot(curPos1, smoothAcc)
            scatter(curPos1(locs),smoothAcc(locs),'o')
        end
    end
end

Deceleration.velPSTH = [];
Deceleration.accPSTH = [];
Deceleration.posX = [];
Deceleration.posY = [];
Deceleration.decType = [];
Deceleration.timetoLick1 = [];
timetoLick2 = [];
v = tracking.position.vy;
acc = gradient(tracking.position.v)./0.033;

portPos = [53 59; 76 83; 100 108; 116 125];

for tt  = 1:length(Deceleration.ts)
    [~,idx] = min(abs(tracking.timestamps-Deceleration.ts(tt)));
    Deceleration.velPSTH = [Deceleration.velPSTH; v(idx-60:idx+150)'];
    Deceleration.accPSTH = [Deceleration.accPSTH; acc(idx-60:idx+150)'];
    Deceleration.posX = [Deceleration.posX tracking.position.x(idx)];
    Deceleration.posY = [Deceleration.posY tracking.position.y(idx)];
    % Calculate time to the closest lick
    [~,idxLick] = min(abs(behavTrials.timestamps(:,2)-Deceleration.ts(tt)));
    Deceleration.timetoLick1 = [Deceleration.timetoLick1 Deceleration.ts(tt)-behavTrials.timestamps(idxLick,2)];
    curlickloc = behavTrials.lickLoc(idxLick);
    % if Deceleration.timetoLick1(tt)>=-0.5
    %     Deceleration.decType = [Deceleration.decType 1];
    if ((Deceleration.posY(tt)>28 && Deceleration.posY(tt) <34) || (Deceleration.posY(tt)>53 && Deceleration.posY(tt)<59) || ...
            (Deceleration.posY(tt)>76 && Deceleration.posY(tt) <83) || (Deceleration.posY(tt) >100 && Deceleration.posY(tt) <108) ||...
            (Deceleration.posY(tt) >116)) && ...
            (Deceleration.posY(tt)> portPos(curlickloc-1,1) && Deceleration.posY(tt)< portPos(curlickloc-1,2)) && ...
            abs(Deceleration.ts(tt)-behavTrials.timestamps(idxLick,2))<0.5% The deceleration is near the eventual lick port 
        Deceleration.decType = [Deceleration.decType 1];
    elseif Deceleration.posX(tt) <=3 && ((Deceleration.posY(tt)>28 && Deceleration.posY(tt) <34) || (Deceleration.posY(tt)>53 && Deceleration.posY(tt)<59) || ...
            (Deceleration.posY(tt)>76 && Deceleration.posY(tt) <83) || (Deceleration.posY(tt) >100 && Deceleration.posY(tt) <108) ||...
            (Deceleration.posY(tt) >116)) && ...
            ~(Deceleration.posY(tt)> portPos(curlickloc-1,1) && Deceleration.posY(tt)< portPos(curlickloc-1,2)) % The deceleration is not near the eventual lick port 
        Deceleration.decType = [Deceleration.decType 2];   
    else
        Deceleration.decType = [Deceleration.decType 3];
    end

end

for tt = 1:length(behavTrials.lickLoc)
    [~,idxLick] = min(abs(Deceleration.ts - behavTrials.timestamps(tt,2)));
    timetoLick2 = [timetoLick2 behavTrials.timestamps(tt,2)-Deceleration.ts(idxLick)];
end

if plotfig
    fig2= figure;
    tAxis = linspace(-2,4,181);

    plotAvgStd(Deceleration.velPSTH,1,3,1,fig2,tAxis',[0 0 1])
    title('Speed')

    plotAvgStd(Deceleration.accPSTH,1,3,2,fig2,tAxis',[0 0 1])
    title('Acceleration')

    subplot(1,3,3)
    histogram(Deceleration.timetoLick1(Deceleration.decType==1),-1:0.02:3)
end

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    stdpsth = nanstd(array,1)./sqrt(size(array,1));
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end