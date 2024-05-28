function periLickSpeed

% sesstoAnalyze = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
%     'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
%     'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
%     'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
%     'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
%     'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
%     'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
%     'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
%     'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
%     'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
%     'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
%     'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
%     'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 
%     'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
%     'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
%     'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
%     'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28',... 
%     }; 

sesstoAnalyze = {'IZ47\Final\IZ47_230710_sess25'};

plotfig = 1;
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';

col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255;...
    0 0 0];

win = 5*30;

for ll = 1:7
    speedAvg{ll} = [];
end

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name);
    
    file = dir('*.Tracking.Behavior.mat');
    load(file.name);

    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.correct),1) = 0;
    end

    if ~isfield(behavTrials,'stim')
        behavTrials.stim(1:length(behavTrials.correct),1) = 0;
    end

    % Get ts by lickport
    for ll = 1:6
        if ll==7
            idx = behavTrials.linTrial==1 & behavTrials.stim==0 & behavTrials.probe==0 ;
        else
            idx = behavTrials.linTrial==0 & behavTrials.stim==0 & behavTrials.probe==0  & behavTrials.lickLoc==ll-1 ;
        end

        lickTS = behavTrials.timestamps(idx,2);
        if ~isempty(lickTS)
            speedProfile = [];
            for tt = 1:(length(lickTS)-1)
                [~,idx] = min(abs(tracking.timestamps-(lickTS(tt)-0.033)));
                speedProfile(tt,:) = tracking.position.y(idx-win:idx+win);
            end
            speedAvg{ll}(ii,:) = nanmean(speedProfile,1);
        end
    end
end

f1 = figure;
set(gcf,'Renderer','painters')

xAxis = linspace(-5,5,301);

for ll = 1:7
    hold on
    plotAvgStd(speedAvg{ll},1,1,1,f1,xAxis',col(ll,:),0)
end
line([0 0],[0 45])
line([-1 -1],[0 45])
end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col, useMedian)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);
    if ~useMedian
        meanpsth = nanmean(array,1);
        stdpsth = nanstd(array,1)./sqrt(size(array,1));
        lArr  = meanpsth-stdpsth;
        uArr = meanpsth+stdpsth;
    else
        meanpsth = nanmedian(array,1);
        a = quantile(array,4,1);
        lArr  = meanpsth-a(2,:);
        uArr = meanpsth+a(3,:);
    end
    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end
