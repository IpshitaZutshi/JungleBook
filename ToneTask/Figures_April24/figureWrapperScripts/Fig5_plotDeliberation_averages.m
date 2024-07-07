function Fig5_plotDeliberation_averages

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[463 2 1261 947]);

numrows = 3;
numcol = 4;

%% Now, in panel A, show the firing of "tone cells"
%Summary = linkToneCellFiringtoDeceleration('plotfig',false);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\DeliberationPSTHSummary.mat')

spec=cbrewer('seq', 'Blues', 20);
spec(spec>1) = 1;
spec(spec<0) = 0;
col = {'b','k','m'};
timeaxis = linspace(-1,1,size(Summary.psthReward{1},2));
idxT = timeaxis<=0.3 & timeaxis>=-0.3;

[maxRate,maxRateIdx] = max(Summary.psthReward{1},[],2);
[~,idxmax] = sort(maxRateIdx,'ascend'); 
avgRate = [];

for ii = 1:3
    ax1 = subplot(numrows,numcol,ii);
    norm = zscore(Summary.psthReward{ii},[],2);
    imagesc(timeaxis, 1:size(Summary.psthReward{ii},1),norm(idxmax,:))
    colormap(ax1,spec)
    caxis([-1 2])
    xlim([timeaxis(1) timeaxis(end)])

    a = Summary.psthReward{ii}(:,idxT);
    avgRate(:,ii) = mean(a,2);

    plotAvgStd(Summary.psthReward{ii},numrows,numcol,4,fig2,timeaxis',col{ii})
    xlim([timeaxis(1) timeaxis(end)])
end

Stats.toneCells = groupStats([{avgRate(:,1)},{avgRate(:,2)},{avgRate(:,3)}],[],'doPlot',false);

%% Show examples of a single session (IZ47_sess15) of overlays on the UMAP

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\IZ47_230626_sess15.position_behavior_speed_1_smooth_5_bin_0.1.mat';
A = -1.17; E = -31.94;
framelag = [-30 0 10];

sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
cd(sess)
file = dir('*.Tracking.Behavior.mat');
load(file.name)

Dec = findDecelerationPoints('plotfig',false);

for ii = 1:3
    if ii == 3
        decTS = Dec.ts(Dec.decType>=ii);
    else
        decTS = Dec.ts(Dec.decType==ii);
    end
    plotwithinTrialUMAP(umap_name,behav_file,A, E, decTS,tracking,fig2,numrows,numcol,4+ii,framelag)
end

%% Show averages of goal decoding jumps between the different conditions.
%AllDec = linkDecodingtoDeceleration('plotfig',false);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\DeliberationDecodingSummary.mat')
col = {'b','k','r','m'};

t = linspace(-1,1,201);
for dt = 1:4
    plotAvgStd(AllDec.pos{dt},numrows,numcol,9,fig2,t',col{dt})
    title('Position')  
    ylim([-7 3])
    if dt == 1
        xlim([-1 0.3])
    else
        xlim([-1 1])
    end

    plotAvgStd(AllDec.goal{dt},numrows,numcol,10,fig2,t',col{dt})
    title('Goal')  
    ylim([-1 1])
    if dt == 1
        xlim([-1 0.3])
    else
        xlim([-1 1])
    end
end

end

function plotwithinTrialUMAP(umap_name,behav_file,A, E, decTS,tracking,fighandle,numrows,numcol,plotloc,framelag)

TRIAL_TYPE = 0:5;
% load Umap result
Umap_results = readtable(umap_name);
Umap_results = table2array(Umap_results);
    
% load position direction and other information
load(behav_file);

plot_ind = [];
    
for tt = 1:length(TRIAL_TYPE)
    plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & probe_ds==0 & position_y_all>4 & speed_all'>3)];   
end

% Plot manifold in gray
ax1 = subplot(numrows,numcol,plotloc,'Parent',fighandle);
scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
view(A,E)
grid off;
axis off;
axis tight
hold on   

% Now plot the timepoints of interest
for dt = 1:length(decTS)
    [~,frameidx] = min(abs(tracking.timestamps-decTS(dt)));
    startTime = tracking.timestamps(frameidx+framelag(1));
    endTime = tracking.timestamps(frameidx+framelag(3));
    
    [~,startidx] = min(abs(timestamp_beh-startTime));
    [~,endidx] = min(abs(timestamp_beh-endTime));
    plot_ind_pos = zeros(1,length(position_y_all));
    plot_ind_pos(startidx:1:endidx) = 1;
    plot_ind_final = plot_ind_pos & position_y_all>1 & speed_all' >1;
    tsAxis = linspace(startTime,endTime,sum(plot_ind_final));
    scatter3(Umap_results(plot_ind_final,1),Umap_results(plot_ind_final,2),Umap_results(plot_ind_final,3),5,tsAxis','filled');    
    colormap(ax1,'magma');    
end
colorbar

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    rowsWithNaN = all(isnan(array), 2);
    numRowsWithoutNaN = sum(~rowsWithNaN);

    stdpsth = nanstd(array,1)./sqrt(numRowsWithoutNaN);
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    a = ~isnan(lArr);

    fill([xAxis(a); flipud(xAxis(a))],[lArr(a)'; flipud(uArr(a)')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis(a),meanpsth(a),'LineWidth',1,'Color',col);

end