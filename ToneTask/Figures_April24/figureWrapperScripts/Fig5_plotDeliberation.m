function Fig5_plotDeliberation

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[463 2 1261 947]);

numrows = 6;
numcol = 8;

sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
cd(sess)
file = dir('*.Tracking.Behavior.mat');
load(file.name)
file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)
file = dir('*.spikeData.cellinfo.mat');
load(file.name)
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);

acc = gradient(tracking.position.v)./0.033;

cellIDs = [86 166];%50];% 54];
framelag = [-30 0 10];
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

%% Examples of two tone cells
for cc = 1:length(cellIDs)
    dataMat = [];
    dataMatTone = [];
    b = linspace(0,125,50);
    a = linspace(2000,22000,50);
    for kk = 1:6
        subplot(numrows, numcol, cc);
        hold on        
        plot(b,firingMaps.forward.rateMaps{cellIDs(cc)}{1+kk},'Color',col(kk,:),'LineWidth',1);
        xlim([0 122])
        set(gca,'xtick',[])
        box off    
        dataMat = [dataMat;firingMaps.forward.rateMaps{cellIDs(cc)}{1+kk}];
        atm = fillmissing(firingMaps.tone.rateMaps{cellIDs(cc)}{1+kk},'linear');
        dataMatTone = [dataMatTone;atm];  
    end

    BuPu=cbrewer('seq', 'BuPu', 11);
    ax2 = subplot(numrows, numcol, cc+8);
    imagesc(a,1,nanmean(dataMatTone,1))
    colormap(ax2, BuPu)
    box off
    set(gca,'ytick',[])
    xscale log
end

%% Panel A, stills of the approach to port 4, trial 104
col ={'r','c'};

decTS = 7539.328533;
[~,frameidx] = min(abs(tracking.timestamps-decTS));
extractCameraFrames(frameidx,fig2,numrows,numcol,3,framelag)

% plot speed and acceleration profile and superimpose the spikes
subplot(numrows,numcol,[11 12 13])
plot(tracking.timestamps(frameidx+framelag(1):frameidx+framelag(3)),acc(frameidx+framelag(1):frameidx+framelag(3)),'Color','k','LineWidth',1.5);
box off
axis tight
hold on
line([tracking.timestamps(frameidx) tracking.timestamps(frameidx)],[-40 30],'Color','k','LineWidth',1.5)
ylim([-40 30])
for cc = 1:length(cellIDs)
    a = ismember(spikeData.posIdx{cellIDs(cc)}, frameidx+framelag(1):frameidx+framelag(3));
    idxplot = spikeData.posIdx{cellIDs(cc)}(a);
    scatter(tracking.timestamps(idxplot ),acc(idxplot ),15,col{cc},'filled')   
end
line([behavTrials.timestamps(104,2) behavTrials.timestamps(104,2)],[-40 30])

% plot position decoding
plotDecoding(decTS,tracking,fig2,numrows,numcol,19,framelag,behavTrials)

% plot within trial umap
plotwithinTrialUMAP(decTS,tracking,fig2,numrows,numcol,35,framelag)

%% Panel B, stills of the approach to port 4 deliberation, trial 104
decTS = 7538.195166;
[~,frameidx] = min(abs(tracking.timestamps-decTS));
extractCameraFrames(frameidx,fig2,numrows,numcol,6,framelag)

% plot speed and acceleration profile and superimpose the spikes
subplot(numrows,numcol,[14 15 16])
plot(tracking.timestamps(frameidx+framelag(1):frameidx+framelag(3)),acc(frameidx+framelag(1):frameidx+framelag(3)),'Color','k','LineWidth',1.5);
box off
axis tight
hold on
line([tracking.timestamps(frameidx) tracking.timestamps(frameidx)],[-40 30],'Color','k','LineWidth',1.5)
ylim([-40 30])
for cc = 1:length(cellIDs)
    a = ismember(spikeData.posIdx{cellIDs(cc)}, frameidx+framelag(1):frameidx+framelag(3));
    idxplot = spikeData.posIdx{cellIDs(cc)}(a);
    scatter(tracking.timestamps(idxplot ),acc(idxplot ),15,col{cc},'filled')   
end

% plot position decoding
plotDecoding(decTS,tracking,fig2,numrows,numcol,22,framelag,behavTrials)

% plot within trial umap
plotwithinTrialUMAP(decTS,tracking,fig2,numrows,numcol,38,framelag)


%% Now, Panel C show all the deceleration points
combDec = combineDecelerationPoints('plotfig',false);

tAxis = linspace(-1,1,61);
col = {'b','k','m'};
for ii = 1:3
    if ii == 3
        plotAvgStd(combDec.velPSTH(combDec.decType>=ii,:),numrows,numcol,17,fig2,tAxis',col{ii})
        hold on
        title('Speed')    
    
        plotAvgStd(combDec.accPSTH(combDec.decType>=ii,:),numrows,numcol,18,fig2,tAxis',col{ii})  
        hold on
        title('Acceleration')
    else
        plotAvgStd(combDec.velPSTH(combDec.decType==ii,:),numrows,numcol,17,fig2,tAxis',col{ii})
        hold on
        title('Speed')    
    
        plotAvgStd(combDec.accPSTH(combDec.decType==ii,:),numrows,numcol,18,fig2,tAxis',col{ii})  
        hold on
        title('Acceleration')
    end
end

subplot(numrows,numcol,[25 26])
scatter(combDec.timetoLick1(combDec.decType==1), combDec.posY(combDec.decType==1), 5,'b.')
hold on
scatter(combDec.timetoLick1(combDec.decType==2), combDec.posY(combDec.decType==2), 5,'k.')
scatter(combDec.timetoLick1(combDec.decType>=3), combDec.posY(combDec.decType>=3), 5,'m.')
xlabel('time')
ylabel('Position')

subplot(numrows,numcol,[33 34 41 42])
scatter(combDec.posX(combDec.decType==1), combDec.posY(combDec.decType==1), 5,'b.')
hold on
scatter(combDec.posX(combDec.decType==2), combDec.posY(combDec.decType==2), 5,'k.')
scatter(combDec.posX(combDec.decType>=3), combDec.posY(combDec.decType>=3), 5,'m.')
xlabel('x position')
ylabel('y position')

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure5A_deliberationExample.png'));
saveas(gcf,strcat(expPath,'Figure5A_deliberationExample.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure5A_deliberationExample.fig'));

end

function plotwithinTrialUMAP(decTS,tracking,fighandle,numrows,numcol,plotloc,framelag)

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\IZ47_230626_sess15.position_behavior_speed_1_smooth_5_bin_0.1.mat';

A = -1.17; E = -31.94;

TRIAL_TYPE = [0:5];
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
ax1 = subplot(numrows,numcol,[plotloc plotloc+1 plotloc+2 plotloc+8 plotloc+9 plotloc+10],'Parent',fighandle);
scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
view(A,E)
grid off;
axis off;
axis tight
hold on   

% Now plot the timepoints of interest

[~,frameidx] = min(abs(tracking.timestamps-decTS));
startTime = tracking.timestamps(frameidx+framelag(1));
endTime = tracking.timestamps(frameidx+framelag(3));

[~,startidx] = min(abs(timestamp_beh-startTime));
[~,endidx] = min(abs(timestamp_beh-endTime));
plot_ind_pos = zeros(1,length(position_y_all));
plot_ind_pos(startidx:1:endidx) = 1;
plot_ind_final = plot_ind_pos & position_y_all>1 & speed_all' >1;
tsAxis = linspace(startTime,endTime,sum(plot_ind_final));
scatter3(Umap_results(plot_ind_final,1),Umap_results(plot_ind_final,2),Umap_results(plot_ind_final,3),15,tsAxis','filled');    
colormap(ax1,'magma');
colorbar

end

function plotDecoding(decTS,tracking,fighandle,numrows,numcol,plotloc,framelag, behavTrials)

RdPu=cbrewer('seq', 'RdPu', 11);
Purples=cbrewer('seq', 'PuBu', 11);

decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\IZ47\Final\IZ47_230626_sess15\py_data\theta_decoding_lickLoc_y';
file_name = 'up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';

posterior_goal =  ncread(strcat(decodingPath,'\',file_name),'x_position') ;
posterior_pos =  ncread(strcat(decodingPath,'\',file_name),'y_position') ;
post_time =  ncread(strcat(decodingPath,'\',file_name),'time') ;
post_pos =  ncread(strcat(decodingPath,'\',file_name),'y_position_value') ;
post_goal =  ncread(strcat(decodingPath,'\',file_name),'x_position_value') ;

file_name = 'change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';
load(strcat(decodingPath,'\',file_name));

[~,idxstartTrue] = min(abs(post_time-behavTrials.timestamps(104,1)));   

[~,frameidx] = min(abs(tracking.timestamps-decTS));
startTime = tracking.timestamps(frameidx+framelag(1));
endTime = tracking.timestamps(frameidx+framelag(3));

[~,idxstart] = min(abs(post_time-startTime));      
[~,idxend] = min(abs(post_time-endTime));

ax1 = subplot(numrows,numcol,[plotloc plotloc+1 plotloc+2],'Parent',fighandle);
imagesc(post_time(idxstart:idxend),post_pos,posterior_pos(:,idxstart:idxend))
set(gca,'YDir','normal')
colormap(ax1,RdPu)
box off
clim([0 0.4])
xlim([startTime endTime])
hold on
plot(tracking.timestamps(frameidx+framelag(1):frameidx+framelag(3)),tracking.position.y(frameidx+framelag(1):frameidx+framelag(3)),'Color','k')

ax1 = subplot(numrows,numcol,[plotloc+8 plotloc+9 plotloc+10],'Parent',fighandle);
imagesc(post_time(idxstart:idxend),post_goal,posterior_goal(:,idxstart:idxend))
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on

curChanges = change_point{trial==(103)};
for cr = 1:length(curChanges)
    line([post_time(idxstartTrue+curChanges(cr)-1) post_time(idxstartTrue+curChanges(cr)-1)],[0 6],'Color','r','LineWidth',1.5);  
end
xlim([startTime endTime])

end

function extractCameraFrames(ii,fighandle,numrow,numcol,plotloc,idxlag)

frontVid = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15\IZ47_230626_160943\front2023-06-26T16_11_06';
videoObjFront = VideoReader([frontVid '.avi']);   % get video

for ll = 1:length(idxlag)
    subplot(numrow,numcol,plotloc+(ll-1),'Parent',fighandle)
    frame = read(videoObjFront,ii+idxlag(ll));
    frame1 = frame(1:400,200:500,:);
    % Mask the edges
    vertices = [40,400;...
                145, 13;...
                145, 0;...
                183,0;...
                184,13;...
                272, 400];

    mask = poly2mask(vertices(:, 1), vertices(:, 2), size(frame1, 1), size(frame1, 2));
    mask = ~mask;
    frame1(repmat(mask, [1, 1, size(frame1, 3)])) = 0;
    imagesc(frame1)   
    axis off
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