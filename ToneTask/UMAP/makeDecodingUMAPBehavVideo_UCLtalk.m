function makeDecodingUMAPBehavVideo_UCLtalk

sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
cd(sess)
file = dir('*.thetaLFP.mat');
load(file.name)
file = dir('*.spikeData.cellInfo.mat');
load(file.name)
file = dir('*.Tracking.Behavior.mat');
load(file.name)
file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)

RdPu=cbrewer('seq', 'RdPu', 11);
Purples=cbrewer('seq', 'PuBu', 11);

decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\IZ47\Final\IZ47_230626_sess15\py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';

posterior_goal =  ncread(decodingPath,'x_position') ;
posterior_pos =  ncread(decodingPath,'y_position') ;
post_time =  ncread(decodingPath,'time') ;
post_pos =  ncread(decodingPath,'y_position_value') ;
post_goal =  ncread(decodingPath,'x_position_value') ;

cell1 = 86;%107;%45;
cell2 = 166;%54;%166;%275;%218;

trials = [8];%4 8

frontVid = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15\IZ47_230626_160943\front2023-06-26T16_11_06';
videoObjFront = VideoReader([frontVid '.avi']);   % get video

saveLoc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\';
saveName = 'ucl_IZ47_230626_sess15';

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\IZ47_230626_sess15.position_behavior_speed_1_smooth_5_bin_0.1.mat';

A = -1.17; E = -31.94;

numrows = 6;
numcol  = 10;

for tt = 1:length(trials) % Make separate videos for each trial

    writerObj = VideoWriter(strcat(saveLoc,saveName,'_trial',num2str(trials(tt))),'MPEG-4');
    writerObj.FrameRate = 5;
    open(writerObj);

    [~,startFrame] = min(abs(tracking.timestamps-behavTrials.timestamps(trials(tt),1)));
    %Only start the frame once the mouse starts moving    
    startFrame = startFrame+find(tracking.position.vy(startFrame:end)>10,1,'first');
    [~,endFrame] = min(abs(tracking.timestamps-behavTrials.timestamps(trials(tt),2)));
    tsWin = [tracking.timestamps(startFrame) tracking.timestamps(endFrame)];

    fig = figure;
    set(fig,'Position',[436 100 1377 600])
    set(fig,'Color','k')

    subplot(numrows,numcol,[1 2 3 4 11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44 51 52 53 54]) %Front video
    set(gca,'visible','off')
    
    subplot(numrows,numcol,[5 6 15 16 25 26]) % Position UMAP - static
    set(gca,'visible','off')

    subplot(numrows,numcol,[7 8 17 18 27 28]) % Goal UMAP - static
    set(gca,'visible','off')
    
    subplot(numrows,numcol,[9 10 19 20 29 30]) % Single trial UMAP
    set(gca,'visible','off')

    ax = subplot(numrows,numcol,35:40); % Bayesian posterior
    colormap(ax, 'magma');
    caxis([0 0.5])
    box off
    ax = gca;
    % Set the background color of the axes to black
    ax.Color = 'black';
    ax.XColor = 'none';
    ax.YColor = 'white';
    % Set the color of the tick marks to white
    ax.TickDir = 'out';  % optional: set tick direction
    ax.TickLength = [0.02 0.02];  % optional: set tick length   
    % Set the color of the axis labels to white
    ylabel('Position (cm)', 'Color', 'white');
    set(ax,'XLim',tsWin);
    set(ax,'YLim',[min(tracking.position.y) max(tracking.position.y)]);
    set(gca,'Ydir','reverse')
    hold on

    ax = subplot(numrows,numcol,45:50); % Goal posterior
    colormap(ax, 'magma');
    caxis([0 1.1])
    box off
    ax = gca;
    % Set the background color of the axes to black
    ax.Color = 'black';
    ax.XColor = 'none';
    ax.YColor = 'white';
    % Set the color of the tick marks to white
    ax.TickDir = 'out';  % optional: set tick direction
    ax.TickLength = [0.02 0.02];  % optional: set tick length   
    % Set the color of the axis labels to white
    %xlabel('Time (s)', 'Color', 'white');
    ylabel('Goal (port)', 'Color', 'white');
    set(ax,'XLim',tsWin);
    set(ax,'YLim',[0.5 6.5]);
    set(gca,'Ydir','reverse')
    hold on

    ax1 = subplot(numrows,numcol,55:60); % Position, spikes, speed
    box off
    ax = gca;
    % Set the background color of the axes to black
    ax.Color = 'black';    
    % Set the color of the axis lines and labels to white
    ax.XColor = 'white';
    ax.YColor = 'white';
    % Set the color of the tick marks to white
    ax.TickDir = 'out';  % optional: set tick direction
    ax.TickLength = [0.02 0.02];  % optional: set tick length   
    ylabel('Position (cm)', 'Color', 'white');
    set(ax1,'XLim',tsWin);
    set(ax1,'YLim',[min(tracking.position.y) 125]);
    line([tsWin(1) tsWin(2)],[18 18],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')   
    hold on
    line([tsWin(1) tsWin(2)],[35 35],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')        
    line([tsWin(1) tsWin(2)],[61 61],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
    line([tsWin(1) tsWin(2)],[82 82],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
    line([tsWin(1) tsWin(2)],[108 108],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
    line([tsWin(1) tsWin(2)],[122 122],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
    xlabel('Time (s)', 'Color', 'white');
    set(gca,'Ydir','reverse')
    h2 = animatedline('Color',[0.8 0.8 0.8],'LineWidth',1.5);

    for ii = startFrame:endFrame
        
        %% Current front camera
        subplot(numrows,numcol,[1 2 3 4 11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44 51 52 53 54]);
        frame = read(videoObjFront,ii);
        frame1 = frame(1:400,200:500,:);
        %frame1 = frame(50:280,280:430,:);
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

        %% Plot  manifold
        if ii == startFrame
            plotAvgManifold(umap_name,behav_file,numrows,numcol,[5 6 15 16 25 26],A,E,fig,1)
            plotAvgManifold(umap_name,behav_file,numrows,numcol,[7 8 17 18 27 28],A,E,fig,2)
        end
            
        videoManifold(umap_name,behav_file,numrows,numcol,[9 10 19 20 29 30],A,E,fig,[tracking.timestamps(startFrame) tracking.timestamps(ii)],ii)                

        %% Plot posterior
        ax1 = subplot(numrows,numcol,35:40);
        [~,tsfirst] = min(abs(post_time-tracking.timestamps(startFrame)));
        [~,tscur] = min(abs(post_time-tracking.timestamps(ii)));
        if ii>startFrame
            imagesc(post_time(tsfirst:tscur),post_pos,posterior_pos(:,tsfirst:tscur))
        end  
        line([tsWin(1) tsWin(2)],[9 9],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')   
        hold on
        line([tsWin(1) tsWin(2)],[35 35],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')        
        line([tsWin(1) tsWin(2)],[61 61],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[82 82],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[105 108],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[122 122],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')

        ax1 = subplot(numrows,numcol,45:50);
        if ii>startFrame
            imagesc(post_time(tsfirst:tscur),post_goal+1,posterior_goal(:,tsfirst:tscur))
        end        
        line([tsWin(1) tsWin(2)],[1.5 1.5],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        hold on
        line([tsWin(1) tsWin(2)],[2.5 2.5],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[3.5 3.5],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[4.5 4.5],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        line([tsWin(1) tsWin(2)],[5.5 5.5],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
        hold on;
        
        %% Plot spikes
        addpoints(h2,tracking.timestamps(ii),tracking.position.y(ii));
        subplot(numrows,numcol,55:60);
        set(gca,'YDir','reverse')
        if sum(ismember(spikeData.posIdx{cell1}, ii))>0
            ax1 = subplot(numrows,numcol,55:60);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.y(ii),30,'m',"filled")
        end
        if sum(ismember(spikeData.posIdx{cell2}, ii))>0
            ax1 = subplot(numrows,numcol,55:60);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.y(ii),30,'w',"filled")
        end


        M = getframe(fig);   
        writeVideo(writerObj, M);

    end   
   close(writerObj)
end

end

function plotAvgManifold(umap_name,behav_file,numrows,numcol,plotloc,A,E,figHandle,plotcol)
    
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

%% load Umap result
Umap_results = readtable(umap_name);
Umap_results = table2array(Umap_results);

%% load position direction and other information
load(behav_file);
TRIAL_TYPE = 0:5;
plot_ind = [];

for tt = 1:length(TRIAL_TYPE)
    plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==1 & probe_ds==0 & position_y_all>7 & speed_all'>2)];   
end
lick_plot = lick_loc_ds; 
lick_plot(isnan(lick_plot))=0; % deal with nan
lick_plot = lick_plot(plot_ind);

pos_plot = position_y_all;
pos_plot(isnan(pos_plot))=0; % deal with nan
pos_plot = pos_plot(plot_ind);

ax2 = subplot(numrows,numcol,plotloc,'Parent', figHandle);
hold on;
if plotcol == 1
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,pos_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,'jet');
    cb = colorbar('Color','w');
    cb.Label.String = 'Position';
else
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,lick_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,col);
    cb = colorbar('Color','w');
    cb.Label.String = 'Chosen port';
end
view(A,E)
grid off;
axis off;

end

function videoManifold(umap_name,behav_file,numrows,numcol,plotloc,A,E,figHandle,win,ii)

%% load Umap result
Umap_results = readtable(umap_name);
Umap_results = table2array(Umap_results);

%% load position direction and other information
load(behav_file);
TRIAL_TYPE = 0:5;
plot_ind = [];
gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
freqExp = log10(22000/1000);

for tt = 1:length(TRIAL_TYPE)
    plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==1 & probe_ds==0 & position_y_all>7 & speed_all'>1)];   
end

for tt = 1:length(position_y_all)
    if trial_type_ds(tt)<6
        freq = (position_y_all(tt)*gain(trial_type_ds(tt)+1))/122;
        tonepos_all(tt) = 1000*(10.^(freqExp*freq));
    else
        tonepos_all(tt) = 0;
    end
end

freq_plot = tonepos_all;
freq_plot(isnan(freq_plot))=0; % deal with nan


subplot(numrows,numcol,plotloc,'Parent', figHandle);
scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.03);

hold on
ax1 = subplot(numrows,numcol,plotloc,'Parent', figHandle);
[~,idxStart] = min(abs(timestamp_beh-win(1)));
[~,idxEnd] = min(abs(timestamp_beh-win(2)));
scatter3(Umap_results(idxStart:idxEnd,1),Umap_results(idxStart:idxEnd,2),Umap_results(idxStart:idxEnd,3),30,freq_plot(idxStart:idxEnd),'filled');  
colormap(ax1,viridis);
caxis([2000 23000])
view(A,E)
grid off;
axis off
cb = colorbar('Color','w');
cb.Label.String = 'Frequency';
cb.Label.Color = 'w'; % Ensure label color is set to white
set(gca,'ColorScale','log')

end