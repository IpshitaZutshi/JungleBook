function makeDecodingUMAPBehavVideo_BerkeleyTalk

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


trials = [104];

frontVid = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15\IZ47_230626_160943\front2023-06-26T16_11_06';
videoObjFront = VideoReader([frontVid '.avi']);   % get video

saveLoc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\';
saveName = 'berkeley_IZ47_230626_sess15';

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\IZ47_230626_sess15.position_behavior_speed_1_smooth_5_bin_0.1.mat';

A = -1.17; E = -31.94;

numrows = 5;
numcol  = 3;

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
    set(fig,'Position',[723 51 956 923])
    set(fig,'Color','k')

    subplot(numrows,numcol,[1 2 4 5 7 8 10 11]) %Front video
    set(gca,'visible','off')

    subplot(numrows,numcol,[3 6]) % Goal UMAP - static
    set(gca,'visible','off')
    
    subplot(numrows,numcol,[9 12]) % Single trial UMAP
    set(gca,'visible','off')

    for ii = startFrame:endFrame
        
        %% Current front camera
        subplot(numrows,numcol,[1 2 4 5 7 8 10 11]);
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
            plotAvgManifold(umap_name,behav_file,numrows,numcol,[3 6],A,E,fig,2)
        end
            
        videoManifold(umap_name,behav_file,numrows,numcol,[9 12],A,E,fig,[tracking.timestamps(startFrame) tracking.timestamps(ii)],ii)                

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
lick_plot = lick_plot(plot_ind)+1;

pos_plot = position_y_all;
pos_plot(isnan(pos_plot))=0; % deal with nan
pos_plot = pos_plot(plot_ind);

ax2 = subplot(numrows,numcol,plotloc,'Parent', figHandle);
hold on;
if plotcol == 1
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,pos_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,'jet');
    cb = colorbar('Color','w','position',[0.65 0.69 0.01 0.14]);
    cb.Label.String = 'Position';
else
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,lick_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,col);
    cb = colorbar('Color','w','position',[0.65 0.71 0.01 0.14]);
    cb.Label.String = 'Position';
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
if ii == 64964
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.03);
end

hold on
ax1 = subplot(numrows,numcol,plotloc,'Parent', figHandle);
[~,idxStart] = min(abs(timestamp_beh-win(1)));
[~,idxEnd] = min(abs(timestamp_beh-win(2)));

plot_ind_pos = zeros(1,length(freq_plot));
plot_ind_pos(idxStart:1:idxEnd) = 1;
plot_ind = plot_ind_pos & position_y_all>7 & speed_all' >1;
tsAxis = linspace(win(1),win(2),sum(plot_ind));

scatter3(Umap_results(idxStart:idxEnd,1),Umap_results(idxStart:idxEnd,2),Umap_results(idxStart:idxEnd,3),30,tsAxis','filled');  
colormap(ax1,viridis);
caxis([win(1) 7539.3])
view(A,E)
grid off;
axis off
cb = colorbar('Color','w','position',[0.65 0.38 0.01 0.14]);
cb.Label.String = 'Time';
cb.Label.Color = 'w'; % Ensure label color is set to white

end