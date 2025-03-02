function manifoldPlot(varargin)

p = inputParser;
addParameter(p,'umap_path',pwd);
addParameter(p,'umap_name','behavior_speed1'); % 'behavior_speed_1_smooth_1_bin_0.1'
addParameter(p,'behav_file',[]);
addParameter(p,'A',-7);
addParameter(p,'E',40);
addParameter(p,'numrow',1);
addParameter(p,'addFreq',false);
addParameter(p,'addPosPlot',true);
addParameter(p,'numcol',1);
addParameter(p,'rowloc',1);
addParameter(p,'colloc',1);
addParameter(p,'col','jet');
addParameter(p,'singleTrial',false);
addParameter(p,'tsWin',[]);
addParameter(p,'dotSize',8);
addParameter(p,'error',1);
addParameter(p,'dim1',1);
addParameter(p,'dim2',2);
addParameter(p,'dim3',3);
addParameter(p,'TRIAL_TYPE',[5 6 8]);
addParameter(p,'figHandle',[]);
addParameter(p,'poscutOff',0);
addParameter(p,'speedThresh',1)
addParameter(p,'plotcolorbar',false)

parse(p,varargin{:});
umap_path = p.Results.umap_path;
umap_name = p.Results.umap_name;
behav_file = p.Results.behav_file;
A = p.Results.A;
E = p.Results.E;
addFreq = p.Results.addFreq;
numrow = p.Results.numrow;
numcol = p.Results.numcol;
rowloc = p.Results.rowloc;
colloc = p.Results.colloc;
col = p.Results.col;
dim1 = p.Results.dim1;
dim2 = p.Results.dim2;
dim3 = p.Results.dim3;
TRIAL_TYPE = p.Results.TRIAL_TYPE;
figHandle = p.Results.figHandle;
singleTrial = p.Results.singleTrial;
tsWin = p.Results.tsWin;
dotSize = p.Results.dotSize;
error = p.Results.error;
addPosPlot = p.Results.addPosPlot;
poscutOff = p.Results.poscutOff;
speedThresh = p.Results.speedThresh;
plotcolorbar = p.Results.plotcolorbar;

%% load Umap result
%umap_path = dir('Umap*');
%Umap_results = readtable(umap_path.name);
Umap_results = readtable([umap_path, '\Umap_',umap_name,'.csv']);
Umap_results = table2array(Umap_results);

%% load position direction and other information
load(behav_file);

% Correct the second block of linear trials
a = find(trial_type_ds==6);
b = diff(a);
[~, idx] = max(b);
% indices of second half of lin trials
if b(idx)>250 % If there is a second unlabelled block
    c = a(idx+1:end);
    trial_type_ds(c) = 8;
end

lick_loc_ds(trial_type_ds==6) = 6;
lick_loc_ds(trial_type_ds==7) = 7;
lick_loc_ds(trial_type_ds==8) = 8;

if length(TRIAL_TYPE)==9
    lick_loc_ds(lick_loc_ds~=7) = 1;
    TRIAL_TYPE = [1 7];
end

gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
freqExp = log10(22000/1000);

plot_ind = [];
for tt = 1:length(TRIAL_TYPE)
    if TRIAL_TYPE(tt)<6 && exist('correct_ds')
        if exist('probe_ds')
            plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==error & probe_ds==0 & position_y_all>poscutOff & speed_all' >speedThresh)];
        else
            plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==error & position_y_all>poscutOff & speed_all' >speedThresh)]; 
        end
    else
        plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & position_y_all>poscutOff & speed_all' >speedThresh)];
    end
%    plot_ind =  [plot_ind,find(trial_type_ds==TRIAL_TYPE(tt))]; 
end
% plot_ind = plot_ind(plot_ind<9700);
% a(1:9700)= 1;
% a(plot_ind) = 0;
% plot_ind_rest = find(a);
%Take a subset
%plot_ind = plot_ind(1:10:length(plot_ind));
for tt = 1:length(position_y_all)
    if trial_type_ds(tt)<6
        freq = (position_y_all(tt)*gain(trial_type_ds(tt)+1))/122;
        tonepos_all(tt) = 1000*(10.^(freqExp*freq));
    else
        tonepos_all(tt) = 0;
    end
end
pos_plot = position_y_all;
pos_plot(isnan(pos_plot))=0; % deal with nan
pos_plot = pos_plot(plot_ind);

%tonepos_all(tonepos_all>25000) = 0;
freq_plot = tonepos_all;
freq_plot(isnan(freq_plot))=0; % deal with nan
freq_plot = freq_plot(plot_ind);

arm_plot = trial_type_ds;
arm_plot(isnan(arm_plot))=0; % deal with nan
arm_plot = arm_plot(plot_ind);

lick_plot = lick_loc_ds; 
lick_plot(isnan(lick_plot))=0; % deal with nan
lick_plot = lick_plot(plot_ind);
    
if ~singleTrial

    if addPosPlot
        % color by position
        BuPu=cbrewer('div', 'PiYG', 500);
        ax1 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc,'Parent', figHandle);
        set(figHandle,'color','w');
        %scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
        hold on;
        scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),dotSize,pos_plot,'filled');
        %colorbar
        colormap(ax1,'jet');
        view(A,E)
        grid off;
        axis tight
        axis off;
        if plotcolorbar
            colorbar
        end
    end
    
    % color by frequency
    if length(TRIAL_TYPE)==6 && addFreq
        BuPu=viridis;   
        ax1 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc+2,'Parent', figHandle);
        set(figHandle,'color','w');
        %scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
        hold on;
        scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),dotSize,freq_plot,'filled');
        %colorbar
        colormap(ax1,BuPu);
        view(A,E)
        grid off;
        axis off;
        axis tight
        if plotcolorbar
            colorbar
        end
        caxis([2000 23000])
        set(gca,'colorscale','log')
    end

    % color by trialType
    if addPosPlot
        ax2 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc+1,'Parent', figHandle);
    else
        ax2 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc,'Parent', figHandle);
    end
    %scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
    hold on;
    scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),dotSize,lick_plot,'filled');
    %colorbar
    colormap(ax2,col);
    view(A,E)
    grid off;
    axis off;
    axis tight
    if plotcolorbar
        colorbar
    end
    
    if addPosPlot
        Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
    end

else
%% Plot the entire manifold in gray
    ax1 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc,'Parent', figHandle);
    set(figHandle,'color','w');
    hold on;
    scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),2,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
    %colorbar
    view(A,E)
    grid off;
    axis off;
    axis tight
    %colorbar

%% Now plot the trial of interest
    freq_plot = tonepos_all;
    freq_plot(isnan(freq_plot))=0; % deal with nan

    [~,startidx] = min(abs(timestamp_beh-tsWin(1)));
    [~,endidx] = min(abs(timestamp_beh-tsWin(2)));
    plot_ind_pos = zeros(1,length(freq_plot));
    plot_ind_pos(startidx:1:endidx) = 1;
    plot_ind = plot_ind_pos & position_y_all>poscutOff & speed_all' >speedThresh;
    tsAxis = linspace(tsWin(1),tsWin(2),sum(plot_ind));
    %scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),dotSize,freq_plot(plot_ind),'filled');    
    scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),dotSize,tsAxis','filled');    
    %colormap jet
    %colormap(ax1,'viridis');
    colormap(ax1,'magma');
    colorbar
    % caxis([2000 23000])
    % set(gca,'colorscale','log')
end
end

