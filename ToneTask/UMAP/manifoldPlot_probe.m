function manifoldPlot_probe(varargin)

p = inputParser;
addParameter(p,'umap_path',pwd);
addParameter(p,'umap_name','behavior_speed_1_smooth_5_0.1');
addParameter(p,'behav_file',[]);
addParameter(p,'A',-7);
addParameter(p,'E',40);
addParameter(p,'numrow',1);
addParameter(p,'addFreq',false);
addParameter(p,'numcol',1);
addParameter(p,'rowloc',1);
addParameter(p,'colloc',1);
addParameter(p,'probe',false);
addParameter(p,'col','jet');

addParameter(p,'dim1',1);
addParameter(p,'dim2',2);
addParameter(p,'dim3',3);
addParameter(p,'TRIAL_TYPE',[5 6 8]);
addParameter(p,'figHandle',[]);

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
probe = p.Results.probe;
col = p.Results.col;
dim1 = p.Results.dim1;
dim2 = p.Results.dim2;
dim3 = p.Results.dim3;
TRIAL_TYPE = p.Results.TRIAL_TYPE;
figHandle = p.Results.figHandle;

%% load Umap result
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

gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
freqExp = log10(22000/1000);

plot_ind = [];
for tt = 1:length(TRIAL_TYPE)
    if TRIAL_TYPE(tt)<6 && exist('correct_ds')
        if exist('probe_ds')
            %plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==1 & probe_ds==probe)];
            plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & probe_ds==probe)];
        else
            %plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==1)]; 
            plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt))]; 
        end
    else
        plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt))];
    end
%    plot_ind =  [plot_ind,find(trial_type_ds==TRIAL_TYPE(tt))]; 
end

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

tonepos_all(tonepos_all>25000) = 0;
freq_plot = tonepos_all;
freq_plot(isnan(freq_plot))=0; % deal with nan
freq_plot = freq_plot(plot_ind);

arm_plot = trial_type_ds;
arm_plot(isnan(arm_plot))=0; % deal with nan
arm_plot = arm_plot(plot_ind);

lick_plot = lick_loc_ds; 
lick_plot(isnan(lick_plot))=0; % deal with nan
lick_plot = lick_plot(plot_ind);

% color by position
BuPu=cbrewer('div', 'PiYG', 500);
ax1 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc,'Parent', figHandle);
set(figHandle,'color','w');
%scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
hold on;
scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),8,pos_plot,'filled');
%colorbar
colormap(ax1,'jet');
view(A,E)
grid off;
title('Neural Manifold (position)')
xlabel(['Dim' num2str(dim1)]);
ylabel(['Dim' num2str(dim2)]);
zlabel(['Dim' num2str(dim3)]);
axis off;
%colorbar

% color by frequency
if length(TRIAL_TYPE)==6 && addFreq
    BuPu=viridis;   
    ax1 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc+2,'Parent', figHandle);
    set(figHandle,'color','w');
    %scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
    hold on;
    scatter3(Umap_results(plot_ind,dim1),Umap_results(plot_ind,dim2),Umap_results(plot_ind,dim3),8,freq_plot,'filled');
    %colorbar
    colormap(ax1,BuPu);
    view(A,E)
    grid off;
    title('Neural Manifold (position)')
    xlabel(['Dim' num2str(dim1)]);
    ylabel(['Dim' num2str(dim2)]);
    zlabel(['Dim' num2str(dim3)]);
    axis off;
   % colorbar
end
% color by trialType
ax2 = subplot(numrow,numcol,numcol*(rowloc-1)+colloc+1,'Parent', figHandle);
%scatter3(Umap_results(plot_ind_rest,dim1),Umap_results(plot_ind_rest,dim2),Umap_results(plot_ind_rest,dim3),5,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.1);
hold on;
for ii = 1:6
    res = Umap_results(plot_ind,:);
    lick_id = lick_plot==(ii-1);
    scatter3(res(lick_id,dim1),res(lick_id,dim2),res(lick_id,dim3),8,col(ii,:),'filled');
end
view(A,E)
grid off;
title('Neural Manifold (position)')
xlabel(['Dim' num2str(dim1)]);
ylabel(['Dim' num2str(dim2)]);
zlabel(['Dim' num2str(dim3)]);
axis off;
%colorbar

Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);
end

