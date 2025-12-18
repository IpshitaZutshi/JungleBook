function manifoldPlot_stim_overlay(varargin)

p = inputParser;
addParameter(p,'umap_path',pwd);
addParameter(p,'umap_name','behavior_speed_1_smooth_5_0.1');
addParameter(p,'behav_file',[]);
addParameter(p,'A',-7);
addParameter(p,'E',40);
addParameter(p,'poscutOff',0);
addParameter(p,'speedThresh',1); %#ok<NASGU>  % kept for compatibility

addParameter(p,'dim1',1);
addParameter(p,'dim2',2);
addParameter(p,'dim3',3);
addParameter(p,'TRIAL_TYPE',[0:5]);
addParameter(p,'figHandle',[]);
addParameter(p,'addFreq',false); %#ok<NASGU>  % kept for compatibility

% optional styling
addParameter(p,'ptSize',8);
addParameter(p,'alphaBase',1.0);
addParameter(p,'alphaStim',1.0);
addParameter(p,'col','jet');

parse(p,varargin{:});
umap_path  = p.Results.umap_path;
umap_name  = p.Results.umap_name;
behav_file = p.Results.behav_file;
A          = p.Results.A;
E          = p.Results.E;
dim1       = p.Results.dim1;
dim2       = p.Results.dim2;
dim3       = p.Results.dim3;
TRIAL_TYPE = p.Results.TRIAL_TYPE;
figHandle  = p.Results.figHandle;
poscutOff  = p.Results.poscutOff;

ptSize     = p.Results.ptSize;
alphaBase  = p.Results.alphaBase;
alphaStim  = p.Results.alphaStim;
col = p.Results.col;

%% load Umap result
Umap_results = readtable(fullfile(umap_path, ['Umap_', umap_name, '.csv']));
Umap_results = table2array(Umap_results);

%% load position direction and other information
load(behav_file);

% ----- Correct the second block of linear trials (your original logic) -----
a = find(trial_type_ds==6);
b = diff(a);
[~, idx] = max(b);
if ~isempty(b) && b(idx)>250
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

%% build plot_ind (your original selection logic, cleaned)
plot_ind = [];

for tt = 1:length(TRIAL_TYPE)
    if TRIAL_TYPE(tt) < 6
        plot_ind = [plot_ind, find(lick_loc_ds==TRIAL_TYPE(tt) & position_y_all>poscutOff & trial_num_ds>2)]; %#ok<AGROW>
    else
        plot_ind = [plot_ind, find(lick_loc_ds==TRIAL_TYPE(tt))]; %#ok<AGROW>
    end
end

% keep unique + stable order
plot_ind = unique(plot_ind, 'stable');

%% prepare variables aligned to plot_ind
pos_plot = position_y_all;
pos_plot(isnan(pos_plot)) = 0;
pos_plot = pos_plot(plot_ind);

lick_plot = lick_loc_ds;
lick_plot(isnan(lick_plot)) = 0;
lick_plot = lick_plot(plot_ind);

stim_plot = stim_ds(plot_ind);   % 0 baseline, 1 stim
correct_plot = correct_ds(plot_ind);          % 1 correct, 0 incorrect

res = Umap_results(plot_ind,:);

%% figure setup
if isempty(figHandle) || ~ishandle(figHandle)
    figHandle = figure;
else
    figure(figHandle);
end
set(figHandle, 'color', 'w');

ax1 = subplot(2,3,1,'Parent',figHandle);
ax2 = subplot(2,3,2,'Parent',figHandle);
ax3 = subplot(2,3,3,'Parent',figHandle);
ax4 = subplot(2,3,4,'Parent',figHandle);
ax5 = subplot(2,3,5,'Parent',figHandle);

%% ---------------- 1) color by position ----------------
axes(ax1); %#ok<LAXES>
scatter3(res(:,dim1), res(:,dim2), res(:,dim3), ptSize, pos_plot, 'filled');
colormap(ax1, 'jet');
view(A,E); grid off; axis off;
title('Neural Manifold (position)');
xlabel(['Dim' num2str(dim1)]); ylabel(['Dim' num2str(dim2)]); zlabel(['Dim' num2str(dim3)]);

%% ---------------- 2) color by lick location ----------------
axes(ax2); %#ok<LAXES>
hold on;

lick_ids = unique(lick_plot);
lick_ids = lick_ids(:)';  % row

for ii = 1:numel(lick_ids)
    id = lick_ids(ii);
    idx = (lick_plot == id);
    scatter3(res(idx,dim1), res(idx,dim2), res(idx,dim3), ptSize, col(ii,:), 'filled');
end

view(A,E); grid off; axis off;
title('Neural Manifold (lick location)');
xlabel(['Dim' num2str(dim1)]); ylabel(['Dim' num2str(dim2)]); zlabel(['Dim' num2str(dim3)]);

%% ---------------- 3) baseline gray vs stim blue ----------------
axes(ax3); %#ok<LAXES>
hold on;

base_idx = (stim_plot == 0);
stim_idx = (stim_plot == 1);

% baseline (gray, transparent)
scatter3(res(base_idx,dim1), res(base_idx,dim2), res(base_idx,dim3), ...
    ptSize, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', alphaBase);

% stim (blue, opaque)
scatter3(res(stim_idx,dim1), res(stim_idx,dim2), res(stim_idx,dim3), ...
    ptSize, [0.2 0.4 1.0], 'filled', 'MarkerFaceAlpha', alphaStim);

view(A,E); grid off; axis off;
title('Baseline (gray) vs Stim (blue)');
xlabel(['Dim' num2str(dim1)]); ylabel(['Dim' num2str(dim2)]); zlabel(['Dim' num2str(dim3)]);

%% ---------------- 4) nostim: correct (gray) vs incorrect (magenta) ----------------
axes(ax4); %#ok<LAXES>
hold on;

nostim_idx = (stim_plot == 0);
nostim_corr = nostim_idx &  correct_plot;
nostim_err  = nostim_idx & ~correct_plot;

% correct (gray)
scatter3(res(nostim_corr,dim1), res(nostim_corr,dim2), res(nostim_corr,dim3), ...
    ptSize, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.8);

% incorrect (magenta)
scatter3(res(nostim_err,dim1), res(nostim_err,dim2), res(nostim_err,dim3), ...
    ptSize, [1.0 0.0 1.0], 'filled', 'MarkerFaceAlpha', 0.8);

view(A,E); grid off; axis off;
title('No-stim: correct (gray) vs incorrect (magenta)');
xlabel(['Dim' num2str(dim1)]); ylabel(['Dim' num2str(dim2)]); zlabel(['Dim' num2str(dim3)]);

%% ---------------- 5) stim: correct (blue) vs incorrect (magenta) ----------------
axes(ax5); %#ok<LAXES>
hold on;

stim_idx2 = (stim_plot == 1);
stim_corr = stim_idx2 &  correct_plot;
stim_err  = stim_idx2 & ~correct_plot;

% correct (blue)
scatter3(res(stim_corr,dim1), res(stim_corr,dim2), res(stim_corr,dim3), ...
    ptSize, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.8);

% incorrect (magenta)
scatter3(res(stim_err,dim1), res(stim_err,dim2), res(stim_err,dim3), ...
    ptSize, [1.0 0.0 1.0], 'filled', 'MarkerFaceAlpha', 0.8);

view(A,E); grid off; axis off;
title('Stim: correct (blue) vs incorrect (magenta)');
xlabel(['Dim' num2str(dim1)]); ylabel(['Dim' num2str(dim2)]); zlabel(['Dim' num2str(dim3)]);


%% link camera so rotating one rotates all
Link = linkprop([ax1 ax2 ax3 ax4 ax5], ...
    {'CameraUpVector','CameraPosition','CameraTarget','XLim','YLim','ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

end
