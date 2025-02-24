function manifoldPlot_Cue

Umap_results = readtable('Umap_CA1.csv');
Umap_results_CA1 = table2array(Umap_results);

Umap_results = readtable('Umap_V1.csv');
Umap_results_V1 = table2array(Umap_results);

load('Umap_behavior.mat');

TRIAL_TYPE = [1];

plot_ind = [];

for tt = 1:length(TRIAL_TYPE)
   plot_ind =  [plot_ind,find(direction_ds==TRIAL_TYPE(tt))]; 
end

posy_plot = pos_y_ds;
posy_plot = posy_plot(plot_ind);

posx_plot = pos_x_ds;
posx_plot = posx_plot(plot_ind);

cue_plot = cue_ds;
cue_plot = cue_plot(plot_ind);

choice_plot = choice_ds;
choice_plot = choice_plot(plot_ind);

col = [1 0 1;154/255 205/255 50/255];
fig2 = figure;
set(gcf,'Color','w')
ax1 = subplot(2,2,1);
scatter3(Umap_results_CA1(plot_ind,1),Umap_results_CA1(plot_ind,2),Umap_results_CA1(plot_ind,3),5,posy_plot)
grid off
axis off
title('Y position')
view(82.7387,90)
colormap(ax1,"jet")

ax2 = subplot(2,2,2);
scatter3(Umap_results_CA1(plot_ind,1),Umap_results_CA1(plot_ind,2),Umap_results_CA1(plot_ind,3),5,posx_plot)
clim([85 105])
grid off
axis off
title('X position')
view(82.7387,90)
colormap(ax2,"jet")

ax3 = subplot(2,2,3);
scatter3(Umap_results_CA1(plot_ind,1),Umap_results_CA1(plot_ind,2),Umap_results_CA1(plot_ind,3),5,cue_plot)
grid off
axis off
title('Cue')
view(82.7387,90)
colormap(ax3,col)

ax4 = subplot(2,2,4);
scatter3(Umap_results_CA1(plot_ind,1),Umap_results_CA1(plot_ind,2),Umap_results_CA1(plot_ind,3),5,choice_plot)
grid off
axis off
title('Choice')
view(82.7387,90)
colormap(ax4,col)

Link = linkprop([ax1, ax2, ax3,ax4],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

fig3 = figure;
set(gcf,'Color','w')
ax5 = subplot(2,2,1);
scatter3(Umap_results_V1(plot_ind,1),Umap_results_V1(plot_ind,2),Umap_results_V1(plot_ind,3),5,posy_plot)
grid off
axis off
title('Y position')
colormap(ax5,"jet")

ax6 = subplot(2,2,2);
scatter3(Umap_results_V1(plot_ind,1),Umap_results_V1(plot_ind,2),Umap_results_V1(plot_ind,3),5,posx_plot)
grid off
axis off
title('X position')
colormap(ax6,"jet")

ax7 = subplot(2,2,3);
scatter3(Umap_results_V1(plot_ind,1),Umap_results_V1(plot_ind,2),Umap_results_V1(plot_ind,3),5,cue_plot)
grid off
axis off
colormap(ax7,col)
title('Cue')

ax8=subplot(2,2,4);
scatter3(Umap_results_V1(plot_ind,1),Umap_results_V1(plot_ind,2),Umap_results_V1(plot_ind,3),5,choice_plot)
grid off
axis off
colormap(ax8,col)
title('Choice')

Link = linkprop([ax5,ax6,ax7,ax8],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

figure
set(gcf,'Color','w')
a = subplot(2,2,1);
scatter(posx_plot,posy_plot,5,posy_plot)
colormap(a,"jet")
colorbar
title('Y position')
b = subplot(2,2,2);
scatter(posx_plot,posy_plot,5,posx_plot)
title('X position')
colorbar
colormap(b,"jet")
c = subplot(2,2,3);
scatter(posx_plot,posy_plot,5,cue_plot)
title('Cue')
colorbar
colormap(c,col)
d = subplot(2,2,4);
scatter(posx_plot,posy_plot,5,choice_plot)
title('Choice')
colorbar
colormap(d,col)
% 
% subplot(1,3,2)
% scatter3(Umap_results_V1(:,1),Umap_results_V1(:,2),Umap_results_V1(:,3),'.')
end