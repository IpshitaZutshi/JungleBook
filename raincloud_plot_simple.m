function h = raincloud_plot_simple(X,cl)

[a,b] = ksdensity(X);

wdth = 0.8; % width of boxplot
% TODO, should probably be some percentage of max.height of kernel density plot

% density plot
% h{1} = area(a,b); hold on
h{1} = plot(a,b); hold on
%set(h{1}, 'FaceColor', cl);
set(h{1}, 'Color', cl);
%set(h{1}, 'EdgeColor', [0.1 0.1 0.1]);
set(h{1}, 'LineWidth', 2);

% make some space under the density plot for the boxplot
xl = get(gca,'XLim');
set(gca,'XLim',[-2 xl(2)]);

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
Y = quantile(X,[0.25 0.75 0.5 0.02 0.98]);

% 'box' of 'boxplot'
h{2} = rectangle('Position',[-1-(wdth*0.5) Y(1) wdth Y(2)-Y(1)]);
set(h{2},'EdgeColor',cl)
set(h{2},'LineWidth',2);
% could also set 'FaceColor' here as Micah does, but I prefer without

% mean line
h{3} = line([-1.2 -0.8],[Y(3) Y(3)],'col','k','LineWidth',2);

% whiskers
h{4} = line([-1 -1],[Y(2) Y(5)],'col','k','LineWidth',2);
h{5} = line([-1 -1],[Y(1) Y(4)],'col','k','LineWidth',2);

% raindrops
h{3} = scatter(jit - 1,X);
h{3}.SizeData = 5;
h{3}.MarkerFaceColor = cl;
h{3}.MarkerEdgeColor = 'none';
