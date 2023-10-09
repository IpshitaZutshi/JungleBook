function plotTrajectorySpatial_SupFigure6

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);
%rows = 22, columns = 6
numrows = 18;
numcol = 9;


%% Panel F: 2-D maps of the same neurons from A-C
% Cell 1
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 7;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);

subplot(numrows,numcol,[1 2])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off
title(strcat(num2str(maxRate),' Hz'))

subplot(numrows,numcol,[1+1*numcol 2+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off

subplot(numrows,numcol,[1+2*numcol 2+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off

corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

% Cell 2
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 75;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);

subplot(numrows,numcol,[3 4])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
caxis([0 maxRate-5])
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
set(gca,'YDir','normal')
axis off
title(strcat(num2str(maxRate),' Hz'))

subplot(numrows,numcol,[3+1*numcol 4+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off

subplot(numrows,numcol,[3+2*numcol 4+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

% Cell 3
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 114;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);
axis off

subplot(numrows,numcol,[5 6])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
title(strcat(num2str(maxRate),' Hz'))
axis off

subplot(numrows,numcol,[5+1*numcol 6+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
axis off

subplot(numrows,numcol,[5+2*numcol 6+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
axis off
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

%% Correlation of 2D maps
Summary = compileMice2DPlaceFields;
col = [52/243 52/243 52/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243;...
    180/243 180/243 180/243];
    
subplot(numrows,numcol,[7 8 7+1*numcol 8+1*numcol 7+2*numcol 8+2*numcol])
SupFig6Stats.Maps2D = groupStats({Summary.linCorr,Summary.toneNoToneCorr,Summary.tonelinEndCorr,Summary.linlinEndCorr},[],...
    'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_trajectTuning.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_trajectTuning.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_trajectTuning.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_trajectTuning.mat'),'SupFig6Stats'); 

end