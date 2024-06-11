function FigS4_ThetaSequences

%% Plot three examples of trials with place and tone cells interspersed with theta oscillations
%plotThetaSequences

%% Plot examples of pairwise CCGs
fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[1950 40 900 950]);

numrows = 11;
numcol = 6;

%% Plot examples of cell sequences that combine place and "tone cells"
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotPosPhaseCCG(5,83,103,numrows, numcol, 1, fig2,40)
plotPosPhaseCCG(134,186,[],numrows, numcol, 4, fig2,30)
%plotPosPhaseCCG(134,64,[],numrows, numcol, 7, fig2,30)

%% Tone cells only
sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
cd(sess)
plotPosPhaseCCG(86,166,[],numrows, numcol, 7, fig2,40)

compiledthetaComp = compileThetaCompression;
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\thetaCompression.mat')

%% Theta sequences in the spatial domain
subplot(numrows,numcol,[9*numcol+1 9*numcol+2 10*numcol+1 10*numcol+2])
[R,P] = corrcoef(compiledthetaComp.place_offset,compiledthetaComp.ccgs_time_offset,'rows','complete');
scatter(compiledthetaComp.place_offset,compiledthetaComp.ccgs_time_offset,'.')
hold on
lsline
title(strcat(num2str(R(1,2)),'|',num2str(P(1,2))));

%% Theta sequences in the  "progression to choice" domain
subplot(numrows,numcol,[9*numcol+3 9*numcol+4 10*numcol+3 10*numcol+4])
scatter(compiledthetaComp.tone_offset,compiledthetaComp.tone_ccgs_time_offset,'.')
hold on
lsline
[R,P] = corrcoef(compiledthetaComp.tone_offset,compiledthetaComp.tone_ccgs_time_offset,'rows','complete');
title(strcat(num2str(R(1,2)),'|',num2str(P(1,2))));

%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure4B_ThetaSequences.png'));
saveas(gcf,strcat(expPath,'SupFigure4B_ThetaSequences.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure4B_ThetaSequences.fig'));
end