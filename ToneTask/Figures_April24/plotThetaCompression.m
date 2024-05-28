function plotThetaCompression

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[30 275 1800 550]);

numrows = 6;
numcol = 6;

%% Plot examples of cell sequences that combine place and "tone cells"
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
plotPosPhaseCCG(5,83,103,numrows, numcol, 1, fig2)
% 
% plotPosPhaseCCG(79,186,[],numrows, numcol, 4, fig2)
%     134,186,64,numrows, numcol, 4, fig2)

end