function FigS16_projectChangePointUMAP

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[463 2 1204 745]);

numrows = 3;
numcol = 5;

sess = {'IZ43_220828_sess4','IZ44_220830_sess7','IZ47_230707_sess24','IZ48_230628_sess17'};

col1 = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.png'));
saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.fig'));


end

