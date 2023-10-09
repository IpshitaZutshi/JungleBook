function plotPGAM_Figure6

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrows = 21;
numcol = 21;

pgamPATH = 'C:\Data\PGAMAnalysis\processedData\';

%% Panel A: Plot example cells, and the PGAM fits
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220830_sess6';
cd(sessloc)
cellNum = 41;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(1,1, firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(1,3,results,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220828_sess5';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
cellNum = 38;%32;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(1,8,firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(1,10,results,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 28;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(1,15,firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(1,17,results,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
cellNum = 275;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(7,1,firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(7,3,results,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220901_sess8';
cd(sessloc)
cellNum = 46;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(7,8,firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(7,10,results,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 134;
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCellAvgMaps(7,15,firingMaps,fig2,numrows,numcol,cellNum)
a = strsplit(sessloc,'\');
load(strcat(pgamPATH,a{end},'\results_struct.mat'));
plotPGAMExampleFits(7,17,results,fig2,numrows,numcol,cellNum)


%% Extract PGAM data
Summary = compileTuningProportion('plotfig',false);

%% Plot A: tuning to each variables

subplot(numrows,numcol,[(1:1:15)+(numcol*12) (1:1:15)+(numcol*13)  (1:1:15)+(numcol*14)  (1:1:15)+(numcol*15) (1:1:15)+(numcol*16) (1:1:15)+(numcol*17)])
hold on

% Scatter plots
colors = [...
    1.0000    0.7333    0.6196; ... 
    0.8510    0.3255    0.0980; ... 
    1.0000    0.8902    0.6314; ... 
    0.9294    0.6941    0.1255; ... 
    0.6667    0.5098    0.7020; ... 
    0.4941    0.1843    0.5569; ... 
    0.6980    0.8196    0.5412; ... 
    0.4667    0.6745    0.1882; ... 
    0.5216    0.6824    0.7882; ... 
    0    0.4471    0.7412; ...
    0   0    0; ... 
    0.5    0.5    0.5];

c = 1;

% Extract tuning proportions per mouse and sessions
names = {'IZ39','IZ40','IZ43','IZ44','IZ47','IZ48'};

for mm = 1:6
    x = 0.25;
    propSigMouse = Summary.propSigAll(Summary.mouseIDProp==mm,:);
    propSigAvg = nanmean(propSigMouse,1);

    color = colors(c, :);
    colorL = colors(c+1, :);
    for v = 1:size(propSigMouse,2)
        s(mm) = scatter(ones(size(propSigMouse(:,v)))*x, propSigMouse(:,v), 14, color, 'filled');
        xvals = [x-0.03, x+0.03];
        yvals = repmat(propSigAvg(v),1,2);
        plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2)
        x = x + 0.25;   
    end
    legend(s, names, 'AutoUpdate', 'off', 'FontSize', 12)
    c = c + 2;
end

% Boxplots
positions = [0.35, 0.6, 0.85, 1.1, 1.35, 1.6, 1.85, 2.1 2.35];
for v = 1:size(Summary.propSigAll, 2)
   b = boxplot(Summary.propSigAll(:,v), 'Positions', positions(v), 'PlotStyle', 'traditional', ...
        'Colors', [0.5216 0.6824 0.7882], 'Symbol', ".", 'Widths', 0.09);
    set(b, {'LineWidth'}, {1.3})
end

h = findobj('LineStyle','--'); 
set(h, 'LineStyle','-');
g = findobj(gca,'Tag','Box');
for j = 1:length(g)
    patch(get(g(j),'XData'), get(g(j),'YData'), [0.7451 0.8392 0.9020],'FaceAlpha',.5);
end

xticks([0.3, 0.55, 0.8, 1.05, 1.3, 1.55, 1.8, 2.05 2.3])
xticklabels(["y", "yRev","yNoTone","relDistStop", "licks", "y relDistStop", ...
        "y licks","licks relDistStop","y licks relDistStop"])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel', a,'FontSize', 11);
xlim([0.2, 2.5])
ylim([0,1])
box off

%% Plot B: tuning to both relDistStop and y 
c = 1;
subplot(numrows,numcol,[(16:1:21)+(numcol*12) (16:1:21)+(numcol*13) (16:1:21)+(numcol*14) (16:1:21)+(numcol*15) (16:1:21)+(numcol*16) (16:1:21)+(numcol*17)])
hold on
for mm = 1:6
    x = 0.25;
    propSigMouse = Summary.propSigBoth(Summary.mouseIDProp==mm,:);
    propSigAvg = nanmean(propSigMouse,1);

    color = colors(c, :);
    colorL = colors(c+1, :);
    for v = 1:size(propSigMouse,2)
        s(mm) = scatter(ones(size(propSigMouse(:,v)))*x, propSigMouse(:,v), 14, color, 'filled');
        xvals = [x-0.03, x+0.03];
        yvals = repmat(propSigAvg(v),1,2);
        plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2)
        x = x + 0.25;   
    end
    c = c + 2;
end

% Boxplots
positions = [0.35, 0.6, 0.85, 1.1];
for v = 1:size(Summary.propSigBoth, 2)
   b = boxplot(Summary.propSigBoth(:,v), 'Positions', positions(v), 'PlotStyle', 'traditional', ...
        'Colors', [0.5216 0.6824 0.7882], 'Symbol', ".", 'Widths', 0.09);
    set(b, {'LineWidth'}, {1.3})
end

h = findobj('LineStyle','--'); 
set(h, 'LineStyle','-');
g = findobj(gca,'Tag','Box');
for j = 1:length(g)
    patch(get(g(j),'XData'), get(g(j),'YData'), [0.7451 0.8392 0.9020],'FaceAlpha',.5);
end

xticks([0.3, 0.55, 0.8, 1.05])
xticklabels(["y Only", "Dist only", "y > relDistStop", "relDistStop > y"])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel', a,'FontSize', 11);
xlim([0.2, 1.2])
ylim([0,1])
box off

%% Plot C: Distribution of mutual information
idxBoth  =  Summary.sigAll(:,1) & Summary.sigAll(:,4);
idxSpaceOnly  = Summary.sigAll(:,1) & ~Summary.sigAll(:,4); 
idxDistOnly = ~Summary.sigAll(:,1) & Summary.sigAll(:,4); 

infoExtract = Summary.mutInfoAll(idxBoth,:);
InfoBothSpace = infoExtract(:,1)>infoExtract(:,4);
infoBothDist  = infoExtract(:,4)>infoExtract(:,1);

subplot(numrows,numcol,[(1:1:4)+(numcol*18) (1:1:4)+(numcol*19) (1:1:4)+(numcol*20)])
dist{1} = log2(Summary.mutInfoAll(idxSpaceOnly,1)./Summary.mutInfoAll(idxSpaceOnly,4));
dist{2} = log2(Summary.mutInfoAll(idxDistOnly,1)./Summary.mutInfoAll(idxDistOnly,4));
nhist(dist,'probability','samebins')

subplot(numrows,numcol,[(5:1:8)+(numcol*18) (5:1:8)+(numcol*19) (5:1:8)+(numcol*20)])
scatter(Summary.mutInfoAll(idxBoth,1),Summary.mutInfoAll(idxBoth,4),10,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.7)
ylim([0 3.6])
xlim([0 3.6])
refline(1)
xlabel('Mut Info, y')
ylabel('Mut Info, Dist to stop')

subplot(numrows,numcol,[(9:1:12)+(numcol*18) (9:1:12)+(numcol*19) (9:1:12)+(numcol*20)])
data = log2((infoExtract(:,1)./infoExtract(:,4)));%./(infoExtract(:,1)+infoExtract(:,4));
histogram(data,-3:0.3:5,'Normalization','probability')
xlabel('Log(2) of ratio of mutual info')
ylabel('Proportion')
box off

%% Look at lick kernels

idxLick  =  Summary.sigAll(:,10);
kerStrength = Summary.kernelStrengthAll(idxLick==1,5:10);

strengthChoice = kerStrength(:,1);
strengthSpont = mean([kerStrength(:,2) kerStrength(:,4)],2);
strengthHome = kerStrength(:,5);

subplot(numrows,numcol,[(13:1:16)+(numcol*18) (13:1:16)+(numcol*19) (13:1:16)+(numcol*20)])
scatter(strengthChoice, strengthHome,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
lsline
xlabel('Ker Strength, Choice licks')
ylabel('Ker Strength, Home licks')
[R,p] = corrcoef(strengthChoice, strengthHome);
title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))

subplot(numrows,numcol,[(17:1:20)+(numcol*18) (17:1:20)+(numcol*19) (17:1:20)+(numcol*20)])
scatter(strengthChoice, strengthSpont,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
xlim([-1 1])
ylim([-1 1])
lsline
xlabel('Ker Strength, Choice licks')
ylabel('Ker Strength, Spont licks')
[R,p] = corrcoef(strengthChoice, strengthSpont);
title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))


% subplot(3,5,13)
% scatter(strengthSpont, strengthHome,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
% xlim([-1 1])
% ylim([-1 1])
% lsline
% xlabel('Ker Strength, Spont licks')
% ylabel('Ker Strength, Home licks')
% [R,p] = corrcoef(strengthSpont, strengthHome);
% title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure6_PGAM.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure6_PGAM.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure6_PGAM.fig'));
%save(strcat(expPath,'\Compiled\Figures\Figure4_taskTuning',num2str(groupStyle),'.mat'),'Fig4Stats'); 

end
