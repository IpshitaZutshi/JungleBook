function plotToneSpeed_SupFigure6

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[204 27 1500 800]);

BuPu=cbrewer('seq', 'BuPu', 11);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

numrows = 6;
numcol = 9;

col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

%% Panel A: Example lick tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39\Final\IZ39_220624_sess10';
cd(sessloc)
cellNum = 19;

file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);

dataMat = [];
dataMatTone = [];

b = linspace(0,125,50);
a = linspace(2000,22000,50);

for kk = 1:6
    subplot(numrows, numcol, [1 2]);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    

    dataMat = [dataMat;firingMaps.forward.rateMaps{cellNum}{1+kk}];
    atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{1+kk},'linear');
    dataMatTone = [dataMatTone;atm];
    
end

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [1+numcol 2+numcol]);
imagesc(1,b,nanmean(dataMat,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat,1))))

% PSTH 
for ii = 1:6
    idx = behavTrials.lickLoc==(ii-1) & behavTrials.linTrial ==0; 
    st = behavTrials.timestamps(idx,2);
    [stccg, tPSTH] = CCG({spikes.times{cellNum} st},[],'binSize',0.1,'duration',2,'norm','rate');
    subplot(numrows, numcol, [1+2*numcol 2+2*numcol]);
    hold on
    plot(tPSTH,stccg(:,2,1)','Color',col(ii,:),'LineWidth',1);   
end

%% Panel B: Speed relationship of that cell

dtime = 0.2;%mean(diff(tracking.timestamps));
win = [tracking.timestamps(1) tracking.timestamps(end)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'units','rate','win',win);
spkRate  = spkData.data(:,cellNum);
% Interpolate
timeS = tracking.timestamps;
vel = interp1(timeS,tracking.position.v,spkData.timestamps)';
[vel_sorted idxv] = sort(vel);
rate_sorted  = spkRate(idxv);
vel(vel>60) = nan;
vel_sorted(vel_sorted>60) = nan;
[b,e] = discretize(vel_sorted,0:5:40);
meanrate = splitapply(@mean,rate_sorted,b');

Greys=cbrewer('seq', 'Greys', 11);
ax1 = subplot(numrows,numcol,[3 4 3+numcol 4+numcol]);
[values, centers] = hist3([vel' spkRate],[20 20],'Normalization','probability');
values_all = sum(sum(values,1));
imagesc(centers{:}, values.'./values_all)
set(gca,'YDir','Normal')
colormap(ax1,Greys)
caxis([0 0.0008])
box off
colorbar
% scatter(vel,spkRate,5,[0.5 0.5 0.5],"filled")
% ylim([-2 70])
hold on
lsline
plot(e(2:end),meanrate,'LineWidth',2,'Color','r')
[r,p] = corrcoef(vel,spkRate,'Rows','pairwise');
title(strcat('r=', num2str(r(1,2)),' p=', num2str(p(1,2))))


%% Panel A2: Example tone tuned cell 2

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220901_sess8';
cd(sessloc)
cellNum = 46;

file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);

dataMat = [];
dataMatTone = [];

b = linspace(0,125,50);
a = linspace(2000,22000,50);

for kk = 1:6
    subplot(numrows, numcol, [6 7]);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    

    dataMat = [dataMat;firingMaps.forward.rateMaps{cellNum}{1+kk}];
    atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{1+kk},'linear');
    dataMatTone = [dataMatTone;atm];
    
end

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [6+numcol 7+numcol]);
imagesc(1,b,nanmean(dataMat,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat,1))))

% PSTH 
for ii = 1:6
    idx = behavTrials.lickLoc==(ii-1) & behavTrials.linTrial ==0; 
    st = behavTrials.timestamps(idx,2);
    [stccg, tPSTH] = CCG({spikes.times{cellNum} st},[],'binSize',0.1,'duration',2,'norm','rate');
    subplot(numrows, numcol, [6+2*numcol 7+2*numcol]);
    hold on
    plot(tPSTH,stccg(:,2,1)','Color',col(ii,:),'LineWidth',1);   
end

%% Panel B: Speed relationship of that cell

dtime = 0.1;%mean(diff(tracking.timestamps));
win = [tracking.timestamps(1) tracking.timestamps(end)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'units','rate','win',win);
spkRate  = spkData.data(:,cellNum);
% Interpolate
timeS = tracking.timestamps;
vel = interp1(timeS,tracking.position.v,spkData.timestamps)';
[vel_sorted idxv] = sort(vel);
vel(vel>60) = nan;
vel_sorted(vel_sorted>60) = nan;
rate_sorted  = spkRate(idxv);
[b,e] = discretize(vel_sorted,0:5:40);
meanrate = splitapply(@mean,rate_sorted,b');

Greys=cbrewer('seq', 'Greys', 11);
ax1 = subplot(numrows,numcol,[8 9 8+numcol 9+numcol]);
[values, centers] = hist3([vel' spkRate],[20 20],'Normalization','probability');
values_all = sum(sum(values,1));
imagesc(centers{:}, values.'./values_all)
set(gca,'YDir','Normal')
colormap(ax1,Greys)
caxis([0 0.004])
colorbar
box off
hold on
lsline
plot(e(2:length(meanrate)+1),meanrate,'LineWidth',2,'Color','r')
[r,p] = corrcoef(vel,spkRate,'Rows','pairwise');
title(strcat('r=', num2str(r(1,2)),' p=', num2str(p(1,2))))


%% Calculate speed relationship
Summary = compileToneSpeed;

subplot(numrows,numcol,[1+3*numcol 2+3*numcol 1+4*numcol 2+4*numcol])
histogram(Summary.r,'Normalization','probability')
title('Non-spatial cells')
box off 

subplot(numrows,numcol,[4+3*numcol 5+3*numcol 4+4*numcol 5+4*numcol])
histogram(Summary.r_int,'Normalization','probability')
title ('Spatial cells')
box off

subplot(numrows,numcol,[6+3*numcol 7+3*numcol 6+4*numcol 7+4*numcol])
data.tone = Summary.r;
data.space = Summary.r_int;
nhist(data,'proportion','samebins')

subplot(numrows,numcol,[9+3*numcol 9+4*numcol])
SupFigSpeedStats = groupStats([{data.tone} {data.space}],[],'inaxis',true);

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_speedTuning.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_speedTuning.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_speedTuning.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\SupFigure6_speedTuning.mat'),'SupFigSpeedStats'); 

end