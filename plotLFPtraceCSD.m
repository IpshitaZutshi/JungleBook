function M = plotLFPtraceCSD(lfp,events)

% cd('Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ27\Saline\IZ27_864um_201022_sess5')
% [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
% jj = 1;
% lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels,'noPrompts', true);
% lfp = bz_interpolateLFP(lfp);
% analogEv = 64;
% numAnalog = 2;
% for ii = 1:numAnalog
%     analogCh(ii) = (analogEv-1)+ii;
% end
% 
% [pulses] = bz_getAnalogPulses('analogCh',analogCh);
% pulTr = (pulses.stimComb==2);
% events = pulses.intsPeriods(:,pulTr);
% events = round(events*1250);

file = dir(['*Behavior.mat']);
load(file(1).name);

[colormap] = cbrewer('seq','YlGnBu',90);%PuBuGn'BuPu'
colormap(colormap<0) = 0;

for ii = 1:62
    arr(ii) = (ii-1)*1000;
end

preTime = 8;
totalTime = 65;
eventNum = 13;
startPoint = events(1,eventNum)-preTime*1250;

% Actual TS of startPoint
TS = lfp.timestamps(startPoint:startPoint+(totalTime*1250));
stimTS = startPoint:startPoint+(totalTime*1250);
data = double(lfp.data(startPoint:startPoint+(totalTime*1250),:))*20-arr;

x = 1:length(TS);
% 
% %Find closest index during behavior
[~,idxStart] = min(abs(behavior.timestamps-TS(1)));
[~,idxEnd] = min(abs(behavior.timestamps-TS(end)));
bX = interp1(behavior.timestamps(idxStart:idxEnd),behavior.position.x(idxStart:idxEnd),TS);
bY = interp1(behavior.timestamps(idxStart:idxEnd),behavior.position.y(idxStart:idxEnd),TS);

stim = stimTS;
for ii = 1:size(events,2)
    if ismember(events(1,ii),stimTS)
        startEvent = find(stimTS == events(1,ii));
        endEvent = find(stimTS == events(2,ii));
        stim(startEvent:endEvent) = -1;
    end
end
stim(stim >=0) = 0;
stim = abs(stim);


fig = figure;
set(gcf,'Color','w') 
%set(gcf,'renderer','painters')
%set(gcf,'Color','None') 
hold on 
set(gcf, 'Position',[10 10 1400 800])

%% Position definition for LFP trace

axPosn = [];
for ii = 40:4:(62*4)
    axPosn = [axPosn ii+2 ii+3 ii+4];
end
%ax1 = subplot(65,1,5:65);
ax1 = subplot(65,4,axPosn);
set(ax1,'XLim',[1 1250]);
%set(ax1,'YLim',[-6.3^10^4 0]);
set(gca,'visible','off')
for ii = 1:62
    h(ii) = animatedline('Color',colormap(22+ii,:),'LineWidth',1.1);%,'MaximumNumPoints',1250);
end

%% Position definition for stim
%ax2 = subplot(65,1,[2 3 4]);
ax2 = subplot(65,4,[6 7 8 10 11 12 14 15 16 18 19 20]);
set(ax2,'XLim',[1 1250]);
set(gca,'visible','off')
xS = animatedline('Color','b','LineWidth',3);%,'MaximumNumPoints',1250);
ylim([0.8 1.3])
xlim([1 1250])

%% Position definition for behavior

ax3 = subplot(65,4,[101:4:(55*4)]);
img = imread('C:\Users\Ipshita\Desktop\maze.png');
img = imresize(img,[55 75]);
image(img);

set(ax3,'XLim',[0 80]);
set(ax3,'YLim',[0 60]);
set(ax3,'YDir','reverse');
set(gca,'visible','off')
xB = animatedline('Color',[0.5 0.5 0.5],'LineWidth',3);%,'MaximumNumPoints',1250);

%ax4 = subplot(65,1,1);
ax4 = subplot(65,4,[1 5 9 13]);
set(ax4,'XLim',[-1 1]);
set(ax4,'YLim',[-1 1]);
set(gca,'visible','off')

%% initial plot for stim and data
addpoints(xS,x(1:1250),stim(1:1250));

for j = 1:62
    addpoints(h(j),x(1:1250),data(1:1250,j));
    hold on
end

addpoints(xB,bY(1:1250),bX(1:1250));

%% Plot everything
dt = 250;%1250/1000;
frame = 1;
flag = 0;
h1 = text(ax4,0,0,'','Color','blue','FontSize',16);

for k = 1250:dt:(length(x)-dt)
    for j = 1:62
        addpoints(h(j),x(k:k+dt),data(k:k+dt,j));
        hold on
    end
    addpoints(xS,x(k:k+dt),stim(k:k+dt));
    addpoints(xB,bY(k:k+dt),bX(k:k+dt));    
    set(ax2,'XLim',[k-1150 k+100]);
    if sum(stim(k:k+dt)) > 1 && flag == 1
        delete(h2)
        h1 = text(ax4,-1,0,'Stim ON','Color','blue','FontSize',25);
        flag = 0;
    elseif sum(stim(k:k+dt))==0 && flag == 0
        delete(h1)
        h2 = text(ax4,-1,0,'Stim OFF','Color',[0.5 0.5 0.5],'FontSize',25);
        flag = 1;
    end
    
    set(ax1,'XLim',[k-1150 k+100]);
    %set(ax1,'YLim',[-6.3^10^4 0]);

    drawnow %limitrate
%     if b > (1/30)
%         drawnow % update screen every 1/30 seconds
%         a = tic; % reset timer after updating
%     end
%      xlim([k k+1250])
%     ylim([-7*10^4 0])
    M(frame) = getframe(fig);    
    frame = frame+1;
end
close all

V = VideoWriter('C:\Users\Ipshita\Desktop\CSDtrace','MPEG-4');
V.FrameRate = 5;
V.Quality = 100;
open(V)
writeVideo(V,M)
close(V)
% % 
% [h, w, ~] = size(M(1).cdata);  
% hf = figure;
% set(hf, 'position', [150 150 w h]);
% axis off
% movie(hf,M,1,5);

end
%drawnow

% 
% 
% 
% startPoint = 47000;
% x = 1:size(lfp.data(startPoint:(startPoint+1250*5),1),1);
% Dx = 1250/2;
% figure
% hold all
% set(gcf, 'Position',[10 10 800 800])
% for jj = x
%     for ii = 1:64
%         plot(double(lfp.data(startPoint+jj,ii))-arr(ii), 'Color',colormap(22+ii,:),'LineWidth',1.2)   
%     end
%     xlim([x(jj) x(jj+Dx)]);drawnow
% end

