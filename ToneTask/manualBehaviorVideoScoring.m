aviFile = dir(['front*avi']); 
aviFile = erase(aviFile.name,'.avi');
videoObj1 = VideoReader([aviFile '.avi']);   % get video


aviFile = dir(['test*avi']); 
aviFile = erase(aviFile.name,'.avi');
videoObj2 = VideoReader([aviFile '.avi']);

filename = dir(['*Tracking.Behavior.mat']); 

writerObj = VideoWriter('videoManualScore','MPEG-4');
writerObj.FrameRate = 30;
% open the video writer
open(writerObj);

load(filename.name)

numFrames = length(tracking.timestamps);
plotSpikes = 1;
cellNum = 4;

fig = figure;
set(fig,'Position',[1935 53 1890 939])
set(fig,'Color','none')


subplot(3,6,[1 2 7 8 13 14]);
set(gca,'visible','off')

subplot(3,6,[3 9 15]);
set(gca,'visible','off')


ax1 = subplot(3,6,[4 5 6]);
axis off
title(ax1,'Y position','Color','w')
%xlabel(ax1,'Time(s)','Color','w')
set(ax1,'XLim',[tracking.timestamps(1) tracking.timestamps(end)]);
set(ax1,'YLim',[min(tracking.position.y) max(tracking.position.y)]);
h1 = animatedline('Color',[0.9 0.9 0.9],'LineWidth',0.8);

ax2 = subplot(3,6,[10 11 12]);
title('X position','Color','w')
axis off
%set(gca,'visible','off')
%xlabel('Time(s)','Color','w')
set(ax2,'XLim',[tracking.timestamps(1) tracking.timestamps(end)]);
set(ax2,'YLim',[min(tracking.position.x) max(tracking.position.x)]);
h2 = animatedline('Color',[0.9 0.9 0.9],'LineWidth',0.8);

ax3 = subplot(3,6,[16 17 18]);
title('Velocity','Color','w')
axis off
%set(gca,'visible','off')
%xlabel('Time(s)','Color','w')
set(ax3,'XLim',[tracking.timestamps(1) tracking.timestamps(end)]);
set(ax3,'YLim',[min(tracking.position.v) max(tracking.position.v)]);
h3 = animatedline('Color',[0.9 0.9 0.9],'LineWidth',0.8);


for ii = 1385:numFrames
    subplot(3,6,[1 2 7 8 13 14]);
    frame = read(videoObj1,ii-7);
    imagesc(frame(50:end,200:450,:))   
    axis off
    
    subplot(3,6,[3 9 15]);
    frame = read(videoObj2,ii);
    imagesc(frame)   
    axis off
    
    addpoints(h1,tracking.timestamps(ii),tracking.position.y(ii));
    addpoints(h2,tracking.timestamps(ii),tracking.position.x(ii));
    addpoints(h3,tracking.timestamps(ii),tracking.position.v(ii));
    
    drawnow limitrate
    
    if plotSpikes
        if sum(ismember(spikeData.posIdx{cellNum}, ii))>0
            ax1 = subplot(3,6,[4 5 6]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.y(ii),5,'r',"filled")
            
            ax2 = subplot(3,6,[10 11 12]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.x(ii),5,'r',"filled")
            
            ax3 = subplot(3,6,[16 17 18]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.v(ii),5,'r',"filled")
        end
    end
       
    subplot(3,6,[1 2 7 8 13 14]);
    title(strcat('Frame number: ',num2str(ii)),'Color','w')
    M = getframe(fig);   
    writeVideo(writerObj, M);
end

close(writerObj)
