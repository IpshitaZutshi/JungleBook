
trialNumber = 42; %port6:(45, deliberation; 51, no delib); Port 4: (35, 39, 60 - deliberation; 42, 26 - no delib)

% filename = dir(['*Tracking.Behavior.mat']); 
% load(filename(1).name);
% 
% file = dir(['*TrialBehavior.Behavior.mat']);
% load(filename(1).name);
% 
% file = dir(['*spikeData.cellinfo.mat']);
% load(filename(1).name);
    
aviFile = dir(['front*avi']); 
aviFile = erase(aviFile.name,'.avi');
videoObj1 = VideoReader([aviFile '.avi']);   % get video

aviFile = dir(['test*avi']); 
aviFile = erase(aviFile.name,'.avi');
videoObj2 = VideoReader([aviFile '.avi']);

cellNum1 = 28;
cellNum2 = 4;

writerObj = VideoWriter(strcat('toneTrialVideo_',num2str(trialNumber),'_cellNum_',num2str(cellNum1)),'MPEG-4');
writerObj.FrameRate = 30;
% open the video writer
open(writerObj);


[~,startFrame] = min(abs(tracking.timestamps-behavTrials.timestamps(trialNumber,1)));
[~,endFrame] = min(abs(tracking.timestamps-behavTrials.timestamps(trialNumber,2)));

plotSpikes = 1;

fig = figure;
set(fig,'Position',[1935 53 1890 939])
set(fig,'Color','w')

subplot(3,6,[1 2 7 8 13 14]);
set(gca,'visible','off')

subplot(3,6,[3 9 15]);
set(gca,'visible','off')

ax1 = subplot(3,6,[4 5 6]);
axis off
title(ax1,'Y position','Color','k')
%xlabel(ax1,'Time(s)','Color','w')
set(ax1,'XLim',[tracking.timestamps(startFrame) tracking.timestamps(endFrame)]);
set(ax1,'YLim',[min(tracking.position.y) max(tracking.position.y)]);
h1 = animatedline('Color',[0.3 0.3 0.3],'LineWidth',0.8);

ax2 = subplot(3,6,[10 11 12]);
title('X position','Color','k')
axis off
%set(gca,'visible','off')
%xlabel('Time(s)','Color','w')
set(ax2,'XLim',[tracking.timestamps(startFrame) tracking.timestamps(endFrame)]);
set(ax2,'YLim',[min(tracking.position.x) max(tracking.position.x)]);
h2 = animatedline('Color',[0.3 0.3 0.3],'LineWidth',0.8);

ax3 = subplot(3,6,[16 17 18]);
title('Velocity','Color','k')
axis off
%set(gca,'visible','off')
%xlabel('Time(s)','Color','w')
set(ax3,'XLim',[tracking.timestamps(startFrame) tracking.timestamps(endFrame)]);
set(ax3,'YLim',[min(tracking.position.v) max(tracking.position.v)]);
h3 = animatedline('Color',[0.3 0.3 0.3],'LineWidth',0.8);


for ii = startFrame:endFrame
    
    subplot(3,6,[1 2 7 8 13 14]);
    frame = read(videoObj1,ii);
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
        if sum(ismember(spikeData.posIdx{cellNum1}, ii))>0
            ax1 = subplot(3,6,[4 5 6]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.y(ii),15,'r',"filled")
            
            ax2 = subplot(3,6,[10 11 12]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.x(ii),15,'r',"filled")
            
            ax3 = subplot(3,6,[16 17 18]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.v(ii),15,'r',"filled")
        end
        
        if sum(ismember(spikeData.posIdx{cellNum2}, ii))>0
            ax1 = subplot(3,6,[4 5 6]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.y(ii),15,'b',"filled")
            
            ax2 = subplot(3,6,[10 11 12]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.x(ii),15,'b',"filled")
            
            ax3 = subplot(3,6,[16 17 18]);
            hold on
            scatter(tracking.timestamps(ii),tracking.position.v(ii),15,'b',"filled")
        end
        
    end
       
    subplot(3,6,[1 2 7 8 13 14]);
    title(strcat('Frame number: ',num2str(ii)),'Color','w')
    M = getframe(fig);   
    writeVideo(writerObj, M);
end

close(writerObj)
