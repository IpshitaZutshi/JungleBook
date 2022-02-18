function M = temp(lfp)

[~,idx] = min(abs(lfp.timestamps-6337));
data = lfp.data(idx:idx+2000,:);
x = 1:size(data,1);
[colormap] = cbrewer('seq','YlGnBu',90);%PuBuGn
colormap(colormap<0) = 0;

fig = figure;
set(gcf,'Color','w') 
hold on 
set(gcf, 'Position',[10 10 1400 800])
ax1 = subplot(1,1,1);
set(ax1,'XLim',[0 1250])
set(gca,'visible','off')

for ii = 1:63    
    arr(ii) = (ii-1)*1000;
end
data = data(:,[1:7 9:64]);
data = 1.8*double(data)-arr;

for ii = 1:63
    h(ii) = animatedline('Color',colormap(22+ii,:),'LineWidth',1.1);%,'MaximumNumPoints',1250);
end

for j = 1:63
    addpoints(h(j),x(1:1250),data(1:1250,j));
    hold on
end

dt = 2;%1250/1000;
frame = 1;
for k = 1250:dt:(length(x)-dt)
    for j = 1:63
        addpoints(h(j),x(k:k+dt),data(k:k+dt,j));
        hold on
    end
    set(ax1,'XLim',[k-1100 k+100]);
    drawnow 
    xlabel('Time')
    M(frame) = getframe(fig);    
    frame = frame+1;
end

close;

V = VideoWriter('C:\Users\Ipshita\Desktop\trace','MPEG-4');
V.FrameRate = 40;
V.Quality = 100;
open(V)
writeVideo(V,M)
close(V)


%close all