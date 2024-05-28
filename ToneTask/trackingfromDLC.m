Xpos = A(:,[2 11 14])*0.2732;
Prob = A(:,[4 13 16]);
Ypos = A(:,[3 12 15])*0.2732;

Ts = tracking.timestamps;
gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];

%% Build vector of positions prior to lick

timebefore = 1*tracking.samplingRate;
timeafter = 1*tracking.samplingRate;
trial = 1;
XLick = [];
YLick = [];

XLick_cent = [];
YLick_cent = [];
lickLoc = [];
for ii = 1:length(behavTrials.linTrial)
   if behavTrials.linTrial(ii)==0 && behavTrials.lickLoc(ii)<5 && behavTrials.correct(ii)==1
        currTS = behavTrials.timestamps(ii,2);
        [~,idx] = min(abs(Ts-currTS));

        XLick(trial,:) = Xpos(idx-timebefore:idx+timeafter,1);
        YLick(trial,:) = Ypos(idx-timebefore:idx+timeafter,1)- Ypos(idx,1);
        lickLoc(trial) = behavTrials.lickLoc(ii);
        trial = trial+1;

        XLick_cent(trial,:) = tracking.position.x(idx-timebefore:idx+timeafter,1);
        YLick_cent(trial,:) = tracking.position.y(idx-timebefore:idx+timeafter,1)-tracking.position.y(idx);        
   end
end

figure
set(gcf,'Color','w')

for ii = 1:size(XLick,1)

    subplot(1,2,1)
    plot(XLick(ii,:),YLick(ii,:),'Color',[0.7 0.7 0.7])
    hold on

    subplot(1,2,2)
    plot(XLick_cent(ii,:),YLick_cent(ii,:),'Color',[0.7 0.7 0.7])
    hold on
end

subplot(1,2,1)
z = zeros(1,size(YLick,2));
col = linspace(-1, 1, length(z));
surface([nanmedian(XLick,1);nanmedian(XLick,1)],[nanmedian(YLick,1);nanmedian(YLick,1)],[z;z],[col;col],...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
avgX = nanmedian(XLick,1);
avgY = nanmedian(YLick,1);
scatter(avgX(floor(length(avgX)/2)),avgY(floor(length(avgY)/2)),'o','filled','w');
colorbar
colormap hsv

subplot(1,2,2)
z = zeros(1,size(YLick_cent,2));
col = linspace(-1, 1, length(z));
surface([nanmedian(XLick_cent,1);nanmedian(XLick_cent,1)],[nanmedian(YLick_cent,1);nanmedian(YLick_cent,1)],[z;z],[col;col],...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
avgX = nanmedian(XLick_cent,1);
avgY = nanmedian(YLick_cent,1);
scatter(avgX(floor(length(avgX)/2)),avgY(floor(length(avgY)/2)),'o','filled','w');
colorbar
colormap hsv

