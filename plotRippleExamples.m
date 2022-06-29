function plotRippleExamples

prestim = 1;
analogCh = 2;
numAnalog = 2;
shankNum = 1;
peritime = 0.5;
[colormap] = cbrewer('seq','PuBuGn',100);
colormap(colormap<0) = 0;

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

if exist([sessionInfo.FileName '.ripples.events.mat'],'file') 
   load([sessionInfo.FileName '.ripples.events.mat']);
else
   disp('First calculate .ripple file! Skipping');
   return
end

load([sessionInfo.FileName '.DS2.events.mat']);
ripples = DS2;

if exist([sessionInfo.FileName '.sharpwaves.events.mat'],'file') 
   load([sessionInfo.FileName '.sharpwaves.events.mat']);
else
   disp('No sharp waves associated with this session');
end

% Load pulses
disp('Getting analog-in inputs...');
[pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);

%Skip session if too few ripples
if length(ripples.peaks)<5
    return
end

if exist('pulses')
    if analogCh<=numAnalog
        pulTr = (pulses.stimComb==analogCh);
    else
        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==analogCh);
    end
end

events = pulses.intsPeriods(1,pulTr);
%events = events(((events + 5) <= max(ripples.peaks)) & ((events - 5) > 0));

%Generate logicals for ripples in pre versus post
ripple_pre = [];
ripple_post = [];
ripple_pre(1:length(ripples.peaks)) = 0;
ripple_post(1:length(ripples.peaks)) = 0;  

for pp = 1:length(ripples.peaks)
    tempDiff = ripples.peaks(pp) - events;

    if min(abs(tempDiff)) <=4.8 && min(abs(tempDiff)) > 0.2% If a ripple occurs within 5 seconds of a stimulus
       [~,idxmin] =  min(abs(tempDiff));
       if tempDiff(idxmin) > 0
           ripple_post(pp) = 1;
       elseif tempDiff(idxmin) < 0
           ripple_pre(pp) = 1;
       end
    end
end

ripple_pre = logical(ripple_pre);
ripple_post = logical(ripple_post);

if prestim == 1
    ripplesToPlot = ripples.peaks(ripple_pre);
else
    ripplesToPlot = ripples.peaks(ripple_post);
end

% Now you have the timestamps of the relevant ripples. Plot +-500 ms around
% the ripple event
lfp = bz_GetLFP(sessionInfo.AnatGrps(shankNum).Channels,'noPrompts', true);
lfp = bz_interpolateLFP(lfp);
data = lfp.data;
timestamps = lfp.timestamps;

for pp = 1:length(ripplesToPlot)
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'color','w')
    hold on
    for kk = 1:(length(sessionInfo.AnatGrps(shankNum).Channels)-2)
        plot(timestamps(((ripplesToPlot(pp)-peritime)*1250):((ripplesToPlot(pp)+peritime)*1250)),1*(data(((ripplesToPlot(pp)-peritime)*1250):((ripplesToPlot(pp)+peritime)*1250),kk))-(kk-1)*450,'Color',colormap(kk+22,:),'LineWidth',1.1)
    end
    line([ripplesToPlot(pp) ripplesToPlot(pp)],[2500 (-4*10^4)],'Color','red')
end
end