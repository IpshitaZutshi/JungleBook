%load('alignedData.mat');
%% test
ts_photo = downsample(photo_aligned.Time, 1250);
sig_photo = downsample(photo_aligned.Data, 1250);

[xmaxds,locds] = findpeaks(sig_photo,'MinPeakProminence',3,'MinPeakDistance',5);
[xmax,loc] = findpeaks(photo_aligned.Data,'MinPeakProminence',4,'MinPeakDistance',2*1250);
locrand = randsample(length(photo_aligned.Data),length(loc));

figure(1)
subplot(3,1,1)
plot(photo_aligned.Time,photo_aligned.Data)
hold on
plot(photo_aligned.Time(loc), xmax, 'ro')
title('NE dyByf, red circles are detected NE peaks')

figure(2)
subplot(3,1,1)
plot(photo_aligned.Time,photo_aligned.Data)
hold on
plot(photo_aligned.Time(locrand), photo_aligned.Data(locrand), 'ro')
% subplot(4,1,2)
% plot(ts_photo, sig_photo)
% hold on
% plot(ts_photo(locds), xmaxds, 'r.')

movingwin = [5 1];
params.fpass = [1 100];
params.Fs = 1250;
params.tapers = [3 5];

[S,t,f] = mtspecgramc(lfp_aligned.Data,movingwin,params);    
meanPwr = smoothdata(nanmean(S(:,f>=6&f<=12),2));

Wn = [5/1250 12/1250];
[b,a] = butter(3,Wn);
lfpFilt = filtfilt(b,a,lfp_aligned.Data);

Wn = [40/1250 80/1250];
[b,a] = butter(3,Wn);
lfpFilt2 = filtfilt(b,a,lfp_aligned.Data);

psth_lfp = zeros(length(loc),(1250*2)+1);
psth_power = zeros(length(loc),(1250*2)+1);%zeros(length(locds),10+1);

psth_lfprand = zeros(length(loc),(1250*2)+1);
psth_powerrand = zeros(length(loc),(1250*2)+1);%zeros(length(locds),10+1);

%% psth

for i = 1:length(loc)
    if loc(i) < (length(lfp_aligned.Data)-(1*1250)) && (loc(i) > (1*1250))
        psth_lfp(i,:) = lfpFilt((loc(i)-1*1250):(loc(i)+1*1250));
        psth_power(i,:) = lfpFilt2((loc(i)-1*1250):(loc(i)+1*1250));
    else
        psth_lfp(i,:) = nan;
        psth_power(i,:) = nan;
    end
end

for i = 1:length(locrand)
    if locrand(i) < (length(lfp_aligned.Data)-(1*1250)) && (locrand(i) > (1*1250))
        psth_lfprand(i,:) = lfpFilt((locrand(i)-1*1250):(locrand(i)+1*1250));
        psth_powerrand(i,:) = lfpFilt2((locrand(i)-1*1250):(locrand(i)+1*1250));
    else
        psth_lfprand(i,:) = nan;
        psth_powerrand(i,:) = nan;
    end
end

% for i = 1:length(locds)
%     if locds(i) < length(meanPwr)-5  && (locds(i) > 5)
%         psth_power(i,:) = meanPwr((locds(i)-5):(locds(i)+5));
%     else
%         psth_power(i,:) = nan;
%     end
% end

figure(1)
subplot(3,1,2)
shadedErrorBar([],psth_lfp,{@mean,@std})
hold on
line([(1*1250)+1 (1*1250)+1],[-500 1000],'Color','red')
xlim([0 (2*1250)+1])
xlabel('Time(1 second psth)')

subplot(3,1,3)
shadedErrorBar([],psth_power,{@mean,@std})
line([(1*1250)+1 (1*1250)+1],[-300 300],'Color','red')
xlim([0 (2*1250)+1])
xlabel('Time(1 second psth)')

figure(2)
subplot(3,1,2)
shadedErrorBar([],psth_lfprand,{@mean,@std})
hold on
line([(1*1250)+1 (1*1250)+1],[-500 1000],'Color','red')
xlim([0 (2*1250)+1])

subplot(3,1,3)
shadedErrorBar([],psth_powerrand,{@mean,@std})
line([(1*1250)+1 (1*1250)+1],[-300 300],'Color','red')
xlim([0 (2*1250)+1])
    
figure
subplot(2,1,1)
imagesc(psth_lfp)
subplot(2,1,2)
imagesc(psth_power)
