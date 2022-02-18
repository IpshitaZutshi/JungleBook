% Eventually generates perfectly aligned photometry signal, optitrack data,
% and lfp, all sampled at the user defined sampling rate (1250 if using
% buzcode lfp extraction)

clear;clc;

%% path to folder with recordings
basepath = pwd;
photometryFile = [basepath,basepath((end-21):(end-13)),'photometry.mat'];
trackingFile = [basepath,basepath((end-21):(end-14)),'.csv'];
lfpChannel = 11;
Fs = 1250;

%% get photometry signal

% photometry is sampled at 20 kHz, TTLs at 1 kHz

load(photometryFile);

startTS = OpTrTTL.times(find(OpTrTTL.values(:,1)>3,1)); %Find timestamp of first TTL from optitrack to spike2
stopTS = OpTrTTL.times(find(OpTrTTL.values(:,1)>3,1,'last')); %Find timestamp of last TTL from optitrack to spike2

% Align the TTL timestamp with the photometry
[~,startIdx] = min(abs(AChSig1.times-startTS));
[~,stopIdx] = min(abs(AChSig1.times-stopTS));

% Extract the relevant data
fluorTS = AChSig1.times(startIdx:stopIdx);
fluorRaw = AChSig1.values(startIdx:stopIdx); 

% Now make the timeseries start from 0,
fluorTS = fluorTS - fluorTS(1);

% Confirm the sampling frequency and downsample to the specified sampling
% rate
Fs_photo = 1/(fluorTS(2)-fluorTS(1));
fluorTS = downsample(fluorTS,int16(Fs_photo/Fs));
fluorDS = downsample(fluorRaw,int16(Fs_photo/Fs));

%% filter the signal
Wn = 15/(Fs/2); %[Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]
[b,a] = butter(3,Wn);
fluorFilt = filtfilt(b,a,fluorDS);

% detrend signal to correct for bleaching - use instead of df/f?
% dfByf = detrend(fluorFilt,'linear');

%calculate df/f over a sliding window of 130 seconds

for i = 1:length(fluorFilt)
    if i < (length(fluorFilt)-30*Fs)
        f0 = mean(fluorFilt(i:i+30*Fs));
        dfByf(i,1) = ((fluorFilt(i)-f0)/f0)*100;
    else
        f0 = mean(fluorFilt(i:end));
        dfByf(i,1) = ((fluorFilt(i)-f0)/f0)*100;
    end
end

photometrydata  = timeseries(dfByf,fluorTS,'name','dfByf');

%% get velocity data - sampled at 100 hz, later upsampled to match the lfp and photometry
fps = 100; mvt = 1;

optitrack=csvread(trackingFile,7,1);
idx = find(optitrack(:,2)==0);
optitrack(idx,2:8) = nan(size(idx,1),7);
optitrack = fillmissing(optitrack,'linear');
tS = optitrack(:,1);
x = optitrack(:,6);
y = optitrack(:,7);
vel = [0;sqrt((x(2:end) - x(1:end-1)).^2 + (y(2:end) - y(1:end-1)).^2)*fps];
velSmooth = smoothdata(vel,'movmean',fps*mvt); 
speeddata = timeseries(velSmooth,tS,'name','Speed');

%% read intan data

info = read_Intan_RHD2000_file(pwd);

%Read TTLs from optitrack
num_channelsdigin = info.num_board_dig_in_channels;
digital_ch = read_Intan_digitalin_file('digitalin.dat',num_channelsdigin);

%Read eeg signals
lfp = bz_GetLFP(lfpChannel,'noPrompts',true);

%Extract TTLs
digital = digital_ch(:,1);
digital = downsample(digital,int16(frequency_parameters.amplifier_sample_rate/Fs));

startTTL = find(digital,1);
stopTTL = find(digital,1,'last');

LFP = lfp.data(startTTL:stopTTL);
LFP = double(LFP); %.lfp files are in int16 format and the spectrogram needs a double
timestamps = lfp.timestamps(startTTL:stopTTL);
timestamps = timestamps - timestamps(1);
lfpdata  = timeseries(LFP,timestamps,'name','LFP');

%% Perfectly synchronize all the timestamps & save

[lfp_aligned, photo_aligned] = synchronize(lfpdata,photometrydata,'uniform','interval',1/Fs);
vel_aligned = resample(speeddata,photo_aligned.Time);

save([basepath,'\alignedData'],'vel_aligned','lfp_aligned','photo_aligned'); 

%% get lfp power

movingwin = [5 1];
params.fpass = [1 100];
params.Fs = 1250;
params.tapers = [3 5];

[S,t,f] = mtspecgramc(lfp_aligned.Data,movingwin,params);    
meanThetaPwr = smoothdata(nanmean(S(:,f>=6&f<=12),2));
meanGammaPwr = smoothdata(nanmean(S(:,f>=40&f<=90),2));
Slog = log(S);

%% Plot the data
figure

subplot(3,1,1)
plot(photo_aligned);
ylabel('Fluorescence');
xlabel('Time (s)');

subplot(3,1,2)
plot(vel_aligned);
ylabel('Velocity');
xlabel('Time (s)');

% subplot(5,1,3)
% plot(t,meanThetaPwr);
% ylabel('6-12 Hz power');
% xlabel('Time (s)'); 
% 
% subplot(5,1,4)
% plot(t,meanGammaPwr);
% ylabel('40-90 Hz power');
% xlabel('Time (s)'); 

subplot(3,1,3)
imagesc(t,f,Slog');
ylabel('Frequency');
xlabel('Time (s)');
set(gca,'YDir','normal')
colormap(jet)
% 
% saveas(gcf,[basepath,'\Initial'],'epsc');
% saveas(gcf,[basepath,'\Initial'],'fig');
% saveas(gcf,[basepath,'\Initial'],'jpg');



