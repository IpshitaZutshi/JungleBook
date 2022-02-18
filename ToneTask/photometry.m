try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
    board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
    read_Intan_RHD2000_file_2;
catch
    disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end

%[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
%bz_LFPfromDat(pwd,'outFs',1250);

fsample = frequency_parameters.amplifier_sample_rate;
disp('Loading analogin.dat...');
filename = 'analogin.dat';  
num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir(filename);
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v=v*0.000050354; %convert to volt

Fs = 30000;
%Fpass = [396 304];
Wn = 20/(Fs/2);% [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)];%
[b,a] = butter(3,Wn);

F_green=v(4,:);
F_lp = filtfilt(b,a,F_green);
F_lp = downsample(F_lp,Fs/1250);

t = 1:1:length(F_lp);
b = polyfit(1:1:length(F_lp),F_lp,5);
a = polyval(b,t);
dfByF = (F_lp-a)./a;
%F_lp = detrend(F_lp);

file = dir('*TrialBehavior.Events.mat');
load(file.name)
figure
gain = [0 1 2];
timelag = 10;
ts = linspace(-timelag,timelag,(20*1250)+1);

col = [103/243 189/243 170/243;8/243 133/243 161/243;56/243 61/243 150/243];
for ii = 1:3
    idxcorr = find(behavTrials.correct==0 & behavTrials.toneTrial==0 & behavTrials.toneGain==gain(ii) & behavTrials.linTrial==0);
    psth = [];
    for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
        curridx = round(behavTrials.timestamps(idxcorr(aa),1)*1250);
        psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
    end
    subplot(2,8,2*(ii-1)+1)
    %subplot(2,3,1)
    imagesc(ts,1:size(psth,1),zscore(psth,[],2))
    hold on
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.toneGain==0))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==0),25,[103/243 189/243 170/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.lickLoc==1))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==1),25,[8/243 133/243 161/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.lickLoc==2))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==2),25,[56/243 60/243 150/243],'filled');
%     
    caxis([-4 4])
    line([0 0],[0 length(idxcorr)])
    title('Start, incorrect trials')
    subplot(2,8,7)
    %subplot(2,3,3)
    xt = linspace(0,1/1250,(timelag*1250));
    plot(ts,nanmean(zscore(psth,[],2)),'Color',col(ii,:),'LineWidth',1.5)
    hold on
    line([0 0],[-1 3])
    ylim([-1 3])

    idxcorr = find(behavTrials.correct==1 & behavTrials.toneTrial==0 & behavTrials.toneGain==gain(ii)& behavTrials.linTrial==0);
    psth = [];
    for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
        curridx = round(behavTrials.timestamps(idxcorr(aa),1)*1250);
        psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
    end
    subplot(2,8,2*(ii-1)+2)
    %subplot(2,3,2)
    imagesc(ts,1:size(psth,1),zscore(psth,[],2))
    hold on
%     scatter(ones(1,sum(behavTrials.correct==1 & behavTrials.lickLoc==0))*(timelag*2*1250),find(behavTrials.correct==1 & behavTrials.lickLoc==0),25,[103/243 189/243 170/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==1 & behavTrials.lickLoc==1))*(timelag*2*1250),find(behavTrials.correct==1 & behavTrials.lickLoc==1),25,[8/243 133/243 161/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==1 & behavTrials.lickLoc==2))*(timelag*2*1250),find(behavTrials.correct==1 & behavTrials.lickLoc==2),25,[56/243 60/243 150/243],'filled');
     line([0 0],[0 length(idxcorr)])
    caxis([-4 4])
    title('Start, correct trials')
    subplot(2,8,8)
    %subplot(2,3,3)
    plot(ts,nanmean(zscore(psth,[],2)),'Color',col(ii,:),'LineWidth',1.5)
    hold on
    line([0 0],[-1 3])

    idxcorr = find(behavTrials.correct==0 & behavTrials.toneTrial==0 & behavTrials.toneGain==gain(ii) & behavTrials.linTrial==0);
    psth = [];
    for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
        curridx = round(behavTrials.timestamps(idxcorr(aa),2)*1250);
        psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
    end
    subplot(2,8,8+2*(ii-1)+1)
    %subplot(2,3,4)
    imagesc(ts,1:size(psth,1),zscore(psth,[],2))
    line([0 0],[0 length(idxcorr)])
    hold on
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.lickLoc==0))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==0),25,[103/243 189/243 170/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.lickLoc==1))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==1),25,[8/243 133/243 161/243],'filled');
%     scatter(ones(1,sum(behavTrials.correct==0 & behavTrials.lickLoc==2))*(timelag*2*1250),find(behavTrials.correct==0 & behavTrials.lickLoc==2),25,[56/243 60/243 150/243],'filled');
    caxis([-4 4])
    title('End, incorrect trials')
    subplot(2,8,15)
    %subplot(2,3,6)
    plot(ts,nanmean(zscore(psth,[],2)),'Color',col(ii,:),'LineWidth',1.5)
    hold on
    line([0 0],[-1 3])

    idxcorr = find(behavTrials.correct==1 & behavTrials.toneTrial==0 & behavTrials.toneGain==gain(ii) & behavTrials.linTrial==0);
    psth = [];
    for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
        curridx = round(behavTrials.timestamps(idxcorr(aa),2)*1250);
        psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
    end
    subplot(2,8,8+2*(ii-1)+2)
    %subplot(2,3,5)
    imagesc(ts,1:size(psth,1),zscore(psth,[],2))
    line([0 0],[0 length(idxcorr)])
    hold on
    caxis([-4 4])
    title('End, correct trials')
    subplot(2,8,16)
    %subplot(2,3,6)
    plot(ts,nanmean(zscore(psth,[],2)),'Color',col(ii,:),'LineWidth',1.5)
    hold on
    line([0 0],[-1 3])
end

figure
idxcorr = find(behavTrials.linTrial==1);
psth = [];
for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
    curridx = round(behavTrials.timestamps(idxcorr(aa),1)*1250);
    psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
end
subplot(2,2,1)
imagesc(ts,1:size(psth,1),zscore(psth,[],2))
hold on
caxis([-4 4])
title('LTT start')

subplot(2,2,3)
plot(ts,nanmean(zscore(psth,[],2)),'LineWidth',1.5)
hold on
line([0 0],[-1 3])
xlabel('Time(s)')


psth = [];
for aa = 1:(length(idxcorr))%size(behavTrials.timestamps,1) 
    curridx = round(behavTrials.timestamps(idxcorr(aa),2)*1250);
    psth(aa,:) = F_lp(curridx-(timelag*1250):curridx+(timelag*1250));   
end
subplot(2,2,2)
imagesc(ts,1:size(psth,1),zscore(psth,[],2))
hold on
caxis([-4 4])
title('LTT end')

subplot(2,2,4)
plot(ts,nanmean(zscore(psth,[],2)),'LineWidth',1.5)
hold on
line([0 0],[-1 3])
xlabel('Time(s)')




% power profile of pyr channel of all session
% lfpT = bz_GetLFP(16,'noPrompts',true);
% 
% params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
% [S,tLFP,fLFP] = mtspecgramc(single(lfpT.data),[1 0.5],params);
% S = log10(S); % in Db
% S_det= bsxfun(@minus,S,polyval(polyfit(fLFP,mean(S,1),2),fLFP)); % detrending
% thetaPower = mean(S_det(:,(fLFP>6 & fLFP<12)),2);
%         

% F_red=v(5,:);
% F_rlp = filtfilt(b,a,F_red);
% F_rlp = downsample(F_rlp,Fs/1250);
% 
% figure
% subplot(2,1,1)
% plot(F_lp)
% subplot(2,1,2)
% plot(F_rlp)
% 
% % 

% % 
% 
% params.Fs = 30000; params.fpass = [100 400]; params.tapers = [3 5]; params.pad = 1;
% [S,t,f] = mtspecgramc(single(F_green),[2 1],params);
% 
% sig_480=mean(S(:,(f>212 & f<221)),2);
% %a = detrend(sig_480);
% sig_405=mean(S(:,(f>315 & f<322)),2);
% %b = detrend(sig_405);
% 
% figure
% subplot(4,1,1)
% %plot(a)
% plot(t,sig_480)
% title('Ach signal')
% 
% subplot(4,1,2)
% %plot(b)  
% plot(t,sig_405)
% title('control signal')
% 
% subplot(4,1,3)
% %plot(b)    
% plot(tLFP,thetaPower)
% title('theta power')
% 
% subplot(4,1,4)
% %plot(b)    
% imagesc(tLFP,fLFP,S_det')
% 
% set(gcf,'Renderer','painters')
% saveas(gcf,'Power.png');

% [s,f,t]=spectrogram(F_gr_smooth,1*SampleRate,0,1*SampleRate,SampleRate,'yaxis');
% 
% 
% F_green_smooth = downsample(F_green,fsample/1250);
% F_lp = filtfilt(b,a,F_green_smooth);
% 
% F_red=v(5,30*10^5:end);
% F_red_smooth = downsample(F_red,fsample/1250);
% F_red_lp= filtfilt(b,a,F_red_smooth);
% 
% pos=v(2,30*10^5:end);
% pos_ds = downsample(pos,fsample/1250);
% 
% 
% [s,f,t]=spectrogram(F_gr_smooth,1*SampleRate,0,1*SampleRate,SampleRate,'yaxis');
% sig_480=mean(abs(s(find(f>212 & f<221),:)));
% 
% F_red=v(5,30*10^5:end);
% F_r_smooth=smoothdata(F_red,'movmean',5000);
% [s,f,t]=spectrogram(F_r_smooth,1*SampleRate,0,1*SampleRate,SampleRate,'yaxis');
% sig_570=mean(abs(s(find(f>315 & f<323),:)));
% 
% Pos=v(2,:);
% Pos_smooth=smoothdata(Pos,'movmean',5000);
% 
% 
% ACh = downsample(Ach,fsample/1250);
% Achfilt = filtfilt(b,a,ACh);
% fakeACh = downsample(fakeAch,fsample/1250);
% Achfakefilt = filtfilt(b,a,fakeACh);
% 
% mAch  = mean(Achfilt);
% dfbfF = (Achfilt-mAch)./mAch;
% mAchFake  = mean(Achfakefilt);
% dfbfFfake = (Achfakefilt-mAchFake)./mAchFake;
% 
% figure
% subplot(2,1,1)
% plot(dfbfF)
% subplot(2,1,2)
% plot(dfbfFfake)