behav = bz_decodeBehavior_iz(pwd,spikes,tracking,'dim_names',{'y'},'dim_bins',120,'time_smooth',0.5);
ts = behavTrials.timestamps(behavTrials.linTrial==0,2);

bins = 5./median(diff(behav.timestamps));
decodPos = [];
ranPos = [];
for ii = 1:length(ts)
    [~,idx] = min(abs(behav.timestamps-ts(ii)));
    decodPos(ii,:) = behav.errors(idx-bins:idx+bins);   
    ranPos(ii,:) = behav.errors_rand(idx-bins:idx+bins);
end

figure
plot(nanmedian(decodPos,1))
hold on
plot(nanmedian(ranPos,1))