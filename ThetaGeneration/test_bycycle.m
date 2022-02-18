sT = (6000000/1250)+7;
eT = (6000000/1250)+12;
samplingRate = 1250;

figure
subplot(4,1,1)
hold on
plot([sT*samplingRate:eT*samplingRate],data(sT*samplingRate:eT*samplingRate),'k')
plot([sT*samplingRate:eT*samplingRate],lfp_lowpass(sT*samplingRate:eT*samplingRate),'b','LineWidth',2)
plot(shape_features.sample_peak(shape_features.sample_peak>=sT*samplingRate & shape_features.sample_peak<=eT*samplingRate),lfp_lowpass(shape_features.sample_peak(shape_features.sample_peak>=sT*samplingRate & shape_features.sample_peak<=eT*samplingRate)),'r.')
plot(shape_features.sample_last_trough(shape_features.sample_last_trough>=sT*samplingRate & shape_features.sample_last_trough<=eT*samplingRate),lfp_lowpass(shape_features.sample_last_trough(shape_features.sample_last_trough>=sT*samplingRate & shape_features.sample_last_trough<=eT*samplingRate)),'g.')
plot(shape_features.sample_zerox_decay(shape_features.sample_zerox_decay>=sT*samplingRate & shape_features.sample_zerox_decay<=eT*samplingRate),lfp_lowpass(ceil(shape_features.sample_zerox_decay(shape_features.sample_zerox_decay>=sT*samplingRate & shape_features.sample_zerox_decay<=eT*samplingRate))),'m.')
plot(shape_features.sample_zerox_rise(shape_features.sample_zerox_rise>=sT*samplingRate & shape_features.sample_zerox_rise<=eT*samplingRate),lfp_lowpass(ceil(shape_features.sample_zerox_rise(shape_features.sample_zerox_rise>=sT*samplingRate & shape_features.sample_zerox_rise<=eT*samplingRate))),'c.')
xlim([sT*samplingRate eT*samplingRate])

subplot(4,1,2)
plot(shape_features.sample_peak(shape_features.sample_peak>=sT*samplingRate & shape_features.sample_peak<=eT*samplingRate),shape_features.band_amp(shape_features.sample_peak>=sT*samplingRate & shape_features.sample_peak<=eT*samplingRate))
xlim([sT*samplingRate eT*samplingRate])

subplot(4,1,3)
plot([sT*samplingRate:eT*samplingRate],shape_features.is_burst(sT*samplingRate:eT*samplingRate),'k')
xlim([sT*samplingRate eT*samplingRate])
ylim([-1 2])

subplot(4,1,4)
plot(shape_features.is_burst)
%plot([1889695-(1*samplingRate):1889695+(1*samplingRate)],data(1889695-(1*samplingRate):1889695+(1*samplingRate)),'k')