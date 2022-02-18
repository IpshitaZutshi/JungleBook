function filteredLFPData = getFilteredLFPDataCSD(LFPData,fBand,phaseMethod)
if nargin<3 || isempty(phaseMethod)
    phaseMethod = 'hilbert';
end

signal_filtered = bandpassFilterCSD(LFPData.sample, median(LFPData.Fs), fBand);
[signal_phase signal_amp] = instPhaseCSD(signal_filtered);
if ~strcmpi(phaseMethod,'hilbert')
    signal_phase = phaseFromPeaksCSD(signal_filtered,LFPData.timestamps);
end
filteredLFPData.fBand = fBand;
filteredLFPData.timestamps = LFPData.timestamps;
filteredLFPData.filtered = signal_filtered;
filteredLFPData.phase = signal_phase;
filteredLFPData.amp = signal_amp;
filteredLFPData.Fs = LFPData.Fs;