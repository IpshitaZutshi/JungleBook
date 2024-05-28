function addSpikePhaseData(varargin)

p = inputParser;
addParameter(p,'lfpChan',63)

parse(p,varargin{:})
lfpChan = p.Results.lfpChan;

file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);

file = dir(['*.spikes.cellinfo.mat']);
load(file.name);

file = dir(['*.Tracking.Behavior.mat']);
load(file.name);
sessionInfo = bz_getSessionInfo;

lfp = bz_GetLFP(lfpChan, 'noprompts',true);
passband = [6 12];
[b, a] = butter(3,[passband(1)/(lfp.samplingRate/2) passband(2)/(lfp.samplingRate/2)],'bandpass'); % order 3
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
%power = fastrms(filt,ceil(lfp.samplingRate./passband(1)));  % approximate power is frequency band
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);
[idxLFP] = InIntervals(lfp.timestamps,[tracking.timestamps(1) tracking.timestamps(end)]); 
lfpphase_beh = lfpphase(idxLFP);
lfpts_beh = lfp.timestamps(idxLFP);

batch_size = 1000;  % Adjust the batch size as needed
for unit = 1:length(spikes.UID)
    [idx] = InIntervals(spikes.times{unit}, [tracking.timestamps(1), tracking.timestamps(end)]); 
    tsBehav = spikes.times{unit}(idx);
    if isempty(tsBehav)
        spikeData.phaseIdx{unit} = [];
    else
        spikeData.phaseIdx{unit} = zeros(size(tsBehav));
        for i = 1:batch_size:length(tsBehav)
            idx_end = min(i+batch_size-1, length(tsBehav));
            batch_tsBehav = tsBehav(i:idx_end);
            % Find the index of the closest timestamp for the current batch
            [~, closestIndex] = min(abs(bsxfun(@minus, lfpts_beh, batch_tsBehav')));
            spikeData.phaseIdx{unit}(i:idx_end) = closestIndex;
        end
    end
    spikeData.phase{unit} = lfpphase_beh(spikeData.phaseIdx{unit});
end

save([sessionInfo.FileName '.spikeData.cellinfo.mat'],'spikeData'); 
end
