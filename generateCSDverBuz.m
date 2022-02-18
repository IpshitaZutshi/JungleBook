% Current Source Density Analysis
% Theory from Mitzdorf, U, 1985.
%
% Implementation from Pettersen et al, 2006, Ulbert et al, 2001,
% Rappelsberger, et al , 1981
%
% CSD is Im, a scalar quantity of the dimension A/cm3. 
% d.J = Im
% If Im > 0 (outward currents), a current source results, and if Im < 0, a
% current sink results
%
% Using Hammings three point formula for noise reduction,
% -Im = (sigma)[0.23V(z-2h)+0.08V(z-h)-0.62V(z)+0.08V(z+h)+0.23V(z+2h)]/h2
%
% C = -sigma(D)Phi
% where, D is a (N-2) x N matrix 
%
% Implemented by Ipshita Zutshi, UCSD
% 31st July, 2018

% Description:
%  1) Reads Neuralynx data from N-channel silicon probes, and generates the
%  raw LFP signal
%  2) Next, it filters the signal in the theta range (6-14 Hz)
%  3) Next, the filtered signal from a reference channel (Channel 14 in P89 - signal
%  from fissure) is analyzed, and peaks above 2 standard deviations are
%  extracted (theta peak = phase 0)
%  4) The filtered signal 125 ms before and after each peak is extracted
%  for each channel. The signal across all peaks is averaged.
%  5) CSD of the central N-4 channels is performed by using the averaged
%  LFP 
%  Uses simple CSD, with no boundary assumptions.
%  Examines the theta range (6-14 Hz) by default. 

function cscdata = generateCSDverBuz

NumChannels = 16; % Number of channels on the silicon probe. 
RefChannel = 14; %Reference channel to select theta peaks
NumStd = 2; % Number of standard deviations to consider the signal a peak
fBand = [6 14]; % in Hz
timeBin = 125; % in milliseconds
filepath = 'W:\Ipshita\Silicon Probes\2017-07-16_16-09-52\02. d0 no stim';
%filepath = 'W:\Ipshita\Silicon Probes\2017-07-16_16-09-52\05. d0 stim';

stimSess = 0; %0 if no, 1 if yes

%if isempty(cscdata)
    cscdata = readcsclistCSD(1:1:NumChannels,{filepath});
%end

for i = 1:NumChannels
    cscDataNew{i} = formatcscCSD(cscdata(i),'SUNITS', 'uV');
    filteredCSC{i} = getFilteredLFPDataCSD(cscDataNew{i},fBand,'hilbert');
    
    if i == RefChannel
        %Take the absolute of the filtered signal and calculate the standard
        %deviation
        rmssig = abs(filteredCSC{i}.filtered);
        stdsig = std(rmssig);
        %Find the index of peaks that are more than NumStd standard deviations 
        peakIdx = find(filteredCSC{i}.filtered>(NumStd*stdsig) & [0; diff(filteredCSC{RefChannel}.phase>0)]);
    end
end

% Find the time based on the sampling frequency
Fs = filteredCSC{RefChannel}.Fs;
numBins = (Fs*timeBin)/1000;

%Now average the LFP across peaks for each channel
for i = 1:NumChannels
    k = 1;
    for p = 1:length(peakIdx)
        if (peakIdx(p)-numBins>0) && (peakIdx(p)+numBins)<length(filteredCSC{RefChannel}.filtered)
            LFPsig(:,i,k) = filteredCSC{i}.filtered((peakIdx(p)-numBins):(peakIdx(p)+numBins));
            k = k+1;
        end
    end
end
avgLFPsig = nanmean(LFPsig,3)*-1;

TTL = zeros(k-1,(numBins*2)+1);
if stimSess == 1
    [ON_TS] = GetLaserTTLsCSD(filepath);
    ON_TS = ON_TS/(10^6);
    k=1;
    for p = 1:length(peakIdx)
        if (peakIdx(p)-numBins>0) && (peakIdx(p)+numBins)<length(filteredCSC{RefChannel}.filtered)
            ts = filteredCSC{i}.timestamps(peakIdx(p));
            [tempTTL tempidx] = min(abs(ON_TS-ts));
            if (tempTTL*1000) < timeBin
                TTL(k,numBins+round(((ON_TS(tempidx)-ts)*Fs))) = 1;
            end
            k=k+1;
        end
    end
end

%%Generate the CSD
CSD = diff(avgLFPsig,2,2);

figure
if stimSess == 1
    subplot(2,1,1)
end

cmax = max(max(CSD));
contourf(CSD',40,'LineColor','none'); hold on;
colormap jet;
caxis([-cmax cmax]);
plot([251 251],[1 size(CSD,2)],'--k');
set(gca,'YDir','reverse');
ylim([0.5 NumChannels-1.5])

for i = 2:(NumChannels-1)
    plot(((avgLFPsig(:,i)/(2*max(avgLFPsig(:))))+(i-1)),'k')
end
xlabel('Time');
ylabel('Channel');

if stimSess == 1
    subplot(2,1,2)
    TTLsig = nanmean(TTL,1);
    plot(TTLsig,'r')
    xlim([0 501]);
end
end


