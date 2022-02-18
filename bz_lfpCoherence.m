function [coherence] = bz_lfpCoherence(pyrCh,lfp,varargin)
% Computes coherence of the lfp in the theta channel versus all other lfp
% channels
%
%   INPUTS
%
%    pyrCh      Channel to calculate coherence with
%    lfp        a buzcode-formatted lfp structure (use bz_GetLFP)
%               needs fields: lfp.data, lfp.timestamps, lfp.samplingRate, lfp.channels.
%               If empty or no exist, look for lfp in basePath folder
%    intervals  list of time intervals [0 10; 20 30] to read from 
%                   the LFP file (default is [0 inf])
%    fRange     frequency range (in Hz) (default = [0 200])
%    window     size of the silding time window (s, default 2)
%    overlap    overlap between successive windows (default 1)
%    step       step between successive windows (default 1)
%    tapers      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%    pad        FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%    showfig    true/false - show a summary figure of the results
%               (default:false)
%    saveMat    put your basePath here to save/load(default: true)
%    forceDetect (default false)
%
%   OUTPUTS
%   coherence
%       .mean       log10-transformed mean amplitude of the spectrogram
%       .std
%       .median
%       .channels
%       .channels_shank
%       .frange
%       .f
%       .ic95
%
% IZ-BuzsakiLab 2021

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Parsing inputs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
p = inputParser;
addParameter(p,'fRange',[0 200],@isnumeric)
addParameter(p,'window',2,@isscalar)
addParameter(p,'intervals',[0 Inf])
addParameter(p,'overlap',1, @isscalar)
addParameter(p,'step',1, @isscalar)
addParameter(p,'tapers',[3 5], @isnumeric)
addParameter(p,'pad',0, @isscalar)
addParameter(p,'showfig',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceDetect',false,@islogical)

parse(p,varargin{:})
showfig = p.Results.showfig;
saveMat = p.Results.saveMat;
fRange = p.Results.fRange;
window = p.Results.window;
overlap = p.Results.overlap;
step = p.Results.step;
tapers = p.Results.tapers;
pad = p.Results.pad;
forceDetect = p.Results.forceDetect;
intervals = p.Results.intervals;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Resolving inputs   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if exist([sessionInfo.FileName,'.coherence.channelinfo.mat'],'file') && ~forceDetect
    disp(['Coherence already calculated for ', sessionInfo.FileName,'.Loading file.']);
    load([sessionInfo.FileName,'.coherence.channelinfo.mat']);
    return
end
     

%% Calculate coherence
Fs = sessionInfo.rates.lfp; 
lfpPyrCh = double(lfp.data(:,lfp.channels==pyrCh));

disp('Calculating coherograms channelwise')
for ii = 1:length(lfp.channels)
    lfpCurr = double(lfp.data(:,ii));
    [coher,ph,t,f] = bz_MTCoherogram(lfpPyrCh,lfpCurr,'frequency',Fs,'range',fRange,'window',window,...
        'overlap',overlap,'step',step,'tapers',tapers,'pad',pad);
    
    %Interpolate the time stamps to match the LFP timestamps
    assumedLFPtimestamps = [0:length(lfp.data(:,ii))-1]./lfp.samplingRate;
    timestamps = interp1(assumedLFPtimestamps,lfp.timestamps,t,'nearest');
    keeptimes = InIntervals(timestamps,intervals);
    
    % Only keep relevant interval data
    coherogram(:,:,ii) = coher(:,keeptimes);
    phase(:,:,ii) = ph(:,keeptimes);
    timestamps = timestamps(keeptimes);
    
end

coherence.coherogram = coherogram;
coherence.phase = phase;
coherence.timebins = timestamps;
coherence.frequency = f;

%Average across time to show along channels
avgCoherogram = mean(coherogram,2);
avgCoherogram = reshape(avgCoherogram,length(f),length(lfp.channels));
coherence.avgCoherogram = avgCoherogram';

avgPhase = circ_mean(phase,[],2);
avgPhase = reshape(avgPhase,length(f),length(lfp.channels));
coherence.avgPhase = avgPhase';

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Saving the result to basename.coherence.channelinfo.mat

coherence.date = now;
coherence.params.Fs = Fs;
coherence.params.range = fRange;
coherence.params.window = window;
coherence.params.step = step;
coherence.params.overlap = overlap;
coherence.params.tapers = tapers;
coherence.params.pad = pad;

if saveMat
    save([sessionInfo.FileName,'.coherence.channelinfo.mat'],'coherence','-v7.3');
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if showfig
    figure('Position', [20 20 600 600],'Name',sessionInfo.FileName)
    subplot(2,1,1);
    imagesc(f,1:length(lfp.channels),coherence.avgCoherogram)
    xlabel('Frequency(Hz)')
    ylabel('Channel')
    title('Coherogram Amplitude');
    
    subplot(2,1,2);
    imagesc(f,1:length(lfp.channels),coherence.avgPhase)
    xlabel('Frequency(Hz)')
    ylabel('Channel')
    title('Coherogram Phase');
    
    if ~exist('SummaryFigures','dir')
        mkdir('SummaryFigures')
    end
    saveas(gcf,['SummaryFigures\Coherence.png']);
end

end