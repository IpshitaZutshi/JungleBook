function [ MUA ] = MUAfromDat( basePath,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   basePath
%
%   (options)
%       'channels'      (default: all)   (0-indexed, like neuroscope) -> if
%       put all will take only first channel; otherwise specify a channel.
%       update to be able to work with multiple channels. 

%       'saveMat'       (default: true) save a .mat file
%                                   basePath/baseName.MUAGamma.lfp.mat
%       'forceRedetect' (default: false)
%       'peakthresh'    (default: 3.5std)
%       'filterbounds'  (default: [500 5000])
%       'movingNorm'    (default: false)
%       'MUAsmoothwin'  (default: 0.005)
%       'usepeaks'      (default: false)
%       'SHOWFIG'       (default: false)
%       'compareEMG'    (default: false)
%       'peakthresh'    (default: 3.5std)
%
%OUTPUTS

% RS finish description later 

%%
p = inputParser;
addParameter(p,'channels','all');
addParameter(p,'saveMat',true); 
addParameter(p,'filterbounds',[500 5000]);
addParameter(p,'movingNorm',false);
addParameter(p,'MUAsmoothwin',0.005);
addParameter(p,'usepeaks',true);
addParameter(p,'SHOWFIG',true);
addParameter(p,'compareEMG',false);
addParameter(p,'peakthresh',4);
addParameter(p,'forceRedetect',false);

parse(p,varargin{:})

channels = p.Results.channels;
saveMat = p.Results.saveMat;
MUAfilter = p.Results.filterbounds;
movingNorm = p.Results.movingNorm;
MUAsmoothwin = p.Results.MUAsmoothwin;
usepeaks = p.Results.usepeaks;
SHOWFIG = p.Results.SHOWFIG;
compareEMG = p.Results.compareEMG;
peakthresh = p.Results.peakthresh;
forceRedetect = p.Results.forceRedetect;

%% DEV
%basePath = '/mnt/proraidDL/Database/WMProbeData/180211_WT_M1M3_LFP_Layers_Pupil_EMG_Pole/180211_WT_M1M3_LFP_Layers_Pupil_EMG_180211_130605';
%channels = 10;
%note: after 180211 good grounding

%%
sessionInfo = bz_getSessionInfo(basePath, 'noPrompts', true);
datSampleRate = sessionInfo.rates.wideband;
lfpSampleRate = sessionInfo.lfpSampleRate;
%downfactor = datSampleRate./lfpSampleRate;

%%
baseName = bz_BasenameFromBasepath(basePath);

%%
datfilename = fullfile(basePath,[baseName,'.dat']);
savefile = fullfile(basePath,[baseName,'.MUA.lfp.mat']);
if exist(savefile,'file') && ~forceRedetect
   load(savefile)
   return
end

%%
if isequal(channels,'all')
    channels = sessionInfo.channels(1:(length(sessionInfo.channels)-5));
end

%%
lfp = bz_GetLFP(channels(1),'basepath',basePath,'noPrompts',true);

%%
downfactor = 1;
MUA.data = NaN(size(lfp.data,1),length(channels));

if length(channels)>5
    display('Getting MUAs... hold on to your horses')
end
for ii = 1:length(channels)
    disp(['channel' num2str(ii)])
    if length(channels)>5 && mod(ii,5)==1
        display(['Channel ',num2str(ii),' of ',num2str(length(channels))])
    end
    datlfp.data = bz_LoadBinary(datfilename,...
        'frequency',datSampleRate,'nchannels',sessionInfo.nChannels,...
        'channels',channels(ii)+1,'downsample',downfactor);
    
    datlfp.samplingRate = datSampleRate./downfactor;
    datlfp.timestamps = [0:(length(datlfp.data)-1)]'/datlfp.samplingRate;  %To be overwritten later...
%     datlfp.channels = channels;
    
    %%
    MUALFP = bz_Filter(datlfp,'passband',MUAfilter,'filter','fir1','order',3);
%     MUALFP.channels = channels;
    
    if movingNorm
        %Normalize and re-hilbert
        MUALFP.data = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate,'moving',0.2);
        MUALFP.hilb = hilbert(MUALFP.data);
        MUALFP.amp = abs(MUALFP.hilb);
    end
    
    %%
    if usepeaks
        normMUA = NormToInt(MUALFP.data,'modZ',[0 Inf],MUALFP.samplingRate);
        [MUALFP.peakmags,MUALFP.peaktimes] = findpeaks(-normMUA,...
            'MinPeakHeight',peakthresh,'MinPeakDistance',MUALFP.samplingRate./1000);
        MUALFP.peaktimes = MUALFP.timestamps(MUALFP.peaktimes);

        %MUAspks = bz_SpktToSpkmat({MUALFP.peaktimes},'binsize',0.02);

    end
    %% Smooth the MUA
    MUALFP.smoothamp = smooth(MUALFP.amp,round(MUAsmoothwin.*MUALFP.samplingRate),'moving' );

    
    %% Threshold crossings
    % figure
    % plot(MUALFP.smoothamp,MUALFP.data,'.')
    
    
    %% Interpolate to lfp timestamps
    
    MUA.data(:,ii) = interp1(MUALFP.timestamps,MUALFP.smoothamp,lfp.timestamps);
    MUA.peaks.times{ii} = MUALFP.peaktimes;
    MUA.peaks.magnitudes{ii} = MUALFP.peakmags;
end

%% Saving onto
if usepeaks
    for cc = 1:length(channels)
        groups{cc}=channels(cc).*ones(size(MUA.peaks.times{cc}));
    end
    alltimes = cat(1,MUA.peaks.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    MUA.peaks.spindices = [alltimes groups];
end


MUA.timestamps = lfp.timestamps;
MUA.channels = channels;
MUA.samplingRate = lfp.samplingRate;

if saveMat
    save(savefile,'MUA','-v7.3');
end



%%
if SHOWFIG
    showwin = bz_RandomWindowInIntervals(MUALFP.timestamps([1 end]),10);
    timewin = MUALFP.timestamps >showwin(1) & MUALFP.timestamps <showwin(2);
    timewinMUA = MUA.timestamps >showwin(1) & MUA.timestamps <showwin(2);
    
    figure
    subplot(5,1,(4:5))
    plot(MUALFP.timestamps(timewin),MUALFP.data(timewin),'k')
    hold on
    plot(MUALFP.timestamps(timewin),MUALFP.amp(timewin),'b')
    plot(MUALFP.timestamps(timewin),MUALFP.smoothamp(timewin),'r')
   % plot(MUALFP.peaktimes,-MUALFP.peakmags,'.')
    %plot(MUAspks.timestamps,MUAspks.data,'g')
    legend(['MUA (',num2str(MUAfilter(1)),'-',num2str(MUAfilter(2)),'Hz)'],...
        'Amplitude',['Smooth: ',num2str(MUAsmoothwin.*1000),'ms'])
    xlim(showwin)


    subplot(5,1,1)
    plot(datlfp.timestamps(timewin),datlfp.data(timewin),'k')
    hold on
    bz_ScaleBar('s')
    xlim(showwin)
    
    subplot(5,1,2:3)
    plot(MUA.timestamps(timewinMUA),MUA.data(timewinMUA,1))
    hold on
    plot(MUA.timestamps(timewinMUA),sum(MUA.data(timewinMUA,:),2))
end

%% Get the EMG
if compareEMG
    [EMGFromLFP] = bz_EMGFromLFP(basePath,'samplingFrequency',10,'fromDat',true);

    %%
    EMGFromLFP.MUA = interp1(MUALFP.timestamps,MUALFP.smoothamp,EMGFromLFP.timestamps);
    smoothcorr = corr(EMGFromLFP.MUA,EMGFromLFP.data)

    if SHOWFIG
    %%
    showwin = bz_RandomWindowInIntervals(MUALFP.timestamps([1 end]),10);
    timewin = MUALFP.timestamps >showwin(1) & MUALFP.timestamps <showwin(2);

    figure
    subplot(4,1,1)
    plot(EMGFromLFP.timestamps,EMGFromLFP.data)
    xlim(showwin)
    ylabel('EMGfromLFP')

    subplot(4,1,2)
    plot(MUALFP.timestamps(timewin),MUALFP.smoothamp(timewin),'r')
    xlim(showwin)
    ylabel('MUA')

    subplot(4,1,3)
    plot(MUALFP.timestamps(timewin),MUALFP.smoothamp_mov200(timewin),'r')
    xlim(showwin)
    ylabel('MUA')

    subplot(4,1,4)
    plot(MUALFP.timestamps(timewin),MUALFP.smoothamp_mov1000(timewin),'r')
    xlim(showwin)
    ylabel('MUA')

    % subplot(3,1,3)
    % plot(MUAspks.timestamps,MUAspks.data,'r')
    % xlim(showwin)
    % ylabel('MUA Binned')

    end 
end
