function sessRippleSpectogram(varargin)

%Plots the spectogram around baseline and stim ripples, for awake and sleep
%separately. 
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotFig',true,@islogical);
addParameter(p,'twin',0.5,@isnumeric)
addParameter(p,'nfreqs',100,@isnumeric);
addParameter(p,'ncyc',7,@isnumeric);
addParameter(p,'fBand',[2 300],@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
plotFig = p.Results.plotFig;
saveMat = p.Results.saveMat;
nfreqs = p.Results.nfreqs;
twin = p.Results.twin;
ncyc = p.Results.ncyc;
fBand = p.Results.fBand;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

for ii = 1:size(allSess,1)
    fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); 
    file = dir(('*.session.mat'));
    load(file.name);                 
    file = dir(('*.pulses.events.mat'));
    load(file.name);         
    file = dir(('*.ripples.events.mat'));
    load(file.name);
    file = dir(('*.hippocampalLayersCSD.channelinfo.mat'));
    load(file.name);  

    %% Load LFP
    if isnan(hippocampalLayers.ml)
        channelOrder = hippocampalLayers.all(1:(end-1));
    else
        channelOrder = hippocampalLayers.all;
    end

    lfp = bz_GetLFP(channelOrder,'noPrompts', true);

    
    %% Get pulses
    for rr = 2%1:3
        keepIdx = [];       
        ripSpect.Pre = [];
        ripSpect.Post = [];
        
        if exist('pulses')
            if rr <= 2
                pulTr = (pulses.stimComb==rr);
            else
                pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
            end
        end      

        %Only select the pulses that happened in the home cage,
        %i.e., which were 5 seconds long
        homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
        homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
        pulTr = pulTr & homeCagePulseidx;
        events = pulses.intsPeriods(:,pulTr)';
        
        keepIdx(:,1) = InIntervals(ripples.peaks,events);    
        keepIdx(:,2) = InIntervals(ripples.peaks,(events-5));   
        
        % Now build a spectogram for each ripple
        tsPost = ripples.peaks(logical(keepIdx(:,1)));
        tsPre = ripples.peaks(logical(keepIdx(:,2)));
        
        if plotFig
           figure
           set(gcf,'Renderer','painters')
           set(gcf,'Color','w')
        end
        
        for ch = 1:length(channelOrder)
            ripSpectPre =[];
            ripSpectPost =[];
        
            for iRip = 1:length(tsPre)
                [~,idx] = min(abs(lfp.timestamps-tsPre(iRip)));
                intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(twin*1250))];
                wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals,'chanID',channelOrder(ch));
                ripSpectPre(:,:,iRip) = log10(abs(wavelet.data)); 
                freqs = wavelet.freqs;                
            end

            for iRip = 1:length(tsPost)
                [~,idx] = min(abs(lfp.timestamps-tsPost(iRip)));
                intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(twin*1250))];
                wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals,'chanID',channelOrder(ch));
                ripSpectPost(:,:,iRip) = log10(abs(wavelet.data)); 
                freqs = wavelet.freqs;                        
            end
            ripSpect.Pre(:,:,ch) = nanmedian(ripSpectPre,3);
            ripSpect.Post(:,:,ch) = nanmedian(ripSpectPost,3);           

            if plotFig
                subplot(length(channelOrder),2,(ch*2)-1)
                imagesc(1:1:size(ripSpect.Pre,1),freqs,ripSpect.Pre(:,:,ch)')
                set(gca,'YDir','normal')
                set(gca,'YScale','log')
                ylim([2 300])
                caxis([1 4])
                
                subplot(length(channelOrder),2,(ch*2))
                imagesc(1:1:size(ripSpect.Post,1),freqs,ripSpect.Post(:,:,ch)')
                set(gca,'YDir','normal')                
                set(gca,'YScale','log')
                ylim([2 300])
                caxis([1 4])
            end            
        end
        
    end
end

end

