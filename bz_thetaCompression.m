function [thetaCompression] = bz_thetaCompression(spikes,lfp,varargin)

% USAGE
%[thetaCompression] = bz_thetaCompression(spikes,lfp,varargin)
% 
% INPUTS
% spikes        -spike time cellinfo struct
% lfp           -lfp struct with a single channel (to determine theta phase) from bz_GetLFP()
% fieldLoc      -(Optional) Array of place field locations for cells with size equal to spikes
% fieldSize      -(Optional) Array of place field sizes for cells with size equal to spikes
% passband      -frequency range for theta phase [lowHz highHz] form
% intervals     -(optional) may specify timespans over which to calculate 
%               autocorrelation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
% samplingRate  -specifies lfp sampling frequency default=1250
% method        -method selection for how to generate phase, 
%               possibilties are: 'hilbert' (default) or 'wavelet'
% forceReload   -Overwrite pre-calculated data, default false
% plotfig       -Plot summary figure
% savefig       -Save summary figure
% saveMat       -logical to save cellinfo .mat file with results, default=true
%
%
% OUTPUTS
%
% pairIDs                   - nx2 matrix of pair id combinations, where n = (CellnumC2)
% placefield_difference     - nx1 if place field locations are provided,
%                               physical distance between each cell pair
% placefield_size           - n x 2 matrix of the field size of each pair
% ccgs_place_offset         - offset in the long time scale CCG - place
%                           field time offset (in seconds)
% ccgs_time_offset          - offset in the theta time scale CCG (in ms)
% ccgs_phase_offset         - offset in the theta phase space (in radians)
% ccgs_place                - smoothed ccgs in the long time scale
% ccgs_time                 - smoothed ccgs in the theta time scale
% ccgs_phase                - smoothed ccgs in the theta phase scale
% t_place                   - CCG time lags for long time scale
% t_time                    - CCG time lags for theta time scale
% t_phase                   - CCG time lags for theta phase scale
% compression               - 1 x 6 matrix, linear fit of ccgs_time_offset vs
%                           ccgs_place_offset, and ccgs_time_offset vs place field diff [Slope, R, p, Slope, R, p]

% Calculates CCG across cells in the theta timescale, place field
% timescale, and theta phase
%
% Ipshita Zutshi 2021

% NOTE: wavelet estimation gives some errors! >>

warning off
%% defaults
p = inputParser;
% addRequired(p,'spikes',[],@isstruct);
% addParameter(p,'lfp',[],@isstruct)
addParameter(p,'passband',[6 12],@isnumeric)
addParameter(p,'intervals',[0 inf],@isnumeric)
addParameter(p,'samplingRate',1250,@isnumeric)
addParameter(p,'fieldLoc',[],@isnumeric)
addParameter(p,'fieldSize',[],@isnumeric)
addParameter(p,'method','hilbert',@isstr)
addParameter(p,'forceReload',true,@islogical)
addParameter(p,'plotfig',false,@islogical)
addParameter(p,'savefig',true,@islogical)
addParameter(p,'saveMat',false,@islogical)

%for IZ
addParameter(p,'zonenum',[],@isnumeric)

parse(p,varargin{:})

passband = p.Results.passband;
intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
fieldLoc = p.Results.fieldLoc;
fieldSize = p.Results.fieldSize;
method = p.Results.method; 
forceReload = p.Results.forceReload;
plotfig = p.Results.plotfig;
savefig = p.Results.savefig;
saveMat = p.Results.saveMat;
%for IZ
zonenum = p.Results.zonenum;

if exist('*.thetaCompression.cellinfo.mat','file') && forceReload == false
    disp('loading theta compression data from cellinfo file..')
    load([spikes.sessionName '.thetaCompression.cellinfo.mat'])    
else

    %% Get phase for every time point in LFP
    switch lower(method)
        case ('hilbert')
            [b, a] = butter(3,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass'); % order 3
            filt = FiltFiltM(b,a,double(lfp.data(:,1)));
            power = fastrms(filt,ceil(samplingRate./passband(1)));  % approximate power is frequency band
            hilb = hilbert(filt);
            lfpphase = mod(angle(hilb),2*pi);
            clear fil
        case ('wavelet')% Use Wavelet transform to calulate the signal phases
            [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(double(lfp.data(:,1)),samplingRate,passband(1),passband(2),8,0);
            [~,mIdx]=max(wave);%get index max power for each timepiont
            pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
            lfpphases=wphases(pIdx);%get phase of max amplitude wave at each timepoint
            lfpphase = mod(lfpphases,2*pi);%covert to 0-2pi rather than -pi:pi
    end


    %% Get phases for each spike for each cell   
    pSp = cell(1,length(spikes.times));
    tSp = cell(1,length(spikes.times));
    
    for a = 1:length(spikes.times)    
        bools = InIntervals(spikes.times{a},intervals);
        tSp{a} =spikes.times{a}(bools);
        [tSp{a}, I] = sort(tSp{a});

        if isempty(tSp{a}) 
            pSp{a} = [];            
        else
            pSp{a} = lfpphase(ceil(tSp{a}*samplingRate));
            pSp{a} = unwrap(pSp{a});
            pSp{a} = pSp{a}(I);
        end
    end
    
    % Calculate peak-peak estimation between place fields
    [ccg2,t2] = CCG(tSp,[],'binSize',0.005,'duration',2);
    t2 = t2*1000;
    
    % Theta-cycle time difference between place fields
    [ccg,t] = CCG(tSp,[],'binSize',0.001,'duration',1.2);
    %[ccg,t] = CCG(tSp,[],'binSize',0.01,'duration',2);
    t = t * 1000; % from sec to msec

    % Phase space 
    [ccg3,t3] = CCG(pSp,[],'binSize',pi/50,'duration',10*pi);    
    %[ccg3,t3] = CCG(pSp,[],'binSize',pi/50,'duration',2*pi);    
    
    %% Initialize    
    pairIDs = NaN(nchoosek(length(spikes.times),2),2);
    
    placefield_difference = NaN(nchoosek(length(spikes.times),2),1);
    placefield_size = NaN(nchoosek(length(spikes.times),2),2);
    placefield_center = NaN(nchoosek(length(spikes.times),2),1);
    
    ccgs_place_offset = NaN(nchoosek(length(spikes.times),2),1);
    ccgs_time_offset = NaN(nchoosek(length(spikes.times),2),1);
    ccgs_phase_offset = NaN(nchoosek(length(spikes.times),2),1);

    ccgs_place = NaN(nchoosek(length(spikes.times),2),length(t2));
    ccgs_time = NaN(nchoosek(length(spikes.times),2),length(t));
    ccgs_phase = NaN(nchoosek(length(spikes.times),2),length(t3));
    
    compression = NaN(1,14);
    ccgs_time_peaks = NaN(nchoosek(length(spikes.times),2),12);
    ccgs_phase_peaks = NaN(nchoosek(length(spikes.times),2),6);
    
    t_place = t2;
    t_time = t;
    t_phase = t3;
    
    kk = 1;
    
    %% Get pairwise ccg information
    for i = 1:(size(ccg,3)-1)
        for j = (i+1):size(ccg,3)
            
            ccg_trace = nanconv(ccg(:,i,j)'./max(ccg(:,i,j)),gausswin(80)','edge'); % fast timescale 
            ccg_trace2 = nanconv(ccg2(:,i,j)'./max(ccg2(:,i,j)),gausswin(160)','edge'); % field peak-peak timescale
            ccg_trace3 = nanconv(ccg3(:,i,j)'./max(ccg3(:,i,j)),gausswin(80)','edge'); % theta phase timescale
            
            pairIDs(kk,:) = [i j];
            
            if ~isempty(fieldLoc)
                placefield_difference(kk) = fieldLoc(j) - fieldLoc(i); 
%                 % only for IZ_ remove later
%                 if fieldLoc(i) > 55 && fieldLoc(j)> 55
%                     placefield_center(kk) = 1;
%                 else
%                     placefield_center(kk) = 0;
%                 end
            end
            
            if ~isempty(fieldSize)
                placefield_size(kk,:) = [fieldSize(j) fieldSize(i)]; 
            else
                placefield_size(kk,:) = nan;
            end
            
            % Get place field time offset
            [~,indx] = max(ccg_trace2);
            if max(ccg_trace2)> (mean(ccg_trace2)+(1.5*std(ccg_trace2)))
                ccgs_place_offset(kk) = t2(indx);
            end
            ccgs_place(kk,:) = ccg_trace2;
            
            % Get theta-cycle time offset
            [~,locs] = findpeaks(ccg_trace,'MinPeakHeight',mean(ccg_trace),'MinPeakDistance',100);
            if isempty(locs)~=1
                [~,J]=sort(abs(t(locs)),'ascend');
                ccgs_time_offset(kk)=t(locs(J(1)));
            end   
            
            ccgs_time_peaks(kk,1:length(locs)) = t(locs);             
            ccgs_time(kk,:) = ccg_trace;             

           
            % Phase space
            [~,locs] = findpeaks(ccg_trace3,'MinPeakHeight',mean(ccg_trace3),'MinPeakDistance',90);
            if isempty(locs)~=1
                [~,J]=sort(abs(t3(locs)),'ascend');
                if abs(t3(locs(J(1)))) < 1.5*pi
                    ccgs_phase_offset(kk)=t3(locs(J(1)));
                end
            end            
            ccgs_phase_peaks(kk,1:length(locs)) = t3(locs);   
            ccgs_phase(kk,:) = ccg_trace3; 
            
            kk = kk+1;
        end
    end
end


%% for IZ only! fitting the slope across all peaks 

% Fit for theta phase
idxtoKeep = placefield_center==1 & ...
        (ccgs_place_offset>-600 & ccgs_place_offset<600);   
StartPoint = [0.001,5,400,20,0.00001,0,0];
LowerLimits = [0.00005,0,100,5,0,-1,-10];
UpperLimits = [0.003,10,700,200,0.01,1,100];
data1 = repmat(ccgs_place_offset(idxtoKeep),[1 size(ccgs_phase_peaks,2)]);
data2 = ccgs_phase_peaks(idxtoKeep,:);
data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
x_out = data1;
y_out = data2;  
x_bins = -400:8:400;
y_bins = -20:0.2:20;
if sum(idxtoKeep) > 5
    [fit_params_phase,~] = customFit(x_out,y_out,x_bins,y_bins,StartPoint,LowerLimits,UpperLimits);
else
    fit_params_phase(1:7) = nan;
end

% Fit for theta timescale
StartPoint =  [0.2,   100,  400,   800,  0.00001,  0.6,     0];
LowerLimits = [0.001,   50,  100,   400, 0.000001, 0.2,  -500];
UpperLimits = [0.8,  150,  700,  1500,    0.001,  1.5,   500];
data1 = repmat(ccgs_place_offset(idxtoKeep),[1 size(ccgs_time_peaks,2)]);
data2 = ccgs_time_peaks(idxtoKeep,:);
data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
x_out = data1;
y_out = data2;    
x_bins = -400:8:400;
y_bins = -490:10:490; 
if sum(idxtoKeep) > 5
    [fit_params_time,~] = customFit(x_out,y_out,x_bins,y_bins,StartPoint,LowerLimits,UpperLimits);    
else
    fit_params_time(1:7) = nan;
end
if plotfig
    % for IZ only!
    idxtoKeep =  placefield_center==1 & ...
        (ccgs_place_offset>-600 & ccgs_place_offset<600);
    
    %Plot fig for inspection
    figure
    data1 = repmat(ccgs_place_offset(idxtoKeep),[1 size(ccgs_time_peaks,2)]);
    data2 = ccgs_time_peaks(idxtoKeep,:);
    data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
    data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
    
    set(gcf,'Position',[600 600 600 300])
    subplot(2,3,[1 2])
    s = scatter(data2,data1,2,'filled');
    s.MarkerFaceAlpha = 0.4;
    xlabel('Theta timescale (ms)')
    ylabel('Place field dist (ms)')
    x_bins = -400:8:400;
    hold on
    for i = -5:5
        plot(x_bins*fit_params_time(1)+i*fit_params_time(2),x_bins,'--r')
    end
    ylim([-400,400]),xlim([-400,400])    
    title(strcat('Slope:',num2str(fit_params_time(1)),' Period:',num2str(fit_params_time(2))))
    
    data1 = repmat(ccgs_place_offset(idxtoKeep),[1 size(ccgs_phase_peaks,2)]);
    data2 = ccgs_phase_peaks(idxtoKeep,:);
    data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
    data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
    data2(data2==0) = nan;
    subplot(2,3,[4 5])
    s = scatter(data2,data1,2,'filled');
    s.MarkerFaceAlpha = 0.4;
    xlabel('Theta phase (rad)')
    ylabel('Place field dist (ms)')
    x_bins = -400:8:400;
    hold on
    for i = -5:5
        plot(x_bins*fit_params_phase(1)+i*fit_params_phase(2),x_bins,'--r')
    end
    ylim([-400,400]),xlim([-16,16])    
    title(strcat('Slope:',num2str(fit_params_phase(1)),' Period:',num2str(fit_params_phase(2))))
          
    subplot(2,3,[3 6])
    scatter(placefield_difference(idxtoKeep)*1.75,ccgs_place_offset(idxtoKeep),25,'.')
    if sum(idxtoKeep) > 10
        [R,pVal] = corr(placefield_difference(idxtoKeep)*1.75,ccgs_place_offset(idxtoKeep),'Rows','pairwise','type','Spearman');
        lsline
        title(strcat('R:',num2str(R),' p:',num2str(pVal)))    
    end
    xlabel('Place field dist (cm)')
    ylabel('Place field dist (ms)')

    if savefig
       saveas(gcf,strcat('thetaCompression',num2str(zonenum),'.png'));
       saveas(gcf,strcat('thetaCompression',num2str(zonenum),'.fig'));              
    end
end

close all

thetaCompression.pairIDs = pairIDs;

thetaCompression.placefield_difference = placefield_difference;
thetaCompression.placefield_size = placefield_size;
thetaCompression.placefield_center = placefield_center;

thetaCompression.ccgs_place_offset = ccgs_place_offset;
thetaCompression.ccgs_time_offset = ccgs_time_offset ;
thetaCompression.ccgs_phase_offset = ccgs_phase_offset;

thetaCompression.ccgs_place = ccgs_place;
thetaCompression.ccgs_time = ccgs_time;
thetaCompression.ccgs_phase = ccgs_phase;

thetaCompression.ccgs_time_peaks = ccgs_time_peaks;
thetaCompression.ccgs_phase_peaks = ccgs_phase_peaks;

thetaCompression.t_place = t2;
thetaCompression.t_time = t;
thetaCompression.t_phase = t3;

thetaCompression.compression = [fit_params_time fit_params_phase];

if saveMat
    save([spikes.sessionName '.thetaCompression.cellinfo.mat'],'thetaCompression');    
end

end
