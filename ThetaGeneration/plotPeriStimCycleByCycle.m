lfpCh = 32;%[21 17 16 19]; %[18 31 39 34]for IZ4, [61 60 53 50][59 63 53 37] for IZ5

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

if ~exist('analogEv','var')
    analogEv = 64;
    numAnalog = 2;

    if ~isempty(analogEv)
        for ii = 1:numAnalog
            analogCh(ii) = (analogEv-1)+ii;
        end
    end
end

lfp = bz_GetLFP(lfpCh(1),'noPrompts', true);
samplingRate = lfp.samplingRate;
t = lfp.timestamps;

[shape_features] = bz_cycleBycycle('channels',lfpCh);


[pulses] = bz_getAnalogPulses('analogCh',analogCh);

figure   
for i =1:3
    if i<=numAnalog
        pulTr = (pulses.stimComb==i);
    else
        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
    end
    events = pulses.intsPeriods(1,pulTr);
    %events = round(events*samplingRate);
    events = events(((events + 5) <= max(t)) & ((events - 5) > 0));
    
    voltage.Pre = [];
    voltage.Post = [];
    time.Pre = [];
    time.Post = [];
    risedecay.Pre = [];
    risedecay.Post = [];
    peaktrough.Pre = [];
    peaktrough.Post = [];
    
    idxes = shape_features.sample_peak(1:(end-1));
    for pp = 1:length(events)
        
        events_tmp = abs(t-events(pp));
        [~,idx_temp] = min(events_tmp);
        voltage.Pre = [voltage.Pre; shape_features.volt_amp((idxes<idx_temp-0.5*samplingRate) & idxes>(idx_temp-5*samplingRate)&(shape_features.period > 90)& (shape_features.period < 170),:)];
        voltage.Post = [voltage.Post; shape_features.volt_amp((idxes>idx_temp+0.5*samplingRate) & idxes<(idx_temp+5*samplingRate) &(shape_features.period > 90)& (shape_features.period < 170),:)];
        
        time.Pre = [time.Pre; shape_features.period((idxes<idx_temp-0.5*samplingRate) & idxes>(idx_temp-5*samplingRate)&(shape_features.period > 90)& (shape_features.period < 170),:)];
        time.Post = [time.Post; shape_features.period((idxes>idx_temp+0.5*samplingRate) & idxes<(idx_temp+5*samplingRate)&(shape_features.period > 90)& (shape_features.period < 170),:)];
        
        risedecay.Pre = [risedecay.Pre; shape_features.time_rdsym((idxes<idx_temp-0.5*samplingRate) & idxes>(idx_temp-5*samplingRate) &(shape_features.period > 90)& (shape_features.period < 170),:)];
        risedecay.Post = [risedecay.Post; shape_features.time_rdsym((idxes>idx_temp+0.5*samplingRate) & idxes<(idx_temp+5*samplingRate) &(shape_features.period > 90)& (shape_features.period < 170),:)];
        
        peaktrough.Pre = [peaktrough.Pre; shape_features.time_ptsym((idxes<idx_temp-0.5*samplingRate) & idxes>(idx_temp-5*samplingRate)&(shape_features.period > 90)& (shape_features.period < 170),:)];
        peaktrough.Post = [peaktrough.Post; shape_features.time_ptsym((idxes>idx_temp+0.5*samplingRate) & idxes<(idx_temp+5*samplingRate)&(shape_features.period > 90)& (shape_features.period < 170),:)];
    end
    
    
    subplot(3,3,(3*(i-1)+1))
    nhist(voltage,'binfactor',2,'linewidth',1,'median')
%     
%     subplot(3,4,(4*(i-1)+2))
%     nhist(time,'binfactor',2,'linewidth',1)
    
    subplot(3,3,(3*(i-1)+2))
    nhist(risedecay,'binfactor',2,'linewidth',1,'median')
    
    subplot(3,3,(3*(i-1)+3))
    nhist(peaktrough,'binfactor',2,'linewidth',1,'median')
end
    
    