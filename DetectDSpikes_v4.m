function [DS2] = DetectDSpikes_v4(fname,ch_hilus,ch_molecular,res_per,varargin)

% 
% 20150419 Yuta Senzai v3: consider the post-DS1 value in order to get rid of light evoked LFP deflection (for YM33) 
%                          this effect seems to be only in molecular layer. This may be related to
%                          granule cell inhibition in POMC-Cre::Arch animal
%
% 20150602 Yuta Senzai v4: add m_minValue > 
% res_per = default [0 Inf]
%

% default 
show = 'on';

DS2_lowThreshold = 4500; % 3000*0.38uV = 1.14mV 
DS1_lowThreshold = 3000; % 3000*0.38uV = 1.14mV 
DS2fil_highThreshold = 4500;
DS1fil_highThreshold = 1500;
DS1wb_highThreshold = 3000;
DS2_mol_threshold = 500; 

% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
    end
    switch((varargin{i})),
        case 'DS2_lowThreshold',
            DS2_lowThreshold = varargin{i+1};
            if ~isdscalar(DS2_lowThreshold,'>0'),
                error('Incorrect value for property ''DS2_lowThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS1_lowThreshold',
            DS1_lowThreshold = varargin{i+1};
            if ~isdscalar(DS1_lowThreshold,'>0'),
                error('Incorrect value for property ''DS2_lowThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS2fil_highThreshold',
            DS2fil_highThreshold = varargin{i+1};
            if ~isdscalar(DS2fil_highThreshold,'>0'),
                error('Incorrect value for property ''DS2fil_highThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
            
        case 'DS2_mol_threshold',
            DS2_mol_threshold = varargin{i+1};
            if ~isdscalar(DS2_mol_threshold,'>0'),
                error('Incorrect value for property ''DS2_mol_threshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end    

        case 'DS1fil_highThreshold',
            DS1fil_highThreshold = varargin{i+1};
            if ~isdscalar(DS1fil_highThreshold,'>0'),
                error('Incorrect value for property ''DS1fil_highThreshold'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'DS1wb_highThreshold',
            DS1wb_highThreshold = varargin{i+1};
            if ~isdscalar(DS1wb_highThreshold,'>0'),
                error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'show',
            show = varargin{i+1};
            if ~isstring(show,'on','off'),
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
    end
end


%res_lfp_h = GetEEG(ch_hilus,'restrict',res_per);
%res_lfp_m = GetEEG(ch_molecular,'restrict',res_per);
res_lfp_h = bz_GetLFP(ch_hilus,'intervals',res_per,'noPrompts',true);
res_lfp_m = bz_GetLFP(ch_molecular,'intervals',res_per,'noPrompts',true);

fil_res_lfp_h = bz_Filter(res_lfp_h,'passband',[2 50]);
fil_res_lfp_m = bz_Filter(res_lfp_m,'passband',[2 50]);

res_lfp_h = [res_lfp_h.timestamps,double(res_lfp_h.data)];
res_lfp_m = [res_lfp_m.timestamps,double(res_lfp_m.data)];

fil_res_lfp_h = [fil_res_lfp_h.timestamps,double(fil_res_lfp_h.data)];
fil_res_lfp_m = [fil_res_lfp_m.timestamps,double(fil_res_lfp_m.data)];

% fil_res_lfp_h = FilterLFP(res_lfp_h,...
%     'nyquist',625,'passband',[2 50]);
% fil_res_lfp_m = FilterLFP(res_lfp_m,...
%     'nyquist',625,'passband',[2 50]);

time = fil_res_lfp_h(:,1);
hm_dif = fil_res_lfp_h(:,2) - fil_res_lfp_m(:,2);
%time = fil_res_lfp_h.timestamps;
%hm_dif = fil_res_lfp_h.data - fil_res_lfp_m.data;
lfp_dif = [time, hm_dif];


%% DS2 detection 

DS2_thresholded = hm_dif > DS2_lowThreshold;

% Also do a theta power filter

DS2_start = find(diff(DS2_thresholded)>0);
DS2_stop = find(diff(DS2_thresholded)<0);
% Exclude last DS if it is incomplete
if length(DS2_stop) == length(DS2_start)-1,
    DS2_start = DS2_start(1:end-1);
end
% Exclude first DS if it is incomplete
if length(DS2_stop)-1 == length(DS2_start),
    DS2_stop = DS2_stop(2:end);
end
% Correct special case when both first and last DS are incomplete
if DS2_start(1) > DS2_stop(1),
    DS2_stop(1) = [];
    DS2_start(end) = [];
end

DS2_firstPass = [DS2_start,DS2_stop];

% Take out DS which are too close to file beginning- or end
DS2_firstPass = DS2_firstPass(DS2_start > 46 & DS2_stop < size(res_lfp_h,1)-45,:);

if isempty(DS2_firstPass),
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(DS2_firstPass)) ' putative DS2 events.']);
end

DS2_secondPass = DS2_firstPass;


% Discard ripples with a peak power < highThreshold for
DS2 = [];
for i = 1:size(DS2_secondPass,1)
    if DS2_secondPass(i,1)>50
        m_meanValue = mean(res_lfp_m([DS2_secondPass(i,1):DS2_secondPass(i,2)],2));
        ds2criteria = mean(res_lfp_m([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20) ],2))-DS2_mol_threshold; % 25/1250sec = 20ms, 40/1250 = 32 ms
        m_minValue =  min(res_lfp_m([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20), (DS2_secondPass(i,2)+20):(DS2_secondPass(i,2)+45)],2));
        h_maxValue =  max(res_lfp_h([(DS2_secondPass(i,1)-45):(DS2_secondPass(i,1)-20), (DS2_secondPass(i,2)+20):(DS2_secondPass(i,2)+45)],2));
%         if (df_maxValue > DS2fil_highThreshold) && (m_meanValue < ds2criteria)
        if (m_meanValue < ds2criteria && m_minValue > -6000 && h_maxValue < 6000)
%             DS2lfp = Restrict(res_lfp_h,time(DS2_secondPass(i,:))');
            DS2lfp = res_lfp_h( DS2_secondPass(i,1): DS2_secondPass(i,2),:);
            [~,idx]=max(DS2lfp(:,2));
            DS2 = [DS2 ; time(DS2_secondPass(i,1)), DS2lfp(idx,1), time(DS2_secondPass(i,2))];
        end
    end
end

if isempty(DS2),
    disp('DS2 Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: DS2 ' num2str(length(DS2)) ' events.']);
end
DS2triad = DS2;

clear DS2

DS2.timestamps = DS2triad(:,[1 3]);
DS2.peaks = DS2triad(:,2);
DS2peakIDCs = round(DS2triad(:,2)*1250);
DS2.amplitudes = (res_lfp_h(DS2peakIDCs,2)-res_lfp_m(DS2peakIDCs,2))*.00038; % .00038 = conversion to mV
% amplitude: amplitude of each event (Px1).
% amplitudeUnits: specify the units of the amplitude vector.
DS2.amplitudeUnits = 'mV';
% eventID: numeric ID for classifying various event types (Px1).
DS2.eventID = ones(size(DS2triad,1),1);
% eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
DS2.eventIDlabels = repmat({'DS2'},size(DS2triad,1),1);
% eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
DS2.eventIDbinary = false(size(DS2triad,1),1);
% center: center time-point of event (in seconds; calculated from timestamps; Px1).
% duration: duration of event (in seconds; calculated from timestamps; Px1).
DS2.duration = DS2triad(:,2)-DS2triad(:,1);
DS2.center = DS2triad(:,1)+DS2.duration;

DS2.detectorinfo.detectorname = 'DetectDSpikes_v4';
DS2.detectorinfo.detectionparms = [];
DS2.detectorinfo.detectionintervals = [0 Inf];
DS2.detectorinfo.detectiondate = datetime('today');
DS2.detectorinfo.ml_channel = ch_molecular;
DS2.detectorinfo.h_channel = ch_hilus;
DS2.detectorinfo.detectionchannel = ch_hilus;
%save([fname '.DS1.events.mat'],'DS1');
save([fname '.DS2.events.mat'],'DS2');


end