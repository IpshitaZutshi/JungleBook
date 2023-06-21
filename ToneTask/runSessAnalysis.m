% sess = {'IZ41\Final\IZ41_220624_sess5','IZ41\Final\IZ41_220629_sess7','IZ41\Final\IZ41_220701_sess9',...
%     'IZ41\Final\IZ41_220704_sess11','IZ41\Final\IZ41_220708_sess13','IZ41\Final\IZ41_220714_sess14','IZ46\IZ46_230406_sess9',...
%     'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
%     'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
%     'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
%     'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
%     'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
%     'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
%     'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
%     'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
%     'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
%     'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
%     'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
%     'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
%     'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
%     }; 

sess = {'IZ45\IZ45_230420_sess19'};
    
% ripChan = [6,9,6,60,52,54,33,55,53,55,57,58,52,58,58,7,54,84];
% badChan = [{},{},{},{},{},{36,59,75,74,81},{36,59,75,74,81},{},{},{110,109,108,107,31,19},{110,109,108,107,31,19},{},{},{},{},...
% {110,109,108,107,31,19},{110,109,108,107,31,19},{36,59,75,74,81}];
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

force = 1;
forceRun = 1;
for ii = 1:length(sess)
   
    cd(strcat(expPath,sess{ii}))
    disp(strcat(expPath,sess{ii}))
    basepath = pwd;    
%     
%     file = dir([basepath filesep '*.session.mat']);
%     load(file.name);     
%     session.channelTags.Ripple.channels = ripChan(ii);
%     session.channelTags.Bad.channels = badChan{ii};
%     session.general.basePath = basepath;
%     session.extracellular.chanCoords.layout = 'staggered';
%     C = strsplit(pwd,'\');
%     save([basepath filesep C{end} '.session.mat'],'session');    
    
%       try
%           if isempty(dir('*.deepSuperficialfromRipple.channelinfo.mat'))
%             classification_DeepSuperficial(session);
%           end
%       catch
%          disp(strcat(expPath,sess{ii}))
%       end
%     file = dir([basepath filesep '*.Tracking.Behavior.mat']);
%     load(file.name);    
%     % Get instantaneous speed
%     [~,~,~,vx,vy,~,~] = KalmanVel(tracking.position.x,tracking.position.y,tracking.timestamps,2);
%     v = sqrt(vx.^2+vy.^2);
%     tracking.position.vx  = vx;
%     tracking.position.vy  = vy;
%     tracking.position.v  = v;
%     C = strsplit(pwd,'\');
%     save([basepath filesep C{end} '.Tracking.Behavior.mat'],'tracking');
%         
    cell_metrics = ProcessCellMetrics('manualAdjustMonoSyn',false,'forceReload',true,'submitToDatabase',false,'showGUI',false);
    
     [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');    
      readjustPositionToneTrack;
     
   if isempty(dir('*.rateMapsTrial.cellinfo.mat'))|| force
       disp('Calculating trial by trial rate maps');
       getToneMaps7Ports('plotfig',true)
   end
%     basepath = pwd;
%     file = dir([basepath filesep '*.TrialBehavior.Behavior.mat']);
%     load(file.name);
%     if ~isfield(behavTrials,'stim')
%         behavTrials.stim = zeros(size(behavTrials.timestamps,1),1);
%         C = strsplit(pwd,'\');
%         save([basepath filesep C{end} '.TrialBehavior.Behavior.mat'],'behavTrials');
%     end
    
    if isempty(dir('*.rateMapsAvg.cellinfo.mat'))|| forceRun
        disp('Calculating average rate maps');
        getToneMapsAvg7Ports('plotfig',true)
    end
%    
    digitalIn = bz_getDigitalIn('all','force',true);
    plotting7ports
%     if isempty(dir('*sessionData.mat'))|| force
%         disp('Calculating data matrix');
%         [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix;
%         sum(eventVar.licks,2)
% %     end
% 
    if isempty(dir('*.lfp'))
        try 
            bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
        catch
            disp('Problems with bz_LFPfromDat, resampling...');
            ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
                sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
        end
    end
end

    