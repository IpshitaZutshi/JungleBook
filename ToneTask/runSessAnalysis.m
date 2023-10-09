
% sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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

sess = {'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
%expPath = 'C:\Data\Kilosort\';

force = 1;
forceRun = 1;

parfor ii = 1:length(sess)
   
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
% %     % Get instantaneous speed
%     [~,~,~,vx,vy,~,~] = KalmanVel(tracking.position.x,tracking.position.y,tracking.timestamps,2);
%     v = sqrt(vx.^2+vy.^2);
%     tracking.position.vx  = vx;
%     tracking.position.vy  = vy;
%     tracking.position.v  = v;
%     C = strsplit(pwd,'\');
%     save([basepath filesep C{end} '.Tracking.Behavior.mat'],'tracking');
% % %         
%       cell_metrics = ProcessCellMetrics('manualAdjustMonoSyn',false,'forceReload',true,'submitToDatabase',false,'showGUI',false);
% % %     
%      [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');         
%      getToneTracking;
%     readjustPositionToneTrack;
   
%    if isempty(dir('*.rateMapsTrial.cellinfo.mat'))|| force
%        disp('Calculating trial by trial rate maps');
%        getToneMaps7Ports('plotfig',true)
%    end
%     basepath = pwd;
%     file = dir([basepath filesep '*.TrialBehavior.Behavior.mat']);
%     load(file.name);
%     if ~isfield(behavTrials,'stim')
%         behavTrials.stim = zeros(size(behavTrials.timestamps,1),1);
%         C = strsplit(pwd,'\');
%         save([basepath filesep C{end} '.TrialBehavior.Behavior.mat'],'behavTrials');
%     end
    
%     if isempty(dir('*.rateMapsAvg.cellinfo.mat'))|| forceRun
%         disp('Calculating average rate maps');
%         getToneMapsAvg7Ports('plotfig',false)
%     end
% 
     getToneMapsAvgLickLocProbe('plotfig',false)
%     if isempty(dir('*.rateMapsAvgError.cellinfo.mat'))|| forceRun
%         disp('Calculating average rate maps');
%         getToneMapsAvg7PortsError('plotfig',false)
%     end
    
%     getToneMapsAvg7Ports2D('plotfig',false)
% %    
%     digitalIn = bz_getDigitalIn('all');
%     plotting7ports
% % 
%     if isempty(dir('*.lfp'))
%         try 
%             bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
%         catch
%             disp('Problems with bz_LFPfromDat, resampling...');
%             ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
%                 sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
%         end
%     end
end

    