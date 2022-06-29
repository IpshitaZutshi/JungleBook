function getTrackingAcrossSess(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'isCA3',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
%isCA3 = p.Results.isCA3;

if isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
%chName = {'pyramidal','radiatum','slm','ml'};

for ii = 1:size(allSess,1)
   % try
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        
        autoDSdetectionIZ
%        getPlaceFieldTemplates;
%        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
%         if exist('badChannels.mat','file')
%             load('badChannels.mat')            
%         else 
%             badChannels = [];
%         end
%         if sessionInfo.nChannels<70
%             badChannels = [badChannels 64:1:(sessionInfo.nChannels-1)];
%         else
%             badChannels = [badChannels 128:1:(sessionInfo.nChannels-1)];
%         end
%         SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels,'noPrompts',true); % try to sleep score
%              
%         save([allSess(ii).name '.hippocampalLayers.channelinfo.mat'],'hippocampalLayers');
        %AnalysisBatchTheta;
        %getSessionLinearize('forceReload',true);
%        getPlaceFields;
%       getPlaceFieldsDownsample('isCA3',isCA3);
%         file = dir(('*.hippocampalLayers.channelInfo.mat'));
%         load(file.name);
%         
%         load([sessionInfo.FileName '.session.mat']);
%         if isfield(session.channelTags,'RippleNoise')
%             noiseCh = session.channelTags.RippleNoise.channels-1;
%         else
%             noiseCh = [];
%         end
%             
%         rippleMasterDetectorIZ('rippleChannel',hippocampalLayers.pyramidal,'SWChannel',hippocampalLayers.radiatum,'noiseCh',noiseCh);
%         close all
% %         getLFPduringtrack('refChannel',[],'pyrChPlus',hippocampalLayers.all,'numtrials',15);
%        getPhasePrecession
%         file = dir(('*.region.mat'));
%         load(file.name);   
%          region.CA1sp = 39;
% %        region.CA1so = 31; 
%         save(strcat(allSess(ii).name,'.region.mat'),'region');

%         badChannels = [46 108 116];
%         save('icaChannels.mat','icaChannels');
             
            
%          getSessionTracking('convFact',0.1149,'roiTracking','manual','roiLED','manual','forceReload',true);
%          getSessionArmChoice('task','alternation','force',true);
%          getSessionLinearize('forceReload',true);
%          getSessionPulses('force',true);
%      catch  
% %         fprintf(' ** Error in session %s... \n',allSess(ii).name);
%      end
end