function compileSessiontemplate(varargin)

% expPath = uigetdir; % select folder
% allpath = strsplit(genpath(expPath),';'); % all folders
% cd(allpath{1});
% allSess = dir('*_sess*');

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\Auditory_Task\',@isfolder);
parse(p,varargin{:});
parentDir = p.Results.parentDir;

mice = {'IZ41\Final'};


% {'IZ11\Final','IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ18\DG_CA3','IZ20\Final','IZ21\Final','IZ23\Final',...
%     'IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ27\Saline', 'IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final'};

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    for ii = 2:size(allSess,1)
        fprintf(' ** Analyzing session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
        session = sessionTemplateIZ(pwd,'noPrompts', true,'showGUI',false);

        %Assign brain regions
        %BuzLin probe 
        if strcmp(mice{m},'IZ18\DG_CA3') == 1 
            if ii == 3
                session.brainRegions.CA1.channels = []; 
                session.brainRegions.CA3.channels = [sessionInfo.AnatGrps(1).Channels(1:21) sessionInfo.AnatGrps(2).Channels(1:21)]+1;
                session.brainRegions.DG.channels = [];  
            else
                session.brainRegions.CA1.channels = []; 
                session.brainRegions.CA3.channels = 1:(sessionInfo.nChannels-5); 
                session.brainRegions.DG.channels = [];                                        
            end                    
        elseif strcmp(mice{m},'IZ20\Final') == 1 
            session.brainRegions.CA1.channels = sessionInfo.AnatGrps(1).Channels(1:32)+1; 
            session.brainRegions.CA3.channels = [];
            session.brainRegions.DG.channels = sessionInfo.AnatGrps(1).Channels(33:end)+1;                           
        elseif strcmp(mice{m},'IZ21\Final') == 1                 
            session.brainRegions.CA1.channels = sessionInfo.AnatGrps(1).Channels(1:32)+1; 
            session.brainRegions.CA3.channels = [];
            session.brainRegions.DG.channels = sessionInfo.AnatGrps(1).Channels(33:end)+1;                   
        elseif strcmp(mice{m},'IZ26\Final') == 1 
            session.brainRegions.CA1.channels = [sessionInfo.AnatGrps(1).Channels sessionInfo.AnatGrps(2).Channels sessionInfo.AnatGrps(3).Channels(1:32)]+1;
            session.brainRegions.CA3.channels = [];
            session.brainRegions.DG.channels = sessionInfo.AnatGrps(3).Channels(33:end)+1;  
        elseif strcmp(mice{m},'IZ27\Final') == 1 || strcmp(mice{m},'IZ27\Saline') == 1
            session.brainRegions.CA1.channels = sessionInfo.AnatGrps(1).Channels(1:32)+1; 
            session.brainRegions.CA3.channels = [];
            session.brainRegions.DG.channels = sessionInfo.AnatGrps(1).Channels(33:end)+1;                    
        elseif strcmp(mice{m},'IZ28\Final') == 1 || strcmp(mice{m},'IZ28\Saline') == 1       
            session.brainRegions.CA1.channels = sessionInfo.AnatGrps(1).Channels(1:32)+1; 
            session.brainRegions.CA3.channels = sessionInfo.AnatGrps(1).Channels(33:end)+1;   
            session.brainRegions.DG.channels = [];     
        elseif strcmp(mice{m},'IZ33\Final') == 1 || strcmp(mice{m},'IZ33\Saline') == 1       
            session.brainRegions.CA1.channels = [sessionInfo.AnatGrps(1).Channels(1:42) sessionInfo.AnatGrps(2).Channels(1:32)]+1; 
            session.brainRegions.CA3.channels = [sessionInfo.AnatGrps(1).Channels(43:end)+1 sessionInfo.AnatGrps(2).Channels(33:end)+1];    
            session.brainRegions.DG.channels =  [];         
        elseif strcmp(mice{m},'IZ34\Final') == 1 || strcmp(mice{m},'IZ34\Saline') == 1       
            session.brainRegions.CA1.channels = sessionInfo.AnatGrps(2).Channels(1:32)+1; 
            session.brainRegions.CA3.channels = sessionInfo.AnatGrps(2).Channels(33:end)+1;   
            session.brainRegions.DG.channels = [];   
        else % all others are only in CA1
            session.brainRegions.CA1.channels = 1:sessionInfo.nChannels; 
            session.brainRegions.CA3.channels = [];
            session.brainRegions.DG.channels = [];
        end

%         %Assign ripple channel
%         file = dir(('*.hippocampalLayers.channelInfo.mat'));
%         load(file.name);    
%         pyrCh = hippocampalLayers.pyramidal;            
%         session.channelTags.Ripple.channels = pyrCh+1; 
%         if strcmp(mice{m}(3:4),'11') == 1 || strcmp(mice{m}(3:4),'15') == 1 || strcmp(mice{m}(3:4),'23') == 1 || strcmp(mice{m}(3:4),'29') == 1 || strcmp(mice{m}(3:4),'30') == 1
%             session.channelTags.RippleNoise.channels = 29+1;
%         end
% 
%         %Assign bad channels
%         if exist('badChannels.mat','file')
%             load('badChannels.mat')
%             session.channelTags.Bad.channels = badChannels+1;
%         else
%             session.channelTags.Bad.channels = [];
%         end
% 
%         %Assign last group (with analog channels) as bad
%         session.channelTags.Bad.electrodeGroups = numel(sessionInfo.AnatGrps); % Bad spike groups (broken shanks)

        %Other
        session.general.experimenters = 'Ipshita';

        %Save session file
        save([sessionInfo.session.name '.session.mat'],'session');

    end
end
end