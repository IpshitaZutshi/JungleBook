function compileMiceSleepScore(varargin)

%% Defaults and parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',true,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
savePlot = p.Results.savePlot;

tag = 'mEC'; 

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ29\Final','IZ30\Final','IZ31\Final','IZ32\Final','IZ33\Final','IZ34\Final'};  % To add: IZ23
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'}; % To add: IZ35
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add: IZ35
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    %
    reg = {'contramEC','ipsimEC','Both'};
end

tar = {'STEM', 'RETURN'};

SSRatio  = [];


for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));

        file = dir(['*.MergePoints.events.mat']);
        load(file.name);        
        
        file = dir(['*.SleepState.states.mat']);
        load(file.name);  
        
        if strcmp(mice{m},'IZ27\Final')==1 || strcmp(mice{m},'IZ28\Final')==1 || strcmp(mice{m},'IZ29\Final')==1 || strcmp(mice{m},'IZ32\Final')==1 || ...
                strcmp(mice{m},'IZ33\Final')==1 || strcmp(mice{m},'IZ34\Final')==1 || strcmp(mice{m},'IZ27\Saline')==1 || strcmp(mice{m},'IZ28\Saline')==1 
            homecageDur = [MergePoints.timestamps(2,1) MergePoints.timestamps(3,2)];
        else
            homecageDur = [MergePoints.timestamps(2,1) MergePoints.timestamps(2,2)];
        end
        
        %Define proportion of time in each sleep state
        timeIdx = SleepState.idx.timestamps>homecageDur(1) & SleepState.idx.timestamps<homecageDur(2);
        WAKE = sum(SleepState.idx.states(timeIdx)==1)./sum(timeIdx);
        NREM = sum(SleepState.idx.states(timeIdx)==3)./sum(timeIdx);
        REM = sum(SleepState.idx.states(timeIdx)==5)./sum(timeIdx);
        
        SSRatio  = [SSRatio; WAKE NREM REM homecageDur(2)-homecageDur(1)];
    end
end

dataId = [ones(size(SSRatio,1),1)*1; ones(size(SSRatio,1),1)*2;ones(size(SSRatio,1),1)*3];
data = [SSRatio(:,1);SSRatio(:,2);SSRatio(:,3)];

figure
subplot(1,2,1)
groupStats(data,dataId,'inAxis',true)

subplot(1,2,2)
boxplot(SSRatio(:,4)/3600)
% 
% if savePlot
%     saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.png'));
%     saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.fig'));
%     saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.eps'),'epsc');
%     save(strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'Stats.mat'),'stats');
% end

end