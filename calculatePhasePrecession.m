function calculatePhasePrecession(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'forceReload',true,@islogical);
parse(p,varargin{:});
parentDir = p.Results.parentDir;
forceReload = p.Results.forceReload;

mice = {'IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ27\Final','IZ20\Final','IZ21\Final',...
        'IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Final',...
       'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Final','IZ32\Saline',...
       'IZ33\Final','IZ33\Saline','IZ34\Final','IZ34\Saline'};
%{'IZ11\Final','IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ18\DG_CA3','IZ20\Final','IZ21\Final','IZ23\Final',...
 %   'IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ27\Saline', 'IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final'};

parfor m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');

    for ii = 1:size(allSess,1)
        fprintf(' ** Analyzing session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        file = dir(['*.SessionPulses.Events.mat']);
        pulses = load(file.name);
        flds = fieldnames(pulses.sessionPulses);
        if ~exist([flds{1} '\' allSess(ii).name '.phasePrecessionAvg.cellinfo.mat'],'file') || forceReload
            getPhasePrecession;        
        end        
    end      
end
end