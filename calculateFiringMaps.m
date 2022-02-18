function calculateFiringMaps(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'forceReload',true,@islogical);
parse(p,varargin{:});
parentDir = p.Results.parentDir;
forceReload = p.Results.forceReload;

mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};%'IZ11\Final','IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final','IZ21\Final',...
   % 'IZ27\Final','IZ27\Saline', 'IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final'};
%'IZ11\Final','IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final',

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');

    for ii = 1:size(allSess,1)
        fprintf(' ** Analyzing session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        file = dir(['*.SessionPulses.Events.mat']);
        load(file.name);
        flds = fieldnames(sessionPulses);
        if ~exist([flds{1} '\' allSess(ii).name '.firingMapsAvg.cellinfo.mat'],'file') || forceReload
            getPlaceFields;        
        end        
    end      
end
end