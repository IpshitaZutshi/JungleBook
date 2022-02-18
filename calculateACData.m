function calculateACData(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'forceReload',true,@islogical);
parse(p,varargin{:});
parentDir = p.Results.parentDir;
forceReload = p.Results.forceReload;

mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};%,...
   %'IZ27\Final','IZ27\Saline', 'IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final'};

%{'IZ11\Final','IZ12\Final\','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ18\DG_CA3','IZ20\Final','IZ21\Final','IZ23\Final',...
 %   'IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ27\Saline', 'IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final'};

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    fprintf(' ** Analyzing mouse %3.i of %3.i... \n',m, length(mice));
    SessPeriStimACFrequency('expPath',pwd,'force',true);
end
end