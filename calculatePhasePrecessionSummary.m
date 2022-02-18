function calculatePhasePrecessionSummary(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'forceReload',true,@islogical);
parse(p,varargin{:});
parentDir = p.Results.parentDir;
forceReload = p.Results.forceReload;

mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
       'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ27\Saline','IZ28\Final',...
       'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Final','IZ32\Saline',...
       'IZ33\Final','IZ33\Saline','IZ34\Saline'};%
% {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
%        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ27\Saline','IZ28\Final',...
%        'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Final','IZ32\Saline',...
%        'IZ33\Final','IZ33\Saline','IZ34\Saline'};%
%{'IZ15\Final','IZ18\Final','IZ30\Final','IZ20\Final','IZ31\Final'}
%{'IZ11\Final','IZ12\Final\','IZ13\Final','IZ17\Final','IZ21\Final'}
%{'IZ23\Final','IZ24\Final','IZ25\Final','IZ26\Final'
%{'IZ27\Final','IZ27\Saline','IZ28\Final','IZ28\Saline','IZ29\Final','IZ29\Saline'};

parfor m = 1:length(mice)

    %cd(strcat(parentDir, mice{m}));
    fprintf(' ** Analyzing mouse %3.i of %3.i... \n',m, length(mice));
    SessBehaviorPhasePrecession('expPath',strcat(parentDir, mice{m}),'force',true,'makePlots',false,'savePlots',false);
    close all
end
end