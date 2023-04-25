function basepaths = compileCellMetrics(varargin)

% expPath = uigetdir; % select folder
% allpath = strsplit(genpath(expPath),';'); % all folders
% cd(allpath{1});
% allSess = dir('*_sess*');

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\Auditory_Task\',@isfolder);
addParameter(p,'calcMetrics',true,@islogical);
addParameter(p,'calcHippocampalLayers',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
calcMetrics = p.Results.calcMetrics;
calcHippocampalLayers = p.Results.calcHippocampalLayers;

mice = {'IZ41\Final\'};
% mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
%        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ27\Saline','IZ28\Final',...
%        'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Final','IZ32\Saline',...
%        'IZ33\Final','IZ33\Saline','IZ34\Final','IZ34\Saline'};


basepaths = [];
kk = 1;

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    for ii = 1:size(allSess,1)
        fprintf(' ** Analyzing session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        basepaths{kk} = strcat(allSess(ii).folder,'\',allSess(ii).name);
        kk = kk+1;
        if calcMetrics   
            [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix;
           % cell_metrics = ProcessCellMetrics('manualAdjustMonoSyn',false,'forceReload',true,'submitToDatabase',false);
           % close all
        end
%         if calcHippocampalLayers 
%           %  try
%                 getHippocampalLayers_IZ;
%            % end
%             close all
%         end
    end      
end
end