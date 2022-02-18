function SessScriptWrapper(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice = {'IZ27\Final','IZ27\Saline','IZ28\Final', 'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ32\Final','IZ32\Saline',...
         'IZ33\Final','IZ33\Saline','IZ34\Final','IZ34\Saline'};

% {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
%     'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ30\Final','IZ31\Final'};

%    'IZ27\Final','IZ27\Saline','IZ28\Final', 'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ32\Final','IZ32\Saline',...
%        'IZ33\Final','IZ33\Saline','IZ34\Final','IZ34\Saline'
%CA1 mice = 'IZ18\Final','IZ20\Final','IZ21\Final','IZ31\Final'

%'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
%        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ27\Saline','IZ28\Final',...
%        'IZ28\Saline','IZ29\Final','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Final','IZ32\Saline',...
%        'IZ33\Final','IZ33\Saline','IZ34\Final','IZ34\Saline'};

for m = 1:length(mice)
    
     %cd(strcat(parentDir, mice{m}));
     %SessBehaviorPowerSpectrumCSD('expPath',strcat(parentDir, mice{m}),'force',true,'isCA3',false)
     %SessBehaviorCSDUpdateAvg('expPath',strcat(parentDir, mice{m}),'makePlots',false,'plotfig',false);
     %SessPeriStimModIndexCSDradiatum('expPath',strcat(parentDir, mice{m}),'force',true);
     %getPlaceFieldBoundaries('expPath',strcat(parentDir, mice{m}),'downsample',true);
     %SessTrialsPlaceFields('expPath',strcat(parentDir, mice{m}),'downsample',true);
     %SessTriallByTrialMaps('expPath',strcat(parentDir, mice{m}),'downsample',true);
     %SessBehaviorThetaCompression('expPath',strcat(parentDir, mice{m}));
     %SessPeriStimACFrequency('expPath',strcat(parentDir, mice{m}));
     %getTrackingAcrossSess('expPath',strcat(parentDir, mice{m}));
     %SessPeriStimPhaseAmpCSD('expPath',strcat(parentDir, mice{m}));
     %SessBehaviorCoherenceCSD('expPath',strcat(parentDir, mice{m}));
     %SessBehaviorAssembliesSeparate('expPath',strcat(parentDir, mice{m}));
     SessBehaviorCorrMatrixCA3('expPath',strcat(parentDir, mice{m}));
     %SessPeriStimPhaseLocking('expPath',strcat(parentDir, mice{m}),'force',true);
     %SessBehaviorPowerSpectrum('expPath',pwd,'force',false);
     %SessBehaviorLFPSpeed('expPath',pwd);
     %SessBehaviorCSDCoherence('expPath',strcat(parentDir, mice{m}),'force',true);
     %SessBehaviorCSD('expPath',strcat(parentDir, mice{m}));
     %calculateICAwrapper('expPath',strcat(parentDir, mice{m}),'limitTime',false,'force',true)
     %SessUnitISI('expPath',strcat(parentDir, mice{m})
     %UpdateSessPhaseLocking('expPath',strcat(parentDir, mice{m}))
     %SessBehaviorPhasePrecession('expPath',strcat(parentDir, mice{m}))
     %SessRippleCSD('expPath',strcat(parentDir, mice{m}));
end

end