% Compile data across all sessions

function compileMiceCSDSinkShift(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'combineSessions',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

tag = 'mEC'; 

mice = {'IZ12\Final','IZ13\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
    'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ28\Final',...
    'IZ31\Final','IZ33\Final'};

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

layerNum = 1;
    
%% Loop through the mice
for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    
    if exist(strcat('Summ\csdData.mat'),'file')
        load(strcat('Summ\csdData.mat'));
    else 
        disp(['CSD not computed for mouse' mice{m}])
        continue;
    end       
    
    % Define this for each mouse
    for rr = 1:3
        for cc = 1:2
           layers{rr,cc} = [];
        end
    end
    channelOrder = [];

    % Have to cycle through each session to pick the channels for each day
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    
    % Start collecting data
    for ss = 1:size(allSess,1)
        cd(strcat(allSess(ss).folder,'\',allSess(ss).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
        
        load([sessionInfo.FileName '.hippocampalLayersCSD.channelinfo.mat']);
        channelOrder = hippocampalLayers.channelOrder;
        
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
            
            layers{region,target} = [layers{region,target}; hippocampalLayers.slm];
        end
    end
    
    for ll = 1:size(layers{2,1},1)
        
        if ~isempty(csdData.orCh{2,1}{2})
            layerInfo = layers{2,1}(ll);
            [~,layerInfoIdx] = ismember(layerInfo,channelOrder);

            startCh = layerInfoIdx-5;
            endCh = min(layerInfoIdx+5,size(csdData.orCh{2,1}{5},2));
            [~,minCSD] = min(csdData.orCh{2,1}{2}(1,startCh:endCh,ll));
            layerInfoIdxBase(layerNum) = startCh+(minCSD-1); 
            
            startCh = layerInfoIdx-5;
            endCh = min(layerInfoIdx+5,size(csdData.orCh{2,1}{5},2));
            [~,minCSD] = min(csdData.orCh{2,1}{5}(1,startCh:endCh,ll));
            layerInfoIdxStim(layerNum) = startCh+(minCSD-1); 
            
            layerNum = layerNum+1;
        end
    end          
    clear layers
end

a = layerInfoIdxBase-layerInfoIdxStim;
b = histcounts(a,[0:1:8]);
figure
bar([0:1:7],b)
xlabel('Diff in CSD peak channel (baseline-stim)')
end
 


