% Compile data across all sessions

function compileMiceRippleCSD(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',false,@islogical);
addParameter(p,'combineSessions',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
savePlot = p.Results.savePlot;
combineSessions = p.Results.combineSessions;

tag = 'CA1'; 

if strcmp(tag,'CA1') == 1
    mice = {'IZ20\Final','IZ21\Final','IZ31\Final','IZ18\Final','IZ15\Final','IZ30\Final'};
    reg = {'CA1','mEC','Both'};% Excluded - 'IZ15\Final','IZ30\Final',
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ17\Final',...
        'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ31\Final','IZ33\Final'}; %'IZ15\Final','IZ29\Saline','IZ32\Final','IZ30\Final',
    reg = {'CA1','mEC','Both'};%
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ33\Final'}; % Excluded - 'IZ29\Final',,'IZ34\Final,'IZ33\Final'
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ33\Saline'}; % Excluded - 'IZ29\Saline','IZ32\Saline',
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'}; 
    reg = {'contramEC','ipsimEC','Both'};
end

target = {'STEM', 'RETURN'};

for ii = 1:length(reg)
    for jj = 1:length(target) 
        csd{ii,jj} = [];
    end
end
    
%% Loop through the mice
for mm = 1:length(mice)
    
    cd(strcat(parentDir, mice{mm},'\Summ'));
    
    if exist(strcat('Summ\csdRippleData.mat'),'file') && strcmp(mice{mm},'IZ26\Final')==1
        load(strcat('Summ\csdRippleData_shank2.mat'));
    elseif exist(strcat('Summ\csdRippleData.mat'),'file')
        load(strcat('Summ\csdRippleData.mat'));
    else 
        disp(['Ripple CSD not computed for mouse' mice{mm}])
        continue;
    end

    % Pick the channels from the last session of recordings
    cd(strcat(parentDir, mice{mm}));
    allSess = dir('*_sess*');
    ss = size(allSess,1);
    
    cd(strcat(allSess(ss).folder,'\',allSess(ss).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
        
    load([sessionInfo.FileName '.hippocampalLayers.channelinfo.mat']);
    channelOrder = hippocampalLayers.channelOrder;

    [~,layerInfoIdx] = ismember(hippocampalLayers.all(2:4),channelOrder);
                        
    if length(channelOrder)>16
        %correct rad channel
        startCh = layerInfoIdx(2)-5;
        endCh = layerInfoIdx(2)+5;
        [~,minCSD] = min(csdData{2,1}(startCh:endCh));
        layerInfoIdx(2) = startCh+(minCSD-1);  

        %correct slm channel
        startCh = layerInfoIdx(3)-8;
        endCh = min(layerInfoIdx(3)+8,size(csdData{2,1},2));
        [~,maxCSD] = max(csdData{2,1}(startCh:endCh));
        layerInfoIdx(3) = startCh+(maxCSD-1);  

        %add mol layer
        startCh = layerInfoIdx(3)+5;
        endCh = layerInfoIdx(3)+12;
        if size(csdData{2,1},2)==64 && ~strcmp(mice{mm},'IZ34\Final')==1
            [~,minCSD] = min(csdData{2,1}(startCh:endCh));
            layerInfoIdx(4) = startCh+(minCSD-1);  
        else
            layerInfoIdx(4) = nan;
        end
    else
        layerInfoIdx(4) = nan;
    end
    
    if strcmp(mice{mm},'IZ26\Final')==1      
       layerInfoIdx = [12 21 31 43];
    end
    
    for ii = 1:length(reg)
        for jj = 1:length(target) 
            if isempty(csdData{ii,jj})
                csd{ii,jj}(mm,1:4) = nan;            
            elseif isnan(layerInfoIdx(4))   
               if strcmp(mice{mm},'IZ26\Final')==1 && ii<3
                   if ii == 1
                       csd{2,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                       csd{2,jj}(mm,2) = -csdData{ii,jj}(layerInfoIdx(2));                     
                       csd{2,jj}(mm,4) = nan;
                   elseif ii == 2
                       csd{1,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                       csd{1,jj}(mm,2) = -csdData{ii,jj}(layerInfoIdx(2));                     
                       csd{1,jj}(mm,4) = nan;                       
                   end
               else
                   csd{ii,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                   csd{ii,jj}(mm,2) = -csdData{ii,jj}(layerInfoIdx(2));                     
                   csd{ii,jj}(mm,4) = nan;
               end
            else
                if strcmp(mice{mm},'IZ26\Final')==1 && ii<3
                     if ii == 1
                        csd{2,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                        csd{2,jj}(mm,[2 4]) = -csdData{ii,jj}(layerInfoIdx([2 4]));    
                     elseif ii == 2
                        csd{1,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                        csd{1,jj}(mm,[2 4]) = -csdData{ii,jj}(layerInfoIdx([2 4]));                             
                     end
                else
                    csd{ii,jj}(mm,[1 3]) = csdData{ii,jj}(layerInfoIdx([1 3]));                   
                    csd{ii,jj}(mm,[2 4]) = -csdData{ii,jj}(layerInfoIdx([2 4]));    
                end
            end
        end                   
    end
end

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
   %ca3, mec, both
   csdDiff(:,:,1) = -((csd{2,1}+1000)-(csd{3,1}+1000)); 
   csdDiff(:,:,2) = -((csd{2,1}+1000)-(csd{2,2}+1000)); 
   csdDiff(:,:,3) = -((csd{2,1}+1000)-(csd{3,2}+1000)); 
else
   csdDiff(:,:,1) = -((csd{1,1}+1000)-(csd{1,2}+1000)); 
   csdDiff(:,:,2) = -((csd{2,1}+1000)-(csd{2,2}+1000)); 
   csdDiff(:,:,3) = -((csd{3,1}+1000)-(csd{3,2}+1000)); 
end

figure
set(gcf,'Position',[400 400 970 420])
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
for ii = 1:3
    subplot(1,3,ii)
    for yy = 1:size(csdDiff(:,:,ii),2)
        boxchart(yy*ones(size(csdDiff(:,yy,ii))), csdDiff(:,yy,ii), 'BoxFaceColor',[8/243 133/243 161/243],'MarkerStyle','none','Orientation','horizontal')
        hold on
        scatter(csdDiff(:,yy,ii),yy*ones(size(csdDiff(:,yy,ii)))-0.3,10,'MarkerFaceColor',[8/243 133/243 161/243],'MarkerEdgeColor','k')
    end
    xlim([-100 75])
    ylim([0 5])
    set(gca,'YDir','rev')
    line([0 0],[0 5],'Color',[0.5 0.5 0.5],'LineStyle','--')
    xlabel('(Stim-Baseline)')
    title(reg{ii})
    
    % Okay now figure out the stats
    csdTable = reshape(csdDiff(:,:,ii),[numel(csdDiff(:,:,ii)),1]);
    layer = [ones(size(csdDiff(:,:,ii),1),1)*1;ones(size(csdDiff(:,:,ii),1),1)*2;ones(size(csdDiff(:,:,ii),1),1)*3;ones(size(csdDiff(:,:,ii),1),1)*4];
    
    [statsCSD{ii}] = groupStats(csdTable,[layer],'doPlot',false);
    
    % Add stats against 0        
    csdTable = csdDiff(:,:,ii);
    for ss = 1:size(csdTable,2)
        [~,statsCSD{ii}.signrank.p(ss),~,statsCSD{ii}.signrank.stats{ss}] = ttest(csdTable(:,ss));%signrank
        statsCSD{ii}.signrank.n(ss) = sum(~isnan(csdTable(:,ss)));
    end
end

saveas(gcf,strcat(parentDir,'Compiled\Ripples\CompiledCSD\CSD',tag,'.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\CompiledCSD\CSD',tag,'.fig'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\CompiledCSD\CSD',tag,'.eps'),'epsc');
save(strcat(parentDir,'Compiled\Ripples\CompiledCSD\CSD',tag,'.mat'),'statsCSD');

end