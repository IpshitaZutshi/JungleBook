% Compile data across all sessions

function compileMiceHippLayersCSD(varargin)

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
    mice = {'IZ18\Final','IZ20\Final','IZ21\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};% Excluded - 'IZ15\Final','IZ30\Final',
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ31\Final','IZ33\Final'};%'IZ30\Final',,'IZ32\Final','IZ29\Final','IZ15\Final','IZ34\Final',
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

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};


for ii = 1:length(reg)
    for jj = 1:length(target)
        for kk = 1:length(zone)      
            layerProfile.csdPyr{ii,jj}{kk} = [];
            layerProfile.csdOr{ii,jj}{kk} = [];
            layerProfile.csdSlm{ii,jj}{kk} = [];
        end
    end
end


    
%% Loop through the mice
for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    
    if exist(strcat('Summ\csdcohData.mat'),'file')
        load(strcat('Summ\csdcohData.mat'));
    else 
        disp(['CSD Coherence not computed for mouse' mice{m}])
        continue;
    end
    
    
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
        
        load([sessionInfo.FileName '.hippocampalLayers.channelinfo.mat']);
        channelOrder = hippocampalLayers.channelOrder;
        
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
            
            layers{region,target} = [layers{region,target}; hippocampalLayers.all'];
        end
    end
        
    for ii = 1:3
        for jj = 1:2
            for kk = 1:6 
                csdPyr = [];
                csdOr = [];
                csdSlm = [];

                % One value per animal
                % Number of sessions for condition - 
                for ll = 1:size(layers{ii,jj},1)
                    if ~isempty(csdData.orCh{ii,jj}{kk})
                        
                        layerInfo = layers{ii,jj}(ll,2:4);

                        [~,layerInfoIdx] = ismember(layerInfo,channelOrder);
                        
                        if length(channelOrder)>16
                            %correct slm channel
                            startCh = layerInfoIdx(3)-8;
                            endCh = min(layerInfoIdx(3)+8,size(csdData.orCh{ii,jj}{kk},2));
                            [~,minCSD] = min(csdData.orCh{ii,jj}{kk}(1,startCh:endCh,ll));
                            layerInfoIdx(3) = startCh+(minCSD-1);  

                            startCh = layerInfoIdx(2)-2;
                            endCh = min(layerInfoIdx(2)+2,layerInfoIdx(3)-3);
                            [~,maxCSD] = max(csdData.orCh{ii,jj}{kk}(1,startCh:endCh,ll));
                            layerInfoIdx(2) = startCh+(maxCSD-1);  

                            %add mol layer
                            startCh = layerInfoIdx(3)+5;
                            endCh = layerInfoIdx(3)+12;
                            if size(csdData.orCh{ii,jj}{kk},2)==64 && ~strcmp(mice{m},'IZ34\Final')==1
                                [~,maxCSD] = max(csdData.orCh{ii,jj}{kk}(1,startCh:endCh,ll));
                                layerInfoIdx(4) = startCh+(maxCSD-1);  
                            else
                                layerInfoIdx(4) = nan;
                            end
                        else
                            layerInfoIdx(4) = nan;
                        end
            
                        if isnan(layerInfoIdx(4))
                            csdPyr(ll,1:2) = csdcohData.csd{ii,jj}{kk}(1,layerInfoIdx(1:2),ll);
                            csdOr(ll,1:2) = csdData.orCh{ii,jj}{kk}(1,layerInfoIdx(1:2),ll);
                            csdSlm(ll,1:2) = csdData.slmCh{ii,jj}{kk}(1,layerInfoIdx(1:2),ll);
                            csdPyr(ll,3) = -csdcohData.csd{ii,jj}{kk}(1,layerInfoIdx(3),ll);
                            csdOr(ll,3) = -csdData.orCh{ii,jj}{kk}(1,layerInfoIdx(3),ll);
                            csdSlm(ll,3) = -csdData.slmCh{ii,jj}{kk}(1,layerInfoIdx(3),ll);                            
                            csdPyr(ll,4) = nan;
                            csdOr(ll,4) = nan;
                            csdSlm(ll,4) = nan;
                        else
                            csdPyr(ll,:) = csdcohData.csd{ii,jj}{kk}(1,layerInfoIdx,ll);
                            csdOr(ll,:) = csdData.orCh{ii,jj}{kk}(1,layerInfoIdx,ll);
                            csdSlm(ll,:) = csdData.slmCh{ii,jj}{kk}(1,layerInfoIdx,ll);
                            csdPyr(ll,3) = -csdPyr(ll,3);
                            csdOr(ll,3) = -csdOr(ll,3);
                            csdSlm(ll,3) = -csdSlm(ll,3);
                        end
                    end
                end
                
                if combineSessions
                    if ~isempty(csdPyr)
                        layerProfile.csdPyr{ii,jj}{kk}(m,:) = nanmedian(csdPyr,1);
                        layerProfile.csdOr{ii,jj}{kk}(m,:) = nanmedian(csdOr,1);
                        layerProfile.csdSlm{ii,jj}{kk}(m,:) = nanmedian(csdSlm,1);                        
                    end
                else
                    layerProfile.csdPyr{ii,jj}{kk} = [layerProfile.csdPyr{ii,jj}{kk}; csdPyr];
                    layerProfile.csdOr{ii,jj}{kk} = [layerProfile.csdOr{ii,jj}{kk}; csdOr];
                    layerProfile.csdSlm{ii,jj}{kk} = [layerProfile.csdSlm{ii,jj}{kk}; csdSlm];
                end
                clear csdPyr csdOr csdSlm
            end
        end
    end
    
    clear layers
end
 
for ii = 1:length(reg)
    for jj = 1:length(target)
        csdPyr{ii,jj} = [];
        csdOr{ii,jj} = [];
        csdSlm{ii,jj} = [];
    end
end

reg = {'CA3','mEC','Both'};
target = {'STEM', 'RETURN'};

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1 
    for jj = 1:length(target)
        for kk = 1:3  % 4 to 6 is the stim condition
            % Ways to normalize           
            % CA3
            csdPyr{1,jj}{kk} = ((layerProfile.csdPyr{3,jj}{kk}+300)-(layerProfile.csdPyr{2,jj}{kk}+300));%./(abs(layerProfile.csdPyr{2,jj}{kk}));%+abs(layerProfile.csdPyr{3,jj}{kk}));
            csdOr{1,jj}{kk} = ((layerProfile.csdOr{3,jj}{kk}+300)-(layerProfile.csdOr{2,jj}{kk}+300));%./(abs(layerProfile.csdOr{2,jj}{kk}));%+abs(layerProfile.csdOr{3,jj}{kk}));
            csdSlm{1,jj}{kk} = -((layerProfile.csdSlm{3,jj}{kk}+300)-(layerProfile.csdSlm{2,jj}{kk}+300));%./(abs(layerProfile.csdSlm{2,jj}{kk}));%+abs(layerProfile.csdSlm{3,jj}{kk}));
            % mEC
            csdPyr{2,jj}{kk} = ((layerProfile.csdPyr{2,jj}{kk+3}+300)-(layerProfile.csdPyr{2,jj}{kk}+300));%./(abs(layerProfile.csdPyr{2,jj}{kk}));%+abs(layerProfile.csdPyr{2,jj}{kk+3}));
            csdOr{2,jj}{kk} = ((layerProfile.csdOr{2,jj}{kk+3}+300)-(layerProfile.csdOr{2,jj}{kk}+300));%./(abs(layerProfile.csdOr{2,jj}{kk}));%+abs(layerProfile.csdOr{2,jj}{kk+3}));
            csdSlm{2,jj}{kk} = -((layerProfile.csdSlm{2,jj}{kk+3}+300)-(layerProfile.csdSlm{2,jj}{kk}+300));%./(abs(layerProfile.csdSlm{2,jj}{kk}));%+abs(layerProfile.csdSlm{2,jj}{kk+3}));
            
            %Both
            csdPyr{3,jj}{kk} = ((layerProfile.csdPyr{3,jj}{kk+3}+300)-(layerProfile.csdPyr{2,jj}{kk}+300));%./(abs(layerProfile.csdPyr{2,jj}{kk}));%+abs(layerProfile.csdPyr{3,jj}{kk+3}));
            csdOr{3,jj}{kk} = ((layerProfile.csdOr{3,jj}{kk+3}+300)-(layerProfile.csdOr{2,jj}{kk}+300));%./(abs(layerProfile.csdOr{2,jj}{kk}));%+abs(layerProfile.csdOr{3,jj}{kk+3}));
            csdSlm{3,jj}{kk} = -((layerProfile.csdSlm{3,jj}{kk+3}+300)-(layerProfile.csdSlm{2,jj}{kk}+300));%./(abs(layerProfile.csdSlm{2,jj}{kk}));%+abs(layerProfile.csdSlm{3,jj}{kk+3}));
% 
%             csdPyr{1,jj}{kk} = ((layerProfile.csdPyr{3,jj}{kk}+300)-(layerProfile.csdPyr{2,jj}{kk}+300))./(abs(layerProfile.csdPyr{2,jj}{kk})+abs(layerProfile.csdPyr{3,jj}{kk}));
%             csdOr{1,jj}{kk} = ((layerProfile.csdOr{3,jj}{kk}+300)-(layerProfile.csdOr{2,jj}{kk}+300))./(abs(layerProfile.csdOr{2,jj}{kk})+abs(layerProfile.csdOr{3,jj}{kk}));
%             csdSlm{1,jj}{kk} = -((layerProfile.csdSlm{3,jj}{kk}+300)-(layerProfile.csdSlm{2,jj}{kk}+300))./(abs(layerProfile.csdSlm{2,jj}{kk})+abs(layerProfile.csdSlm{3,jj}{kk}));
%             % mEC
%             csdPyr{2,jj}{kk} = ((layerProfile.csdPyr{2,jj}{kk+3}+300)-(layerProfile.csdPyr{2,jj}{kk}+300))./(abs(layerProfile.csdPyr{2,jj}{kk})+abs(layerProfile.csdPyr{2,jj}{kk+3}));
%             csdOr{2,jj}{kk} = ((layerProfile.csdOr{2,jj}{kk+3}+300)-(layerProfile.csdOr{2,jj}{kk}+300))./(abs(layerProfile.csdOr{2,jj}{kk})+abs(layerProfile.csdOr{2,jj}{kk+3}));
%             csdSlm{2,jj}{kk} = -((layerProfile.csdSlm{2,jj}{kk+3}+300)-(layerProfile.csdSlm{2,jj}{kk}+300))./(abs(layerProfile.csdSlm{2,jj}{kk})+abs(layerProfile.csdSlm{2,jj}{kk+3}));
%             
%             %Both
%             csdPyr{3,jj}{kk} = ((layerProfile.csdPyr{3,jj}{kk+3}+300)-(layerProfile.csdPyr{2,jj}{kk}+300))./(abs(layerProfile.csdPyr{2,jj}{kk})+abs(layerProfile.csdPyr{3,jj}{kk+3}));
%             csdOr{3,jj}{kk} = ((layerProfile.csdOr{3,jj}{kk+3}+300)-(layerProfile.csdOr{2,jj}{kk}+300))./(abs(layerProfile.csdOr{2,jj}{kk})+abs(layerProfile.csdOr{3,jj}{kk+3}));
%             csdSlm{3,jj}{kk} = -((layerProfile.csdSlm{3,jj}{kk+3}+300)-(layerProfile.csdSlm{2,jj}{kk}+300))./(abs(layerProfile.csdSlm{2,jj}{kk})+abs(layerProfile.csdSlm{3,jj}{kk+3}));            
        end             
    end
else
    for ii = 1:length(reg)
        for jj = 1:length(target)
            for kk = 1:3  % 4 to 6 is the stim condition
               % Ways to normalize
                csdPyr{ii,jj}{kk} = ((layerProfile.csdPyr{ii,jj}{kk+3}+300)-(layerProfile.csdPyr{ii,jj}{kk}+300));%./(abs(layerProfile.csdPyr{ii,jj}{kk}));%+abs(layerProfile.csdPyr{ii,jj}{kk+3}));
                csdOr{ii,jj}{kk} = ((layerProfile.csdOr{ii,jj}{kk+3}+300)-(layerProfile.csdOr{ii,jj}{kk}+300));%./(abs(layerProfile.csdOr{ii,jj}{kk}));%+abs(layerProfile.csdOr{ii,jj}{kk+3}));
                csdSlm{ii,jj}{kk} = -((layerProfile.csdSlm{ii,jj}{kk+3}+300)-(layerProfile.csdSlm{ii,jj}{kk}+300));%./(abs(layerProfile.csdSlm{ii,jj}{kk}));%+abs(layerProfile.csdSlm{ii,jj}{kk+3}));               
                
%                 csdPyr{ii,jj}{kk} = ((layerProfile.csdPyr{ii,jj}{kk+3}+300)-(layerProfile.csdPyr{ii,jj}{kk}+300))./(abs(layerProfile.csdPyr{ii,jj}{kk})+abs(layerProfile.csdPyr{ii,jj}{kk+3}));
%                 csdOr{ii,jj}{kk} = ((layerProfile.csdOr{ii,jj}{kk+3}+300)-(layerProfile.csdOr{ii,jj}{kk}+300))./(abs(layerProfile.csdOr{ii,jj}{kk})+abs(layerProfile.csdOr{ii,jj}{kk+3}));
%                 csdSlm{ii,jj}{kk} = ((layerProfile.csdOr{ii,jj}{kk+3}+300)-(layerProfile.csdOr{ii,jj}{kk}+300))./(abs(layerProfile.csdSlm{ii,jj}{kk})+abs(layerProfile.csdSlm{ii,jj}{kk+3}));
            end                    
        end        
    end 
end

if strcmp(tag,'mEC') == 1
    for jj = 1%:length(target)
        figure
        set(gcf,'Position',[100 100 1500 1000])
        set(gcf,'renderer','painters')
        for kk = 1:2    
            manip = [];
            layer = [];

            datacsdPyr{kk} = [];
            datacsdOr{kk} = [];
            datacsdSlm{kk} = [];
            
            datacsdPyr{kk} = catpad(2,datacsdPyr{kk},csdPyr{2,jj}{kk});
            datacsdOr{kk} = catpad(2,datacsdOr{kk},csdOr{2,jj}{kk});
            datacsdSlm{kk} = catpad(2,datacsdSlm{kk},csdSlm{2,jj}{kk});


            subplot(2,3,3*(kk-1)+1)
            for yy = 1:size(datacsdPyr{kk},2)
                boxchart(yy*ones(size(datacsdPyr{kk}(:,yy))), datacsdPyr{kk}(:,yy), 'BoxFaceColor',[8/243 133/243 161/243],'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdPyr{kk}(:,yy),yy*ones(size(datacsdPyr{kk}(:,yy)))-0.3,10,'MarkerFaceColor',[8/243 133/243 161/243],'MarkerEdgeColor','k')
            end
            xlim([-100 50])
            ylim([0 5])
            set(gca,'YDir','rev')
            line([0 0],[0 5],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD Pyr')

            subplot(2,3,3*(kk-1)+2)
            for yy = 1:size(datacsdOr{kk},2)
                boxchart(yy*ones(size(datacsdOr{kk}(:,yy))), datacsdOr{kk}(:,yy), 'BoxFaceColor',[8/243 133/243 161/243],'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdOr{kk}(:,yy),yy*ones(size(datacsdOr{kk}(:,yy)))-0.3,10,'MarkerFaceColor',[8/243 133/243 161/243],'MarkerEdgeColor','k')
            end
            ylim([0 5])
            xlim([-100 50])
            set(gca,'YDir','rev')
            line([0 0],[0 5],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD Or')

            subplot(2,3,3*(kk-1)+3)
            for yy = 1:size(datacsdSlm{kk},2)
                boxchart(yy*ones(size(datacsdSlm{kk}(:,yy))), datacsdSlm{kk}(:,yy), 'BoxFaceColor',[8/243 133/243 161/243],'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdSlm{kk}(:,yy),yy*ones(size(datacsdSlm{kk}(:,yy)))-0.3,10,'MarkerFaceColor',[8/243 133/243 161/243],'MarkerEdgeColor','k')
            end
            ylim([0 5])
            xlim([-100 50])
            set(gca,'YDir','rev')
            line([0 0],[0 5],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD slm')
                       
        end
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.png'));
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.fig'));
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.eps'),'epsc');
    end 

    % Okay now figure out the stats
    for kk = 1:2

        dataTableCSDPyr = [];
        dataTableCSDOr = [];
        dataTableCSDSlm = [];

        manip = [];
        layer = [];

        dataTableCSDPyr = datacsdPyr{kk};
        dataTableCSDOr = datacsdOr{kk};
        dataTableCSDSlm = datacsdSlm{kk};

        dataTableCSDPyr = reshape(dataTableCSDPyr,[numel(dataTableCSDPyr),1]);
        dataTableCSDOr = reshape(dataTableCSDOr,[numel(dataTableCSDOr),1]);
        dataTableCSDSlm = reshape(dataTableCSDSlm,[numel(dataTableCSDSlm),1]);
        
        layer = [ones(size(datacsdPyr{1},1),1)*1;ones(size(datacsdPyr{1},1),1)*2;ones(size(datacsdPyr{1},1),1)*3;ones(size(datacsdPyr{1},1),1)*4];

        [statsCSDPyr{kk}] = groupStats(dataTableCSDPyr,[layer],'doPlot',false);
        [statsCSDOr{kk}] = groupStats(dataTableCSDOr,[layer],'doPlot',false);
        [statsCSDSlm{kk}] = groupStats(dataTableCSDSlm,[layer],'doPlot',false);
        
        % Add stats against 0        
        dataTableCSDPyr = datacsdPyr{kk};
        dataTableCSDOr = datacsdOr{kk};
        dataTableCSDSlm = datacsdSlm{kk};
        for ss = 1:size(dataTableCSDPyr,2)
            [statsCSDPyr{kk}.signrank.p(ss),~,statsCSDPyr{kk}.signrank.stats{ss}] = signrank(dataTableCSDPyr(:,ss));
            statsCSDPyr{kk}.signrank.n(ss) = size(dataTableCSDPyr,1);
            [statsCSDOr{kk}.signrank.p(ss),~,statsCSDOr{kk}.signrank.stats{ss}] = signrank(dataTableCSDOr(:,ss));
            statsCSDOr{kk}.signrank.n(ss) = size(dataTableCSDOr,1);
            [statsCSDSlm{kk}.signrank.p(ss),~,statsCSDSlm{kk}.signrank.stats{ss}] = signrank(dataTableCSDSlm(:,ss));
            statsCSDSlm{kk}.signrank.n(ss) = size(dataTableCSDSlm,1);
        end
    end
    
else    
    colMat = [56/243 61/243 150/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243];   

    for jj = 1%:length(target)
        figure
        set(gcf,'Position',[100 100 1500 1000])
        set(gcf,'renderer','painters')
        for kk = 1:2    
            manip = [];
            layer = [];

            datacsdPyr{kk} = [];
            datacsdOr{kk} = [];
            datacsdSlm{kk} = [];

            for ii = 1:length(reg)
                
                datacsdPyr{kk} = catpad(2,datacsdPyr{kk},csdPyr{ii,jj}{kk});
                datacsdOr{kk} = catpad(2,datacsdOr{kk},csdOr{ii,jj}{kk});
                datacsdSlm{kk} = catpad(2,datacsdSlm{kk},csdSlm{ii,jj}{kk});
    
            end

            datacsdPyr{kk} = datacsdPyr{kk}(:,[5 1 9 6 2 10 7 3 11 8 4 12]);
            datacsdOr{kk} = datacsdOr{kk}(:,[5 1 9 6 2 10 7 3 11 8 4 12]);
            datacsdSlm{kk} = datacsdSlm{kk}(:,[5 1 9 6 2 10 7 3 11 8 4 12]);

            subplot(2,3,3*(kk-1)+1)
            for yy = 1:size(datacsdPyr{kk},2)
                boxchart(yy*ones(size(datacsdPyr{kk}(:,yy))), datacsdPyr{kk}(:,yy), 'BoxFaceColor',colMat(rem(yy,3)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdPyr{kk}(:,yy),yy*ones(size(datacsdPyr{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,3)+1,:),'MarkerEdgeColor','k')
            end
            ylim([0 13])
            xlim([-100 20])
            set(gca,'YDir','rev')
            line([0 0],[0 13],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD Pyr')

            subplot(2,3,3*(kk-1)+2)
            for yy = 1:size(datacsdOr{kk},2)
                boxchart(yy*ones(size(datacsdOr{kk}(:,yy))), datacsdOr{kk}(:,yy), 'BoxFaceColor',colMat(rem(yy,3)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdOr{kk}(:,yy),yy*ones(size(datacsdOr{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,3)+1,:),'MarkerEdgeColor','k')
            end
            ylim([0 13])
            xlim([-100 20])
            set(gca,'YDir','rev')
            line([0 0],[0 13],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD Or')

            subplot(2,3,3*(kk-1)+3)
            for yy = 1:size(datacsdSlm{kk},2)
                boxchart(yy*ones(size(datacsdSlm{kk}(:,yy))), datacsdSlm{kk}(:,yy), 'BoxFaceColor',colMat(rem(yy,3)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacsdSlm{kk}(:,yy),yy*ones(size(datacsdSlm{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,3)+1,:),'MarkerEdgeColor','k')
            end
            ylim([0 13])
            xlim([-100 20])
            set(gca,'YDir','rev')
            line([0 0],[0 13],'Color',[0.5 0.5 0.5],'LineStyle','--')
            xlabel('(Stim-Baseline)')
            title('CSD slm')
%                         
        end
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.png'));
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.fig'));
        saveas(gcf,strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.eps'),'epsc');
    end 

    % Okay now figure out the stats
    for kk = 1:2

        dataTableCSDPyr = [];
        dataTableCSDOr = [];
        dataTableCSDSlm = [];

        manip = [];
        layer = [];

        dataTableCSDPyr = datacsdPyr{kk};
        dataTableCSDOr = datacsdOr{kk};
        dataTableCSDSlm = datacsdSlm{kk};

        dataTableCSDPyr = reshape(dataTableCSDPyr,[numel(dataTableCSDPyr),1]);
        dataTableCSDOr = reshape(dataTableCSDOr,[numel(dataTableCSDOr),1]);
        dataTableCSDSlm = reshape(dataTableCSDSlm,[numel(dataTableCSDSlm),1]);
        
        if size(datacsdOr{kk},2) == 12
            layer = [ones(3*size(datacsdPyr{1},1),1)*4;ones(3*size(datacsdPyr{1},1),1)*5;...
                ones(3*size(datacsdPyr{1},1),1)*6; ones(3*size(datacsdPyr{1},1),1)*7];
        else
            layer = [ones(3*size(datacsdPyr{1},1),1)*4;ones(3*size(datacsdPyr{1},1),1)*5;ones(3*size(datacsdPyr{1},1),1)*6];
        end
        manip(1:size(datacsdPyr{1},1)) = 1;
        manip(size(datacsdPyr{1},1)+1:2*size(datacsdPyr{1},1)) = 2;
        manip(2*size(datacsdPyr{1},1)+1:3*size(datacsdPyr{1},1)) = 3;
        if size(datacsdOr{kk},2) == 12
            manip = repmat(manip',[4,1]);
        else
            manip = repmat(manip',[3,1]);
        end


        [statsCSDPyr{kk}] = groupStats(dataTableCSDPyr,[manip layer],'doPlot',false);
        [statsCSDOr{kk}] = groupStats(dataTableCSDOr,[manip layer],'doPlot',false);
        [statsCSDSlm{kk}] = groupStats(dataTableCSDSlm,[manip layer],'doPlot',false);
        
        % Add stats against 0
        
        dataTableCSDPyr = datacsdPyr{kk};
        dataTableCSDOr = datacsdOr{kk};
        dataTableCSDSlm = datacsdSlm{kk};
        for ss = 1:size(dataTableCSDPyr,2)
            [statsCSDPyr{kk}.signrank.p(ss),~,statsCSDPyr{kk}.signrank.stats{ss}] = signrank(dataTableCSDPyr(:,ss));
            statsCSDPyr{kk}.signrank.n(ss) = size(size(dataTableCSDPyr,2),1);
            [statsCSDOr{kk}.signrank.p(ss),~,statsCSDOr{kk}.signrank.stats{ss}] = signrank(dataTableCSDOr(:,ss));
            statsCSDOr{kk}.signrank.n(ss) = size(size(dataTableCSDOr,2),1);
            [statsCSDSlm{kk}.signrank.p(ss),~,statsCSDSlm{kk}.signrank.stats{ss}] = signrank(dataTableCSDSlm(:,ss));
            statsCSDSlm{kk}.signrank.n(ss) = size(size(dataTableCSDSlm,2),1);
        end
    end
end
 save(strcat(parentDir,'\Compiled\Power Profile\CompiledProfileLayers\CSD',tag,'.mat'),'statsCSDPyr','statsCSDOr',...
     'statsCSDSlm');
end