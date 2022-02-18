% Compile data across all sessions

function compileMiceCoherence(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'combineSessions',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
savePlot = p.Results.savePlot;
combineSessions = p.Results.combineSessions;

tag = 'CA1'; 

if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ20\Final','IZ21\Final','IZ31\Final'};%,'IZ30\Final','IZ15\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
         'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Saline','IZ28\Saline',...
         'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};  % To add: IZ23, IZ32, IZ33, IZ34
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ33\Final','IZ29\Final','IZ32\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline'}; % To add: IZ32, IZ33, IZ34
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    %
    reg = {'contramEC','ipsimEC','Both'};
end

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};


for ii = 1:length(reg)
    for jj = 1:length(target)
        for kk = 1:length(zone)      
            layerProfile.cohTheta{ii,jj}{kk} = [];
            layerProfile.cohThetaPhase{ii,jj}{kk} = [];
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
                
                cohTheta = [];
                cohThetaPhase = [];

                % One value per animal
                % Number of sessions for condition - 
                for ll = 1:size(layers{ii,jj},1)
                    if ~isempty(csdcohData.cohTheta{ii,jj}{kk})
                        layerInfo = layers{ii,jj}(ll,2:end);
                        [~,layerInfoIdx] = ismember(layerInfo,channelOrder);

                        %correct slm channel
                        %[~,minCSD] = min(csdcohData.csd{ii,jj}{kk}(1,:,ll));
%                             if abs(minCSD-layerInfoIdx(3))<=7
%                                 layerInfoIdx(3) = minCSD;
%                             end
%                             if layerInfoIdx(3) == length(channelOrder)
%                                 layerInfoIdx(3) = length(channelOrder)-1;
%                             end
                        if length(channelOrder)>32
                            %add mol layer
                            startCh = layerInfoIdx(3)+5;
                            endCh = layerInfoIdx(3)+12;
                            if ~strcmp(mice{m},'IZ34\Final')==1
                                [~,maxCSD] = max(csdcohData.csd{ii,jj}{kk}(1,startCh:endCh,ll));
                                layerInfoIdx(4) = startCh+(maxCSD-1);  
                            else
                                layerInfoIdx(4) = nan;
                            end
                        else
                            layerInfoIdx(4) = nan;
                        end
                        
                        if isnan(layerInfoIdx(4))
                            cohTheta(ll,1:3) = csdcohData.cohTheta{ii,jj}{kk}(1,layerInfoIdx(1:3),ll);
                            cohThetaPhase(ll,1:3) = csdcohData.cohPhaseTheta{ii,jj}{kk}(1,layerInfoIdx(1:3),ll);
                            cohTheta(ll,4) = nan;
                            cohThetaPhase(ll,4) = nan;
                        else
                            cohTheta(ll,:) = csdcohData.cohTheta{ii,jj}{kk}(1,layerInfoIdx,ll);
                            cohThetaPhase(ll,:) = csdcohData.cohPhaseTheta{ii,jj}{kk}(1,layerInfoIdx,ll);
                        end
                    end
                end
                
                if combineSessions
                    if ~isempty(cohTheta)
                        layerProfile.cohTheta{ii,jj}{kk}(m,:) = nanmedian(cohTheta,1);
                        layerProfile.cohThetaPhase{ii,jj}{kk}(m,:) = nanmedian(cohThetaPhase,1);
                    end
                else
                    layerProfile.cohTheta{ii,jj}{kk} = [layerProfile.cohTheta{ii,jj}{kk}; cohTheta];
                    layerProfile.cohThetaPhase{ii,jj}{kk} = [layerProfile.cohThetaPhase{ii,jj}{kk}; cohThetaPhase];
                end
                clear cohTheta cohThetaPhase 
            end
        end
    end
    
    clear layers
end
 
for ii = 1:4
    for jj = 1:2
        for kk = 1:3
            cohTheta{ii,jj}{kk} = [];
            cohThetaPhase{ii,jj}{kk} = [];
        end
    end
end

reg = {'CA3','mEC','Both'};
target = {'STEM', 'RETURN'};

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1 
    for jj = 1:length(target)
        for kk = 1:3  % 4 to 6 is the stim condition
            % Ways to normalize           
            % Baseline
            cohTheta{1,jj}{kk} = layerProfile.cohTheta{2,jj}{kk};
            cohThetaPhase{1,jj}{kk} = layerProfile.cohThetaPhase{2,jj}{kk};
            
            % CA3
            cohTheta{2,jj}{kk} = layerProfile.cohTheta{3,jj}{kk};
            cohThetaPhase{2,jj}{kk} = layerProfile.cohThetaPhase{3,jj}{kk};
            
            % mEC
            cohTheta{3,jj}{kk} = layerProfile.cohTheta{2,jj}{kk+3};
            cohThetaPhase{3,jj}{kk} = layerProfile.cohThetaPhase{2,jj}{kk+3};
            
            %Both
            cohTheta{4,jj}{kk} = layerProfile.cohTheta{3,jj}{kk+3};
            cohThetaPhase{4,jj}{kk} = layerProfile.cohThetaPhase{3,jj}{kk+3};
            
        end             
    end
else   
    for ii = 1:length(reg)
        for jj = 1:length(target)
            for kk = 1:3  % 4 to 6 is the stim condition
                
                if strcmp(tag,'mEC') == 1      
                    if ii == 2
                        cohTheta{1,jj}{kk} = [cohTheta{1,jj}{kk};layerProfile.cohTheta{ii,jj}{kk}];
                        cohThetaPhase{1,jj}{kk} = [cohThetaPhase{1,jj}{kk};layerProfile.cohThetaPhase{ii,jj}{kk}];
                    end
                else
                    cohTheta{1,jj}{kk} = [cohTheta{1,jj}{kk};layerProfile.cohTheta{ii,jj}{kk}];
                    cohThetaPhase{1,jj}{kk} = [cohThetaPhase{1,jj}{kk};layerProfile.cohThetaPhase{ii,jj}{kk}];
                end
                
                cohTheta{ii+1,jj}{kk} = layerProfile.cohTheta{ii,jj}{kk+3};
                cohThetaPhase{ii+1,jj}{kk} = layerProfile.cohThetaPhase{ii,jj}{kk+3};
            end                          
        end        
    end 
end

if strcmp(tag,'mEC') == 1
    
    colMat = [8/243 133/243 161/243; 0.5 0.5 0.5];       
    
    for jj = 1%:length(target)
        figure
        set(gcf,'Position',[100 100 1500 1000])
        set(gcf,'renderer','painters')
        for kk = 1:2    
            manip = [];
            layer = [];

            datacohTheta{kk} = [cohTheta{1,jj}{kk} cohTheta{3,jj}{kk}];
            datacohThetaPhase{kk} = [cohThetaPhase{1,jj}{kk} cohThetaPhase{3,jj}{kk}];
            
            %Rearrange
            datacohTheta{kk} = datacohTheta{kk}(:,[1 5 2 6 3 7 4 8]);
            datacohThetaPhase{kk} = datacohThetaPhase{kk}(:,[1 5 2 6 3 7 4 8]);

            subplot(2,2,2*(kk-1)+1)
            for yy = 1:size(datacohTheta{kk},2)
                boxchart(yy*ones(size(datacohTheta{kk}(:,yy))), datacohTheta{kk}(:,yy), 'BoxFaceColor',colMat(rem(yy,2)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacohTheta{kk}(:,yy),yy*ones(size(datacohTheta{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,2)+1,:),'MarkerEdgeColor','k')
            end         
            ylim([0 7])
            xlim([0.4 1])
            set(gca,'YDir','rev')
            xlabel('Coherence')
            title('Coherence')
            
            subplot(2,2,2*(kk-1)+2)
            for yy = 1:size(datacohThetaPhase{kk},2)
                idxtoUse = ~isnan(datacohThetaPhase{kk}(:,yy));
                meantemp = rad2deg(circ_median(datacohThetaPhase{kk}(idxtoUse,yy)));
                line([meantemp(1) meantemp(1)], [yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.4 yy*ones(size(datacohThetaPhase{kk}(:,yy)))+0.4],'Color',colMat(rem(yy,2)+1,:),'LineWidth',3)
                hold on
                scatter(rad2deg(datacohThetaPhase{kk}(:,yy)),yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,2)+1,:),'MarkerEdgeColor','k')
            end         
            ylim([0 7])
            set(gca,'YDir','rev')
            xlabel('Coherence Phase')                     
        end
%         saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.png'));
%         saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.fig'));
%         saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.eps'),'epsc');
    end 

%     % Okay now figure out the stats
%     for kk = 1:2
% 
%         dataTableCSD = [];
%         dataTablethetaPower = [];
%         dataTablesgPower = [];
%         dataTablemgPower = [];
%         dataTablehfoPower = [];
%         manip = [];
%         layer = [];
% 
%         dataTableCSD = datacsd{kk};
%         dataTablethetaPower = datathetaPower{kk};
%         dataTablesgPower = datasgPower{kk};
%         dataTablemgPower = datamgPower{kk};   
%         dataTablehfoPower = datahfoPower{kk};
% 
%         dataTableCSD = reshape(dataTableCSD,[numel(dataTableCSD),1]);
%         dataTablethetaPower = reshape(dataTablethetaPower,[numel(dataTablethetaPower),1]);
%         dataTablesgPower  = reshape(dataTablesgPower,[numel(dataTablesgPower),1]);
%         dataTablemgPower = reshape(dataTablemgPower,[numel(dataTablemgPower),1]);
%         dataTablehfoPower = reshape(dataTablehfoPower,[numel(dataTablehfoPower),1]);
% %         
% %         layer1 = 1;%{'or'};
% %         layer1 = repmat(layer1,[size(datacsd{1},1),1]);
%         layer2 = 1;%{'pyr'};
%         layer2 = repmat(layer2,[size(datacsd{1},1),1]);
%         layer3 = 2;%{'rad'};
%         layer3 = repmat(layer3,[size(datacsd{1},1),1]);
%         layer4 = 3;%{'slm'};
%         layer4 = repmat(layer4,[size(datacsd{1},1),1]);
%         %layer = [layer1; layer2; layer3; layer4];%;layer1; layer2; layer3; layer4];
%         layer = [layer2; layer3; layer4];
% 
%         [statsCSD{kk}] = groupStats(dataTableCSD,[layer],'doPlot',false,'repeatedMeasures',true);
%         [statsthetaPower{kk}] = groupStats(dataTablethetaPower,[layer],'doPlot',false,'repeatedMeasures',true);
%         [statssgPower{kk}] = groupStats(dataTablesgPower,[layer],'doPlot',false,'repeatedMeasures',true);
%         [statsmgPower{kk}] = groupStats(dataTablemgPower,[layer],'doPlot',false,'repeatedMeasures',true);
%         [statshfoPower{kk}] = groupStats(dataTablehfoPower,[layer],'doPlot',false,'repeatedMeasures',true);
%         
%         dataTableCSD = datacsd{kk};
%         dataTablethetaPower = datathetaPower{kk};
%         dataTablesgPower = datasgPower{kk};
%         dataTablemgPower = datamgPower{kk};   
%         dataTablehfoPower = datahfoPower{kk};      
%         
%         for ss = 1:size(dataTablethetaPower,2)
%             [statsCSD{kk}.signrank.p(1,ss), ~, statsCSD{kk}.signrank.stats{ss}] = signrank(dataTableCSD(:,ss));
%             [statsthetaPower{kk}.signrank.p(1,ss), ~, statsthetaPower{kk}.signrank.stats{ss}] = signrank(dataTablethetaPower(:,ss));
%             [statssgPower{kk}.signrank.p(1,ss), ~, statssgPower{kk}.signrank.stats{ss}] = signrank(dataTablesgPower(:,ss));
%             [statsmgPower{kk}.signrank.p(1,ss), ~, statsmgPower{kk}.signrank.stats{ss}] = signrank(dataTablemgPower(:,ss));
%             [statshfoPower{kk}.signrank.p(1,ss), ~, statshfoPower{kk}.signrank.stats{ss}] = signrank(dataTablehfoPower(:,ss));
%         end
%     end
    
else    
    colMat = [56/243 61/243 150/243;...
        0.5 0.5 0.5;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243];   

    for jj = 1%:length(target)
        figure
        set(gcf,'Position',[100 100 1500 1000])
        set(gcf,'renderer','painters')
        for kk = 1:2    
            manip = [];
            layer = [];

            datacohTheta{kk} = [];
            datacohThetaPhase{kk} = [];
            
            for ii = 1:4
                datacohTheta{kk} = catpad(2,datacohTheta{kk},cohTheta{ii,jj}{kk});
                datacohThetaPhase{kk} = catpad(2,datacohThetaPhase{kk},cohThetaPhase{ii,jj}{kk});   
            end

            datacohTheta{kk} = datacohTheta{kk}(:,[1 9 5 13 2 10 6 14 3 11 7 15 4 12 8 16]);
            datacohThetaPhase{kk} = datacohThetaPhase{kk}(:,[1 9 5 13 2 10 6 14 3 11 7 15 4 12 8 16]);
            datacohThetaPhase{kk}(datacohThetaPhase{kk}<-1.5) = datacohThetaPhase{kk}(datacohThetaPhase{kk}<-1.5)+2*pi;

            subplot(2,2,2*(kk-1)+1)
            for yy = 1:size(datacohTheta{kk},2)
                boxchart(yy*ones(size(datacohTheta{kk}(:,yy))), datacohTheta{kk}(:,yy), 'BoxFaceColor', colMat(rem(yy,4)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
                scatter(datacohTheta{kk}(:,yy),yy*ones(size(datacohTheta{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,4)+1,:),'MarkerEdgeColor','k')
            end
            xlim([0.4 1])
            ylim([0 17])
            set(gca,'YDir','rev')
            xlabel('Coherence')
            
            subplot(2,2,2*(kk-1)+2)
            for yy = 1:size(datacohThetaPhase{kk},2)
                
                boxchart(yy*ones(size(datacohThetaPhase{kk}(:,yy))), rad2deg(datacohThetaPhase{kk}(:,yy)), 'BoxFaceColor', colMat(rem(yy,4)+1,:),'MarkerStyle','none','Orientation','horizontal')
                hold on
  
%                 idxtoUse = ~isnan(datacohThetaPhase{kk}(:,yy));
%                 meantemp = rad2deg(circ_median(datacohThetaPhase{kk}(idxtoUse,yy)));
%                 line([meantemp(1) meantemp(1)], [yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.4 yy*ones(size(datacohThetaPhase{kk}(:,yy)))+0.4],'Color',colMat(rem(yy,4)+1,:),'LineWidth',3)                
%                 hold on
%                 line([meantemp(1)+360 meantemp(1)+360], [yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.4 yy*ones(size(datacohThetaPhase{kk}(:,yy)))+0.4],'Color',colMat(rem(yy,4)+1,:),'LineWidth',3)
                scatter(rad2deg(datacohThetaPhase{kk}(:,yy)),yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,4)+1,:),'MarkerEdgeColor','k')
               % scatter(rad2deg(datacohThetaPhase{kk}(:,yy))+360,yy*ones(size(datacohThetaPhase{kk}(:,yy)))-0.3,10,'MarkerFaceColor',colMat(rem(yy,4)+1,:),'MarkerEdgeColor','k')
            end

            xlim([-100 200])
            ylim([0 17])
            set(gca,'YDir','rev')
            xlabel('Coherence Phase')

        end
        saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.png'));
        saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.fig'));
        saveas(gcf,strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'.eps'),'epsc');
    end 

    % Okay now figure out the stats
    for kk = 1:2

        dataTableCoherence = [];
        dataTableCohPhase = [];
        manip = [];
        layer = [];

        dataTableCoherence = datacohTheta{kk};
        dataTableCohPhase = datacohThetaPhase{kk};

        dataTableCoherence = reshape(dataTableCoherence,[numel(dataTableCoherence),1]);
        dataTableCohPhase = reshape(dataTableCohPhase,[numel(dataTableCohPhase),1]);

        manip(1:size(datacohTheta{1},1)) = 1;
        manip(size(datacohTheta{1},1)+1:2*size(datacohTheta{1},1)) = 2;
        manip(2*size(datacohTheta{1},1)+1:3*size(datacohTheta{1},1)) = 3;
        manip(3*size(datacohTheta{1},1)+1:4*size(datacohTheta{1},1)) = 4;
        manip = repmat(manip',[4,1]);

        layer2 = 1;%{'pyr'};
        layer2 = repmat(layer2,[4*size(datacohTheta{1},1),1]);
        layer3 = 2;%{'rad'};
        layer3 = repmat(layer3,[4*size(datacohTheta{1},1),1]);
        layer4 = 3;%{'slm'};
        layer4 = repmat(layer4,[4*size(datacohTheta{1},1),1]);
        layer5 = 4;%{'slm'};
        layer5 = repmat(layer5,[4*size(datacohTheta{1},1),1]);        
        layer = [layer2; layer3; layer4; layer5];

        [statsCoherence{kk}] = groupStats(dataTableCoherence,[manip layer],'doPlot',false);
        [statsCohPhase{kk}] = groupStats(dataTableCohPhase,[manip layer],'doPlot',false);
        statsCoherence{kk}.posthocdim1 = multcompare(statsCoherence{kk}.anova.stats,'Dimension',1,'CType','tukey-kramer','Display','off');
        statsCoherence{kk}.posthocdim2 = multcompare(statsCoherence{kk}.anova.stats,'Dimension',2,'CType','tukey-kramer','Display','off');
        
        statsCohPhase{kk}.posthocdim1 = multcompare(statsCohPhase{kk}.anova.stats,'Dimension',1,'CType','tukey-kramer','Display','off');
        statsCohPhase{kk}.posthocdim2 = multcompare(statsCohPhase{kk}.anova.stats,'Dimension',2,'CType','tukey-kramer','Display','off');
    end
end
save(strcat(parentDir,'\Compiled\Coherence\Coherence',tag,'Stats.mat'),'statsCoherence','statsCohPhase');
end