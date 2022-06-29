function compileMiceRippleRankOrder(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
          'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
          'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ26\Final','IZ24\Final','IZ25\Final'};

ripCutoff = [2 5 10 20 30];

for ii = 1:4
    compiledAssemblies{ii} = [];
    numRipples{ii} = [];
end
for tag = 1:2
    if tag==1
        range = 2;
        rangetag = 1;
    elseif tag==2
        range = [1 2 3];
        rangetag = [2 3 4];
    end

    for m = 1:length(mice{tag})       
        cd(strcat(parentDir, mice{tag}{m},'\Summ'));
        if exist(strcat('Summ\RippleRankOrder_DS.mat'),'file')
            load(strcat('Summ\RippleRankOrder_DS.mat'));
        else 
%         if exist(strcat('Summ\RippleRankOrder.mat'),'file')
%             load(strcat('Summ\RippleRankOrder.mat'));
%         else 
            disp(['Rank Order not computed for mouse' mice{tag}{m}])
            continue;
        end

        for ii = 1:length(range)
            compiledAssemblies{rangetag(ii)} = [compiledAssemblies{rangetag(ii)}; rippleRankOrder(:,:,range(ii))];
            numRipples{rangetag(ii)} = [numRipples{rangetag(ii)}; numRip(:,:,range(ii))];
        end
    end
end

colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243;
    56/243 61/243 150/243];
figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[100 68 740 900])

for rr = 1:length(ripCutoff)
    % First all ipsi mEC manipulations    
    assemblyData = compiledAssemblies{1};
    assemblyRip = numRipples{1};    
    idxPick = assemblyRip(:,2)>=ripCutoff(rr);
    if rr == 1
        idxPick(42) = 0;
    end
    assemblyData = assemblyData(idxPick,:);
    
    subplot(4,length(ripCutoff),rr)
    dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2];%;ones(size(assemblyData,1),1)*3];
    %assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
    assembly = reshape(assemblyData(:,1:2),[size(assemblyData,1)*2,1]);
    stats{rr}.mEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    [stats{rr}.mEC.signrank.p,~,stats{rr}.mEC.signrank.stats] = signrank(assemblyData(:,1),assemblyData(:,2));
    %[~,stats{rr}.mEC.signrank.p,~,stats{rr}.mEC.signrank.stats] = ttest(assemblyData(:,1),assemblyData(:,2));
    title(strcat('ipsi mEC ', num2str(stats{rr}.mEC.signrank.p)))
    %ylim([0 1])

    % contra mEC manipulations
    assemblyData = compiledAssemblies{2};
    assemblyRip = numRipples{2};    
    idxPick = assemblyRip(:,2)>=ripCutoff(rr);
    assemblyData = assemblyData(idxPick,:);
    subplot(4,length(ripCutoff),length(ripCutoff)+rr)
    dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2];%;ones(size(assemblyData,1),1)*3];
    assembly = reshape(assemblyData(:,1:2),[size(assemblyData,1)*2,1]);
    if length(assembly)>6
        stats{rr}.contramEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat,'repeatedMeasures',true,'plotType','BoxLinesSEM');
        %[stats{rr}.contramEC.signrank.p,~,stats{rr}.contramEC.signrank.stats] = signrank(assemblyData(:,1),assemblyData(:,2));
        [~,stats{rr}.contramEC.signrank.p,~,stats{rr}.contramEC.signrank.stats] = ttest(assemblyData(:,1),assemblyData(:,2));
        title(strcat('contra mEC ', num2str(stats{rr}.contramEC.signrank.p)))
    end
    %ylim([0 1])

    % ipsi mEC manipulations
    assemblyData = compiledAssemblies{3};
    assemblyRip = numRipples{3};    
    idxPick = assemblyRip(:,2)>=ripCutoff(rr);
    assemblyData = assemblyData(idxPick,:);
    subplot(4,length(ripCutoff),2*length(ripCutoff)+rr)
    dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2];%;ones(size(assemblyData,1),1)*3];
    assembly = reshape(assemblyData(:,1:2),[size(assemblyData,1)*2,1]);
    if length(assembly)>6
        stats{rr}.ipsimEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat,'repeatedMeasures',true,'plotType','BoxLinesSEM');
        %[stats{rr}.ipsimEC.signrank.p,~,stats{rr}.ipsimEC.signrank.stats] = signrank(assemblyData(:,1),assemblyData(:,2));
        [~,stats{rr}.ipsimEC.signrank.p,~,stats{rr}.ipsimEC.signrank.stats] = ttest(assemblyData(:,1),assemblyData(:,2));
        title(strcat('ipsi mEC ', num2str(stats{rr}.ipsimEC.signrank.p)))
    end
    %ylim([0 1])

    % bilateral mEC manipulations
    assemblyData = compiledAssemblies{4};
    assemblyRip = numRipples{4};    
    idxPick = assemblyRip(:,2)>=ripCutoff(rr);
    assemblyData = assemblyData(idxPick,:);
    subplot(4,length(ripCutoff),3*length(ripCutoff)+rr)
    dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2];%ones(size(assemblyData,1),1)*3];
    assembly = reshape(assemblyData(:,1:2),[size(assemblyData,1)*2,1]);
    if length(assembly)>6    
        stats{rr}.bilteralmEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat,'repeatedMeasures',true,'plotType','BoxLinesSEM');
        %[stats{rr}.bilteralmEC.signrank.p,~,stats{rr}.bilteralmEC.signrank.stats] = signrank(assemblyData(:,1),assemblyData(:,2));
        [~,stats{rr}.bilteralmEC.signrank.p,~,stats{rr}.bilteralmEC.signrank.stats] = ttest(assemblyData(:,1),assemblyData(:,2));
        title(strcat('bilateral mEC ', num2str(stats{rr}.bilteralmEC.signrank.p)))
    end
    %ylim([0 1])
end

saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.fig'));
save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.mat'),'stats');

% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.png'));
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.eps'),'epsc');
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.fig'));
% save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.mat'),'stats');

end
