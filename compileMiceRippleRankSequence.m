function compileMiceRippleRankSequence(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
 
mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
          'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
          'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ25\Final','IZ26\Final','IZ24\Final'};  

for ii = 1:4
    compiledAssemblies{ii} = [];
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
        if exist(strcat('Summ\RippleRankSequence_DS.mat'),'file')
            load(strcat('Summ\RippleRankSequence_DS.mat'));
        else 
%         if exist(strcat('Summ\RippleRankSequence.mat'),'file')
%             load(strcat('Summ\RippleRankSequence.mat'));
%         else 
            disp(['Rank Order not computed for mouse' mice{tag}{m}])
            continue;
        end

        for ii = 1:length(range)
            compiledAssemblies{rangetag(ii)} = [compiledAssemblies{rangetag(ii)}; rippleRankSequence{range(ii)}];
        end
    end
end

colMat = [103/243 189/243 170/243;
    8/243 133/243 161/243;
    56/243 61/243 150/243];
figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

% First all ipsi mEC manipulations
assemblyData = compiledAssemblies{1};
assemblyBase = assemblyData(:,1:4);
assemblyStim = assemblyData(:,5:8);
assemblyDiff = (assemblyStim-assemblyBase);
avgCorr = nanmean(assemblyDiff,1);
stderrCorr = nanstd(assemblyDiff,1)./sqrt(size(assemblyDiff,1));
dataID = [ones(size(assemblyDiff,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3; ...
    ones(size(assemblyDiff,1),1)*4];%;ones(size(assemblyData,1),1)*5];
assembly = reshape(assemblyDiff(:,1:4),[size(assemblyData,1)*4,1]);
subplot(1,4,1)
stats.mEC = groupStats(assembly,dataID,'doPlot',false,'repeatedMeasures',true);
for ii = 1:size(assemblyDiff,2)
    [stats.mEC.signrank.p(ii),~,stats.mEC.signrank.stats{ii}] = signrank(assemblyDiff(:,ii));
end
fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],colMat(2,:),'linestyle','none','FaceAlpha',0.2);                   
hold on                 
plot(avgCorr,'color',colMat(2,:),'LineWidth',1.5);
title('ipsi mEC')
ylim([-0.3 0.1])
line([1 4],[0 0],'color','k')
ylabel('Difference in rank order correlation')
xlabel('Ripple quadrant')

% contra mEC manipulations
assemblyData = compiledAssemblies{2};
assemblyBase = assemblyData(:,1:4);
assemblyStim = assemblyData(:,5:8);
avgCorr = nanmean(assemblyStim-assemblyBase,1);
stderrCorr = nanstd(assemblyStim-assemblyBase,1)./sqrt(size(assemblyBase,1));
assemblyDiff = assemblyStim-assemblyBase;
dataID = [ones(size(assemblyDiff,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3; ...
    ones(size(assemblyDiff,1),1)*4];%;ones(size(assemblyData,1),1)*5];
assembly = reshape(assemblyDiff(:,1:4),[size(assemblyData,1)*4,1]);
subplot(1,4,2)
stats.contramEC = groupStats(assembly,dataID,'doPlot',false,'repeatedMeasures',true);
for ii = 1:size(assemblyDiff,2)
    [stats.contramEC.signrank.p(ii),~,stats.contramEC.signrank.stats{ii}] = signrank(assemblyDiff(:,ii));
end
fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],colMat(1,:),'linestyle','none','FaceAlpha',0.2);                   
hold on                 
plot(avgCorr,'color',colMat(1,:),'LineWidth',1.5);
title('contra mEC')
ylim([-0.3 0.1])
line([1 4],[0 0],'color','k')
ylabel('Difference in rank order correlation')
xlabel('Ripple quadrant')

% ipsi mEC manipulations
assemblyData = compiledAssemblies{3};
assemblyBase = assemblyData(:,1:4);
assemblyStim = assemblyData(:,5:8);
assemblyDiff = assemblyStim-assemblyBase;
dataID = [ones(size(assemblyDiff,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3; ...
    ones(size(assemblyDiff,1),1)*4];%;ones(size(assemblyData,1),1)*5];
assembly = reshape(assemblyDiff(:,1:4),[size(assemblyData,1)*4,1]);
subplot(1,4,3)
stats.ipsimEC = groupStats(assembly,dataID,'doPlot',false,'repeatedMeasures',true);
for ii = 1:size(assemblyDiff,2)
    [stats.ipsimEC.signrank.p(ii),~,stats.ipsimEC.signrank.stats{ii}] = signrank(assemblyDiff(:,ii));
end
avgCorr = nanmean(assemblyStim-assemblyBase,1);
stderrCorr = nanstd(assemblyStim-assemblyBase,1)./sqrt(size(assemblyBase,1));
fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],colMat(2,:),'linestyle','none','FaceAlpha',0.2);                   
hold on                 
plot(avgCorr,'color',colMat(2,:),'LineWidth',1.5);
title('ipsi mEC')
ylim([-0.3 0.1])
line([1 4],[0 0],'color','k')
ylabel('Difference in rank order correlation')
xlabel('Ripple quadrant')

% bilateral mEC manipulations
assemblyData = compiledAssemblies{4};
assemblyBase = assemblyData(:,1:4);
assemblyStim = assemblyData(:,5:8);
avgCorr = nanmean(assemblyStim-assemblyBase,1);
stderrCorr = nanstd(assemblyStim-assemblyBase,1)./sqrt(size(assemblyBase,1));
assemblyDiff = assemblyStim-assemblyBase;
dataID = [ones(size(assemblyDiff,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3; ...
    ones(size(assemblyDiff,1),1)*4];%;ones(size(assemblyData,1),1)*5];
assembly = reshape(assemblyDiff(:,1:4),[size(assemblyData,1)*4,1]);
subplot(1,4,4)
stats.bimEC = groupStats(assembly,dataID,'doPlot',false,'repeatedMeasures',true);
for ii = 1:size(assemblyDiff,2)
    [stats.bimEC.signrank.p(ii),~,stats.bimEC.signrank.stats{ii}] = signrank(assemblyDiff(:,ii));
end
fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],colMat(3,:),'linestyle','none','FaceAlpha',0.2);                   
hold on                 
plot(avgCorr,'color',colMat(3,:),'LineWidth',1.5);
title('bilateral mEC')
ylim([-0.3 0.1])
line([1 4],[0 0],'color','k')
ylabel('Difference in rank order correlation')
xlabel('Ripple quadrant')


% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.png'));
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.eps'),'epsc');
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.fig'));
% save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder_DS.mat'),'stats');

% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.png'));
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.eps'),'epsc');
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.fig'));
% save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.mat'),'stats');

end
