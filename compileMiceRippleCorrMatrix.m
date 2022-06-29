function compileMiceRippleCorrMatrix(varargin)

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
        if exist(strcat('Summ\RippleCorrMatrix.mat'),'file')
            load(strcat('Summ\RippleCorrMatrix.mat'));
        else 
%         if exist(strcat('Summ\RippleClusters.mat'),'file')
%             load(strcat('Summ\RippleClusters.mat'));
%         else 
            disp(['Corr matrix not computed for mouse' mice{m}])
            continue;
        end

        for ii = 1:length(range)
            compiledAssemblies{rangetag(ii)} = [compiledAssemblies{rangetag(ii)}; corrMatrix{range(ii)}];
        end
    end
end

colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243;
    0.5 0.5 0.5];
    %56/243 61/243 150/243];
figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

% First all ipsi mEC manipulations
assemblyData = compiledAssemblies{1};

subplot(1,4,1)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.mEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('ipsi mEC')
%ylim([0 1])

% contra mEC manipulations
assemblyData = compiledAssemblies{2};
subplot(1,4,2)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.contramEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('contra mEC')
%ylim([0 1])

% ipsi mEC manipulations
assemblyData = compiledAssemblies{3};
subplot(1,4,3)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.ipsimEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('ipsi mEC')
%ylim([0 1])

% bilateral mEC manipulations
assemblyData = compiledAssemblies{4};
subplot(1,4,4)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.bilteralmEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('bilateral mEC')
%ylim([0 1])
% 
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCorrMatrix.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCorrMatrix.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCorrMatrix.fig'));
save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCorrMatrix.mat'),'stats');

% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCluster.png'));
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCluster.eps'),'epsc');
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCluster.fig'));
% save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleCluster.mat'),'stats');
end
