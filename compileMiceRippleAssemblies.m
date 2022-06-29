function compileMiceRippleAssemblies(varargin)

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
        if exist(strcat('Summ\RippleOnlyAssemblies.mat'),'file')
            load(strcat('Summ\RippleOnlyAssemblies.mat'));
        else 
            disp(['Assembly not computed for mouse' mice{m}])
            continue;
        end

        for ii = 1:length(range)
            currManip = RippleAssemblies{range(ii)};
            compiledAssemblies{rangetag(ii)} = [compiledAssemblies{rangetag(ii)}; (currManip(:,2)-currManip(:,1))./(currManip(:,2)+currManip(:,1)) ...
                (currManip(:,2)-currManip(:,3))./(currManip(:,2)+currManip(:,3)) ...
                (currManip(:,2)-currManip(:,4))./(currManip(:,2)+currManip(:,4))];% ...
               % (currManip(:,3)-currManip(:,4))./(currManip(:,3)+currManip(:,4))];
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

assemblyData(assemblyData<-10) = nan;
assemblyData(assemblyData>10) = nan;
threshPos = nanmean(assemblyData(:,1))+1*nanstd(assemblyData(:,1));
threshNeg = nanmean(assemblyData(:,1))-1*nanstd(assemblyData(:,1));
ChangeIdx(:,1) = [sum(assemblyData(:,1)<=threshPos & assemblyData(:,1)>=threshNeg);sum(assemblyData(:,1)>threshPos);sum(assemblyData(:,1)<threshNeg)]./sum(~isnan(assemblyData(:,1)));
ChangeIdx(:,2) = [sum(assemblyData(:,2)<=threshPos & assemblyData(:,2)>=threshNeg);sum(assemblyData(:,2)>threshPos);sum(assemblyData(:,2)<threshNeg)]./sum(~isnan(assemblyData(:,2)));
ChangeIdx(:,3) = [sum(assemblyData(:,3)<=threshPos & assemblyData(:,3)>=threshNeg);sum(assemblyData(:,3)>threshPos);sum(assemblyData(:,3)<threshNeg)]./sum(~isnan(assemblyData(:,3)));

subplot(2,3,1)
b = bar(ChangeIdx','stacked','FaceColor','flat');
b(1).CData = [200/243 200/243 200/243];
b(2).CData = [56/243 61/243 150/243];
b(3).CData = [175/243 54/243 60/243];
%legend({'No change','decrease','increase'})
ylabel('Fraction assembly change')
xticklabels({'Pre','During','Post'})
box off

subplot(2,3,2)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.mEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
%title('ipsi mEC')
ylim([0 1])

% contra mEC manipulations
assemblyData = compiledAssemblies{2};
subplot(2,3,4)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.contramEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('contra mEC')
ylim([0 1])

% ipsi mEC manipulations
assemblyData = compiledAssemblies{3};
subplot(2,3,5)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.ipsimEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('ipsi mEC')
ylim([0 1])

% bilateral mEC manipulations
assemblyData = compiledAssemblies{4};
subplot(2,3,6)
dataID = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3];%;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:3),[size(assemblyData,1)*3,1]);
stats.bilteralmEC = groupStats(abs(assembly),dataID,'inAxis',true,'color',colMat);
title('bilateral mEC')
ylim([0 1])

saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleAssemblies.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleAssemblies.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleAssemblies.fig'));
save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleAssemblies.mat'),'stats');
end
