function compileMiceRippleRankTemplate(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
 
mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
          'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
          'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ25\Final','IZ26\Final','IZ24\Final'};  

for ii = 1:4
    compiledAssemblies.type1{ii} = [];
    compiledAssemblies.baseTemplate{ii} = [];
    compiledAssemblies.stimTemplate{ii} = [];
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
        if exist(strcat('Summ\templateRankOrder.mat'),'file')
            load(strcat('Summ\templateRankOrder.mat'));
        else 
            disp(['Template rank order not computed for mouse' mice{tag}{m}])
            continue;
        end

        for ii = 1:length(range)
            compiledAssemblies.type1{rangetag(ii)} = [compiledAssemblies.type1{rangetag(ii)}; templateRankOrder.analysis1{range(ii)}];
            compiledAssemblies.baseTemplate{rangetag(ii)} = [compiledAssemblies.baseTemplate{rangetag(ii)}; templateRankOrder.analysis2{range(ii),1}];
            compiledAssemblies.stimTemplate{rangetag(ii)} = [compiledAssemblies.stimTemplate{rangetag(ii)}; templateRankOrder.analysis2{range(ii),2}];
        end
    end
end

col = [0.5 0.5 0.5;
    103/243 189/243 170/243;
    8/243 133/243 161/243;
    56/243 61/243 150/243];

colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243];

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
tag = {'mEC','contramEC','ipsimEC','bimEC'};

% First all type1 analysis  - baseline and stim templates - replayed in the
% homecage
for ii = 1:4
    subplot(3,4,ii)
    assembly = compiledAssemblies.type1{ii};
    data = reshape(assembly,[size(assembly,1)*2,1]);
    dataID = [ones(size(assembly,1),1)*1;ones(size(assembly,1),1)*2];   
    stats.type1.(tag{ii}) = groupStats(data,dataID,'inAxis',true,'color',colMat,'plotType','BoxLinesSEM','repeatedMeasures',true);
    title(tag{ii})
end

% Next type 2 analysis with a baseline template
assembly = [compiledAssemblies.baseTemplate{1};compiledAssemblies.baseTemplate{2};...
    compiledAssemblies.baseTemplate{3};compiledAssemblies.baseTemplate{4}];
data = reshape(assembly,[size(assembly,1)*6,1]);
dataID = [ones(size(assembly,1),1)*1;ones(size(assembly,1),1)*1;ones(size(assembly,1),1)*1;ones(size(assembly,1),1)*2;...   
    ones(size(assembly,1),1)*3;ones(size(assembly,1),1)*4];
subplot(3,4,[5 6])
stats.baseTemplate = groupStats(data,dataID,'inAxis',true,'color',col);%,'repeatedMeasures',true);

col = [0.5 0.5 0.5;
    103/243 189/243 170/243;
    0.5 0.5 0.5;
    8/243 133/243 161/243;
    0.5 0.5 0.5;
    56/243 61/243 150/243];
% Next type 2 analysis with a stim template
for ii = 1:4
    subplot(3,4,8+ii)
    assembly = compiledAssemblies.stimTemplate{ii};
    data = [assembly(:,1);assembly(:,4);assembly(:,2);assembly(:,5);assembly(:,3);assembly(:,6)];
    dataID = [ones(size(assembly,1),1)*1;ones(size(assembly,1),1)*2;ones(size(assembly,1),1)*3;ones(size(assembly,1),1)*4;...   
        ones(size(assembly,1),1)*5;ones(size(assembly,1),1)*6];
    stats.stimTemplate.(tag{ii}) = groupStats(data,dataID,'inAxis',true,'color',col);%,'repeatedMeasures',true);
    title(tag{ii})
end

% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.png'));
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.eps'),'epsc');
% saveas(gcf,strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.fig'));
% save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleRankOrder.mat'),'stats');

end
