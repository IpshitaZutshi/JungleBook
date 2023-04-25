function compiledMiceRipples = compileMiceRipplesThreshVary

tag = 'mEC';% mEC, CA3, Bilateral mEC
force = 0;
Control = 0;

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','CA1Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ21\Final',...
         'IZ18\Final','IZ20\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
         'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline'}; % To add, IZ16, IZ23'IZ11\Final','IZ34\Saline'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final','IZ27\Final','IZ28\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ33\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Bilateral'};
end


parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

compiledMiceRipples.mice = [];
compiledMiceRipples.rate = [];
dur{1} = [];dur{2} = [];

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m},'\Summ'));  
    if exist('compiledRipplesThresh.mat','file') && ~force
        disp(['Loading ripples for mouse' mice{m}])
        load('compiledRipplesThresh.mat');
    end
    tempRip = [];
    
    for ii = 1:8
        tempRip(:,ii) = compiledRipples.rate{ii};
    end

    numSess = size(tempRip,1);
    
    compiledMiceRipples.rate = [compiledMiceRipples.rate; tempRip];
    compiledMiceRipples.mice = [compiledMiceRipples.mice; ones(numSess,1)*m];
    dur{1} = [dur{1}; compiledRipples.dur{1}];
    dur{2} = [dur{2}; compiledRipples.dur{2}];
end

figure
set(gcf,'Position',[50 50 1200 850])
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
col = [85/243 85/243 85/243; 8/243 133/243 161/243];

subplot(2,3,1)
data{1} = compiledMiceRipples.rate(:,1);
data{2} = compiledMiceRipples.rate(:,2);
stats.rate{1} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','color',col,'inAxis',true);
%ylim([0 0.3])

subplot(2,3,2)
stats.dur = groupStats(dur,[], 'plotType','violinPlot','color',col,'inAxis',true);

subplot(2,3,3)
data{1} = compiledMiceRipples.rate(:,1);
data{2} = compiledMiceRipples.rate(:,2);
rMice  = unique(compiledMiceRipples.mice);
for ii = 1:length(rMice)
    idx = compiledMiceRipples.mice == rMice(ii);
    dat = compiledMiceRipples.rate(idx,:);
    datAv(ii,:) = nanmean(dat,1);
end
data{1} = datAv(:,1);
data{2} = datAv(:,2);
stats.rate{5} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','color',col,'inAxis',true);

subplot(2,3,4)
data{1} = compiledMiceRipples.rate(:,3);
data{2} = compiledMiceRipples.rate(:,4);
stats.rate{2} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','color',col,'inAxis',true);
%ylim([0 0.3])

subplot(2,3,5)
data{1} = compiledMiceRipples.rate(:,5);
data{2} = compiledMiceRipples.rate(:,6);
stats.rate{3} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','color',col,'inAxis',true);
%ylim([0 0.3])

subplot(2,3,6)
data{1} = compiledMiceRipples.rate(:,7);
data{2} = compiledMiceRipples.rate(:,8);
stats.rate{4} = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','color',col,'inAxis',true);
%ylim([0 0.3])
% if ~Control
%     saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.png'));
%     saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.fig'));
%     saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.eps'),'epsc');
%     save(strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.mat'),'stats')
% else
saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\compiledRipplesVaryThresh',tag,'.png'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\compiledRipplesVaryThresh',tag,'.fig'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\Revisions\compiledRipplesVaryThresh',tag,'.eps'),'epsc');
save(strcat(parentDir,'Compiled\Ripples\Revisions\compiledRipplesVaryThresh',tag,'.mat'),'stats')
% end
end