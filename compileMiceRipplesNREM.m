function compiledMiceRipples = compileMiceRipplesNREM

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

numAnalog = 2;

cmap = cbrewer('qual','Pastel2',length(mice));

compiledMiceRipples.mice = mice;

for ii = 1:(numAnalog+1)
    compiledMiceRipples.numWAKE{ii} = [];
    compiledMiceRipples.numNREM{ii} = [];
end

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m},'\Summ'));  
    if exist('compiledRipplesNREM.mat','file') && ~force
        disp(['Loading ripples for mouse' mice{m}])
        load('compiledRipplesNREM.mat');
    else 
        disp(['Ripples not computed for mouse, calculating...' mice{m}])
        compiledRipples = compileRipplesNREM('expPath',strcat(parentDir, mice{m})); 
    end
    
    for ii = 1:(numAnalog+1)
        numWAKE =[];
        numWAKE = catpad(2,compiledRipples{1}.rate_pre(compiledRipples{1}.rateAnalog == ii)',...
            compiledRipples{1}.rate_post(compiledRipples{1}.rateAnalog == ii)');
        compiledMiceRipples.numWAKE{ii} = [compiledMiceRipples.numWAKE{ii}; numWAKE];

        numNREM =[];
        numNREM = catpad(2,compiledRipples{2}.rate_pre(compiledRipples{2}.rateAnalog == ii)',...
            compiledRipples{2}.rate_post(compiledRipples{2}.rateAnalog == ii)');
        compiledMiceRipples.numNREM{ii} = [compiledMiceRipples.numNREM{ii}; numNREM];        
    end
end

figure
set(gcf,'Position',[50 50 1200 300])
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
colMat = [85/243 85/243 85/243; 85/243 85/243 85/243; 8/243 133/243 161/243; 122/243 122/243 122/243];%; 122/243 122/243 122/243];
col = [85/243 85/243 85/243; 8/243 133/243 161/243];
colormap(col)

if strcmp(tag,'mEC') == 1 
    
    set(gcf,'Position',[200 400 350 330])
    subplot(1,2,1)
    dataAll = [];
    dataAll{1} = compiledMiceRipples.numWAKE{2}(:,1);    
    dataAll{2} = compiledMiceRipples.numWAKE{2}(:,2);    
    stats.rippleWAKE = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);      
   [ stats.rippleWAKE.signrank.p,~,stats.rippleWAKE.signrank.stats] = signrank(dataAll{1},dataAll{2});
    title(num2str(stats.rippleWAKE.signrank.p));
    ylabel('Number of ripples WAKE')
    ylim([0 0.8])
       
    subplot(1,2,2)
    dataAll = [];
    dataAll{1} = compiledMiceRipples.numNREM{2}(:,1);    
    dataAll{2} = compiledMiceRipples.numNREM{2}(:,2);    

    stats.rippleNREM = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);      
   [ stats.rippleNREM.signrank.p,~,stats.rippleNREM.signrank.stats] = signrank(dataAll{1},dataAll{2});
    title(num2str(stats.rippleNREM.signrank.p));
    ylabel('Number of ripples NREM')
    ylim([0 0.8])

else
    set(gcf,'Position',[85 150 450 800])
    for ii = 1:(numAnalog+1)
        subplot((numAnalog+1),2,2*(ii-1)+1) 
        hold on
        dataAll = [];
        dataAll{1} = compiledMiceRipples.numWAKE{ii}(:,1);    
        dataAll{2} = compiledMiceRipples.numWAKE{ii}(:,2);    
        stats.rippleWAKE{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);      
       [ stats.rippleWAKE{ii}.signrank.p,~,stats.rippleWAKE{ii}.signrank.stats] = signrank(dataAll{1},dataAll{2});
        title(num2str(stats.rippleWAKE{ii}.signrank.p));
        ylabel('Number of ripples WAKE')
        ylim([0 0.8])

        subplot((numAnalog+1),2,2*(ii-1)+2) 
        dataAll = [];
        dataAll{1} = compiledMiceRipples.numNREM{ii}(:,1);    
        dataAll{2} = compiledMiceRipples.numNREM{ii}(:,2);    
        stats.rippleNREM{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);      
       [ stats.rippleNREM{ii}.signrank.p,~,stats.rippleNREM{ii}.signrank.stats] = signrank(dataAll{1},dataAll{2});
        title(num2str(stats.rippleNREM{ii}.signrank.p));
        ylabel('Number of ripples NREM')
        ylim([0 0.8])
 
    end    

end


saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesNREM',tag,'.png'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesNREM',tag,'.fig'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesNREM',tag,'.eps'),'epsc');
save(strcat(parentDir,'Compiled\Ripples\compiledRipplesNREM',tag,'.mat'),'stats')

end