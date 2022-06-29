function compiledMiceRipples = compileMiceRipplesmEConly


mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
     'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
     'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add, IZ16, IZ23
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
    compiledMiceRipples.num{ii} = [];
    compiledMiceRipples.amplitude{ii} = [];
    compiledMiceRipples.duration{ii} = [];
    compiledMiceRipples.frequency{ii} = [];
    compiledMiceRipples.swMag{ii} = [];
    compiledMiceRipples.swProp{ii} = [];
end

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m},'\Summ'));  
    if exist('compiledRipples.mat','file') && ~force
        disp(['Loading ripples for mouse' mice{m}])
        load('compiledRipples.mat');
    else 
        disp(['Ripples not computed for mouse, calculating...' mice{m}])
        compiledRipples = compileRipples('expPath',strcat(parentDir, mice{m})); 
    end
    %Get number of sessions for the mouse to normalize
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    numSess = size(allSess,1);
    
    for ii = 1:(numAnalog+1)
        idxAnalog = (compiledRipples.psthAnalogEv ==ii);
        num =[];
        num = catpad(2,compiledRipples.rate_prestim(compiledRipples.rateAnalog == ii)',compiledRipples.rate_pre(compiledRipples.rateAnalog == ii)',...
            compiledRipples.rate_post(compiledRipples.rateAnalog == ii)',compiledRipples.rate_poststim(compiledRipples.rateAnalog == ii)');
        compiledMiceRipples.num{ii} = [compiledMiceRipples.num{ii}; num];
        
        amp =[];
        ampRipples = compiledRipples.peakAmplitude;
        amp = catpad(2,ampRipples(compiledRipples.ripple_prestim & idxAnalog)',...
            ampRipples(compiledRipples.ripple_pre & idxAnalog)',...
            ampRipples(compiledRipples.ripple_post & idxAnalog)',...
            ampRipples(compiledRipples.ripple_poststim & idxAnalog)');          
        compiledMiceRipples.amplitude{ii} = [compiledMiceRipples.amplitude{ii}; amp];
        
        dur =[];
        durRipples = compiledRipples.duration;
        dur = catpad(2,durRipples(compiledRipples.ripple_prestim & idxAnalog)',...
            durRipples(compiledRipples.ripple_pre & idxAnalog)',...
            durRipples(compiledRipples.ripple_post & idxAnalog)',...
            durRipples(compiledRipples.ripple_poststim & idxAnalog)');         
        compiledMiceRipples.duration{ii} = [compiledMiceRipples.duration{ii}; dur];
        
        freq =[];
        freqRipples = compiledRipples.peakFrequency;
        freq = catpad(2,freqRipples(compiledRipples.ripple_prestim & idxAnalog)',...
            freqRipples(compiledRipples.ripple_pre & idxAnalog)',...
            freqRipples(compiledRipples.ripple_post & idxAnalog)',...
            freqRipples(compiledRipples.ripple_poststim & idxAnalog)');         
        compiledMiceRipples.frequency{ii} = [compiledMiceRipples.frequency{ii}; freq];
        
        swMag =[];
        swMagRipples = compiledRipples.sWmag;
        swMag = catpad(2,swMagRipples(compiledRipples.ripple_prestim & idxAnalog)',...
            swMagRipples(compiledRipples.ripple_pre & idxAnalog)',...
            swMagRipples(compiledRipples.ripple_post & idxAnalog)',...
            swMagRipples(compiledRipples.ripple_poststim & idxAnalog)'); 
        compiledMiceRipples.swMag{ii} = [compiledMiceRipples.swMag{ii}; swMag];

        
        swProp =[];
        swProp = [sum(~isnan(swMagRipples(compiledRipples.ripple_prestim & idxAnalog)))./sum((compiledRipples.ripple_prestim & idxAnalog)),...
            sum(~isnan(swMagRipples(compiledRipples.ripple_pre & idxAnalog)))./sum((compiledRipples.ripple_pre & idxAnalog)),...
            sum(~isnan(swMagRipples(compiledRipples.ripple_post & idxAnalog)))./sum((compiledRipples.ripple_post & idxAnalog)),...
            sum(~isnan(swMagRipples(compiledRipples.ripple_poststim & idxAnalog)))./sum((compiledRipples.ripple_poststim & idxAnalog))];
        compiledMiceRipples.swProp{ii} = [compiledMiceRipples.swProp{ii}; swProp];        
    end
end

figure
%set(gcf,'Position',[50 50 1200 1000])
set(gcf,'Renderer','painters')

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    
    subplot(2,4,1) 
    dataAll{1} = compiledMiceRipples.num{2}(:,2);
    dataAll{2} = compiledMiceRipples.num{2}(:,3);    
    dataAll{3} = compiledMiceRipples.num{1}(:,3);    
    dataAll{4} = compiledMiceRipples.num{3}(:,3);    
    stats.ripplerate = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Number of ripples')
    
    subplot(2,4,2) 
    dataAll{1} = compiledMiceRipples.amplitude{2}(:,2);
    dataAll{2} = compiledMiceRipples.amplitude{2}(:,3);    
    dataAll{3} = compiledMiceRipples.amplitude{1}(:,3);    
    dataAll{4} = compiledMiceRipples.amplitude{3}(:,3);    
    stats.amplitude = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Ripple amplitude')    
    
    subplot(2,4,3) 
    dataAll{1} = compiledMiceRipples.duration{2}(:,2);
    dataAll{2} = compiledMiceRipples.duration{2}(:,3);    
    dataAll{3} = compiledMiceRipples.duration{1}(:,3);    
    dataAll{4} = compiledMiceRipples.duration{3}(:,3);    
    stats.duration = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Ripple duration')
    
    subplot(2,4,4) 
    dataAll{1} = compiledMiceRipples.frequency{2}(:,2);
    dataAll{2} = compiledMiceRipples.frequency{2}(:,3);    
    dataAll{3} = compiledMiceRipples.frequency{1}(:,3);    
    dataAll{4} = compiledMiceRipples.frequency{3}(:,3);    
    stats.frequency = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Ripple frequency')
    
    subplot(2,4,5) 
    dataAll{1} = compiledMiceRipples.swMag{2}(:,2);
    dataAll{2} = compiledMiceRipples.swMag{2}(:,3);    
    dataAll{3} = compiledMiceRipples.swMag{1}(:,3);    
    dataAll{4} = compiledMiceRipples.swMag{3}(:,3);    
    stats.swMag = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Sharp wave magnitude')
    
    subplot(2,4,6) 
    dataAll{1} = compiledMiceRipples.swProp{2}(:,2);
    dataAll{2} = compiledMiceRipples.swProp{2}(:,3);    
    dataAll{3} = compiledMiceRipples.swProp{1}(:,3);    
    dataAll{4} = compiledMiceRipples.swProp{3}(:,3);    
    stats.swProp = groupStats(dataAll,[],'inAxis',true);          
    ylabel('Sharp wave proportion')
    
else
    for ii = 1:(numAnalog+1)

        subplot((numAnalog+1),6,6*(ii-1)+1) 
        hold on
        
        dataAll{1} = compiledMiceRipples.num{ii}(:,1);
        dataAll{2} = compiledMiceRipples.num{ii}(:,2);
        dataAll{3} = compiledMiceRipples.num{ii}(:,3);
        dataAll{4} = compiledMiceRipples.num{ii}(:,4);       
            
        stats.ripplerate{ii} = groupStats(dataAll,[],'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesStd');        
        [stats.ripplerate{ii}.signrank.p,~,stats.ripplerate{ii}.signrank.stats] = signrank(compiledMiceRipples.num{ii}(:,2),compiledMiceRipples.num{ii}(:,3));
        title(num2str(stats.ripplerate{ii}.signrank.p));
        ylabel('Number of ripples')

        dataAll = [];
        subplot((numAnalog+1),6,6*(ii-1)+2)
        title('Ripple power')
        dataAll{1} = compiledMiceRipples.amplitude{ii}(:,2);
        dataAll{2} = compiledMiceRipples.amplitude{ii}(:,3);
        stats.amplitude{ii} = groupStats(dataAll,[],'inAxis',true);
        %nhist(dataAll,'samebins','smooth','pdf','median','median','sem','linewidth',1.5)
        [stats.amplitude{ii}.signrank.p,~,stats.amplitude{ii}.signrank.stats] = signrank(compiledMiceRipples.amplitude{ii}(:,2),compiledMiceRipples.amplitude{ii}(:,3));
        title(num2str(stats.amplitude{ii}.signrank.p));
        ylabel('Ripple Amplitude')

        dataAll = [];
        subplot((numAnalog+1),6,6*(ii-1)+3)
        title('Ripple duration')
       % dataAll{1} = compiledMiceRipples.duration{ii}(:,1);
        dataAll{1} = compiledMiceRipples.duration{ii}(:,2);
        dataAll{2} = compiledMiceRipples.duration{ii}(:,3);
       % dataAll{4} = compiledMiceRipples.duration{ii}(:,4);
        stats.duration{ii} = groupStats(dataAll,[],'inAxis',true);
        [stats.duration{ii}.signrank.p,~,stats.duration{ii}.signrank.stats] = signrank(compiledMiceRipples.duration{ii}(:,2),compiledMiceRipples.duration{ii}(:,3));
        title(num2str(stats.duration{ii}.signrank.p));        
        %nhist(dataAll,'samebins','smooth','pdf','median','median','sem','linewidth',1.5)
        ylabel('Ripple duration')
       % xlabel('Ripple duration')

        dataAll = [];
        subplot((numAnalog+1),6,6*(ii-1)+4)
        title('Ripple frequency')
        dataAll{1} = compiledMiceRipples.frequency{ii}(:,2);
        dataAll{2} = compiledMiceRipples.frequency{ii}(:,3);
        stats.frequency{ii} = groupStats(dataAll,[],'inAxis',true);
        [stats.frequency{ii}.signrank.p,~,stats.frequency{ii}.signrank.stats] = signrank(compiledMiceRipples.frequency{ii}(:,2),compiledMiceRipples.frequency{ii}(:,3));
        title(num2str(stats.frequency{ii}.signrank.p));  
        ylabel('Ripple frequency')
        
        dataAll = [];
        subplot((numAnalog+1),6,6*(ii-1)+5)
        title('SW magnitude')
        dataAll{1} = compiledMiceRipples.swMag{ii}(:,2);
        dataAll{2} = compiledMiceRipples.swMag{ii}(:,3);
        stats.swMag{ii} = groupStats(dataAll,[],'inAxis',true);
        [stats.swMag{ii}.signrank.p,~,stats.swMag{ii}.signrank.stats] = signrank(compiledMiceRipples.swMag{ii}(:,2),compiledMiceRipples.swMag{ii}(:,3));
        title(num2str(stats.swMag{ii}.signrank.p));  
        ylabel('SW magnitude')

        dataAll = [];
        subplot((numAnalog+1),6,6*(ii-1)+6)
        title('SW proportion')
        dataAll{1} = compiledMiceRipples.swProp{ii}(:,2);
        dataAll{2} = compiledMiceRipples.swProp{ii}(:,3);
        stats.swProp{ii} = groupStats(dataAll,[],'inAxis',true);
        [stats.swProp{ii}.signrank.p,~,stats.swProp{ii}.signrank.stats] = signrank(compiledMiceRipples.swProp{ii}(:,2),compiledMiceRipples.swProp{ii}(:,3));
        title(num2str(stats.swProp{ii}.signrank.p));  
        ylabel('SW proportion')        
    end    
end

if ~Control
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.png'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.fig'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.eps'));
    save(strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.mat'),'stats')
else
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.png'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.fig'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.eps'));
    save(strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.mat'),'stats')
end
end