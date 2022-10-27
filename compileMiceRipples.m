function compiledMiceRipples = compileMiceRipples

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
    compiledMiceRipples.num{ii} = [];
    compiledMiceRipples.numpost{ii} = [];
    compiledMiceRipples.amplitude{ii} = [];
    compiledMiceRipples.duration{ii} = [];
    compiledMiceRipples.ISI{ii} = [];
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
            compiledRipples.rate_post(compiledRipples.rateAnalog == ii)',compiledRipples.rate_poststim(compiledRipples.rateAnalog == ii)',...
            compiledRipples.rate_postpoststim(compiledRipples.rateAnalog == ii)');
        compiledMiceRipples.num{ii} = [compiledMiceRipples.num{ii}; num];
        
        numpost = [];
        numpost = catpad(2,compiledRipples.rate_poststim1(compiledRipples.rateAnalog == ii)',compiledRipples.rate_poststim2(compiledRipples.rateAnalog == ii)',...
            compiledRipples.rate_poststim3(compiledRipples.rateAnalog == ii)',compiledRipples.rate_poststim4(compiledRipples.rateAnalog == ii)',...
            compiledRipples.rate_poststim5(compiledRipples.rateAnalog == ii)');
        compiledMiceRipples.numpost{ii} = [compiledMiceRipples.numpost{ii}; numpost];
        
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

        
        compiledMiceRipples.swProp{ii} = [compiledMiceRipples.swProp{ii}; compiledRipples.sWprop(compiledRipples.rateAnalog == ii,:)];      
        
        ISI =[];
        ISIRipples = compiledRipples.ISI;
        ISI = catpad(2,ISIRipples(compiledRipples.ripple_prestim & idxAnalog)',...
            ISIRipples(compiledRipples.ripple_pre & idxAnalog)',...
            ISIRipples(compiledRipples.ripple_post & idxAnalog)',...
            ISIRipples(compiledRipples.ripple_poststim & idxAnalog)');          
        compiledMiceRipples.ISI{ii} = [compiledMiceRipples.ISI{ii}; ISI];        
    end
end

figure
set(gcf,'Position',[50 50 1200 300])
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
colMat = [85/243 85/243 85/243; 85/243 85/243 85/243; 8/243 133/243 161/243; 122/243 122/243 122/243];%; 122/243 122/243 122/243];
col = [85/243 85/243 85/243; 8/243 133/243 161/243];
colormap(col)

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
    
elseif strcmp(tag,'mEC') == 1 
    
    set(gcf,'Position',[437 70 1385 800])
    subplot(2,5,1) 
    %subplot(1,4,1)
    dataAll = [];
    dataAll{1} = compiledMiceRipples.num{2}(:,1);
    dataAll{2} = compiledMiceRipples.num{2}(:,2);    
    dataAll{3} = compiledMiceRipples.num{2}(:,3);    
    dataAll{4} = compiledMiceRipples.num{2}(:,4);    
    %dataAll{5} = compiledMiceRipples.num{2}(:,5);  
    stats.ripplerate = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);      
   [ stats.ripplerate.signrank.p,~,stats.ripplerate.signrank.stats] = signrank(compiledMiceRipples.num{2}(:,2),compiledMiceRipples.num{2}(:,3));
    title(num2str(stats.ripplerate.signrank.p));
    ylabel('Number of ripples')
    ylim([0 0.8])
       
    %subplot(1,4,2)
    subplot(2,5,2) 
    hold off
    dataAll = [];
    %dataAll{1} = compiledMiceRipples.amplitude{2}(:,1);
    dataAll{1} = compiledMiceRipples.amplitude{2}(:,2);    
    dataAll{2} = compiledMiceRipples.amplitude{2}(:,3);    
    dataAll{1}(dataAll{1}>2500) = nan;
    dataAll{2}(dataAll{2}>2500) = nan;
    h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', .25, 'dot_dodge_amount', .3,...
     'box_col_match', 0);
    set(h1{2}, 'SizeData', 2);
    h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
    set(h2{2}, 'SizeData', 2);
    box off
    [stats.amplitude.signrank.p,~,stats.amplitude.signrank.stats] = ranksum(compiledMiceRipples.amplitude{2}(:,2),compiledMiceRipples.amplitude{2}(:,3));
    title(num2str(stats.amplitude.signrank.p));    
    stats.amplitude = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false); 
    set(gca,'view',[90 -90])
    xlabel('Ripple amplitude')    
    xlim([0 3000])
    
    %subplot(1,4,3)
    subplot(2,5,3) 
    dataAll = [];
    dataAll{1} = compiledMiceRipples.duration{2}(:,2);    
    dataAll{2} = compiledMiceRipples.duration{2}(:,3); 
%     dataAll{1}(dataAll{1}>2500) = nan;
%     dataAll{2}(dataAll{2}>2500) = nan;
    h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', .25, 'dot_dodge_amount', .3,...
     'box_col_match', 0);
    set(h1{2}, 'SizeData', 2);
    h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
    set(h2{2}, 'SizeData', 2);
    box off    
    %dataAll{4} = compiledMiceRipples.duration{2}(:,4);    
    %nhist(dataAll(2:3), 'binfactor',1,'color','colormap','samebins','proportion','serror','linewidth',1.5)
    [stats.duration.signrank.p,~,stats.duration.signrank.stats] = ranksum(compiledMiceRipples.duration{2}(:,2),compiledMiceRipples.duration{2}(:,3));
    title(num2str(stats.duration.signrank.p));        
    stats.duration = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false);          
    set(gca,'view',[90 -90])
    xlabel('Ripple duration')
    xlim([0 0.12])
    
    %subplot(1,4,4)
    subplot(2,5,4) 
    dataAll = [];
    dataAll{1} = compiledMiceRipples.frequency{2}(:,2);    
    dataAll{2} = compiledMiceRipples.frequency{2}(:,3);    
%     dataAll{1}(dataAll{1}>2500) = nan;
%     dataAll{2}(dataAll{2}>2500) = nan;
    h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', .25, 'dot_dodge_amount', .3,...
     'box_col_match', 0);
    set(h1{2}, 'SizeData', 2);
    h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
    set(h2{2}, 'SizeData', 2);
    box off       
    [stats.frequency.signrank.p,~,stats.frequency.signrank.stats] = ranksum(compiledMiceRipples.frequency{2}(:,2),compiledMiceRipples.frequency{2}(:,3));
    title(num2str(stats.frequency.signrank.p));         
    stats.frequency = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false);     
    set(gca,'view',[90 -90])
    xlabel('Ripple frequency')
    xlim([130 200])
       
    subplot(2,5,5) 
    dataAll = [];
    dataAll{1} = compiledMiceRipples.swMag{2}(:,2);    
    dataAll{2} = compiledMiceRipples.swMag{2}(:,3);       
%     dataAll{1}(dataAll{1}>2500) = nan;
%     dataAll{2}(dataAll{2}>2500) = nan;
    h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', .25, 'dot_dodge_amount', .3,...
     'box_col_match', 0);
    set(h1{2}, 'SizeData', 2);
    h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
     'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
    set(h2{2}, 'SizeData', 2);
    box off          
    stats.swMag = groupStats(dataAll,[],'inAxis',true,'plotType','violinPlot','color',col,'doPlot',false);     
    [stats.swMag.signrank.p,~,stats.swMag.signrank.stats] = ranksum(compiledMiceRipples.swMag{2}(:,2),compiledMiceRipples.swMag{2}(:,3));
    title(num2str(stats.swMag.signrank.p));     
    set(gca,'view',[90 -90])
    xlabel('Sharp wave magnitude')
    xlim([-10 0])
    
    subplot(2,5,6) 
    swMag1 = compiledMiceRipples.swMag{2}(:,2);     
    swMag2 = compiledMiceRipples.swMag{2}(:,3);     
    ripAmp1 = compiledMiceRipples.amplitude{2}(:,2);    
    ripAmp2 = compiledMiceRipples.amplitude{2}(:,3);  
   
    scatter(abs(swMag1),ripAmp1,15,col(1,:),'.')
    hold on
    scatter(abs(swMag2),ripAmp2, 15, col(2,:),'.') 
    [stats.swMagripAmp.Pre.R,stats.swMagripAmp.Pre.p]  = corr(abs(swMag1),ripAmp1);
    [stats.swMagripAmp.Post.R,stats.swMagripAmp.Post.p]  = corr(abs(swMag2),ripAmp2);
    ylabel('Ripple amplitude') 
    xlabel('SPW amplitude') 

    subplot(2,5,7) 
    dataAll = [];
    dataAll{1} = compiledMiceRipples.swProp{2}(:,2);    
    dataAll{2} = compiledMiceRipples.swProp{2}(:,3);  
%     dataAll{1}(dataAll{1}<0.0267) = nan;
%     dataAll{2}(dataAll{2}<0.0267) = nan;    
    stats.swProp = groupStats(dataAll,[],'inAxis',true,'plotType','boxplot','color',col); 
    hold on
    plot(compiledMiceRipples.swProp{2}(:,2:3)')
    [stats.swProp.signrank.p,~,stats.swProp.signrank.stats] = signrank(dataAll{1},dataAll{2});
    title(num2str(stats.swProp.signrank.p));      
    xlabel('Sharp wave proportion') 
    
    subplot(2,5,9) 
    dataAll = [];
    dataAll{1} = compiledMiceRipples.numpost{2}(:,1);
    dataAll{2} = compiledMiceRipples.numpost{2}(:,2);    
    dataAll{3} = compiledMiceRipples.numpost{2}(:,3);    
    dataAll{4} = compiledMiceRipples.numpost{2}(:,4);    
    dataAll{5} = compiledMiceRipples.numpost{2}(:,5);  
    stats.rippleratepost = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat(1,:));      
    ylabel('Number of ripples in recovery')
%     
%     dataAll = [];
%     dataAll(:,1) = compiledMiceRipples.numpost{2}(:,1);
%     dataAll(:,2) = compiledMiceRipples.numpost{2}(:,2);    
%     dataAll(:,3) = compiledMiceRipples.numpost{2}(:,3);    
%     dataAll(:,4) = compiledMiceRipples.numpost{2}(:,4);    
%     dataAll(:,5) = compiledMiceRipples.numpost{2}(:,5);  
%     
else
    set(gcf,'Position',[35 150 1800 680])
    for ii = 1:(numAnalog+1)

        subplot((numAnalog+1),8,8*(ii-1)+1) 
        hold on
        dataAll = [];
        dataAll{1} = compiledMiceRipples.num{ii}(:,1);
        dataAll{2} = compiledMiceRipples.num{ii}(:,2);
        dataAll{3} = compiledMiceRipples.num{ii}(:,3);
        dataAll{4} = compiledMiceRipples.num{ii}(:,4);       
            
        stats.ripplerate{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat);        
        [stats.ripplerate{ii}.signrank.p,~,stats.ripplerate{ii}.signrank.stats] = signrank(compiledMiceRipples.num{ii}(:,2),compiledMiceRipples.num{ii}(:,3));
        title(num2str(stats.ripplerate{ii}.signrank.p));
        ylabel('Number of ripples')
        ylim([0 0.8])
        
        subplot((numAnalog+1),8,8*(ii-1)+2)
        title('Ripple power')
        dataAll = [];
        dataAll{1} = compiledMiceRipples.amplitude{ii}(:,2);    
        dataAll{2} = compiledMiceRipples.amplitude{ii}(:,3);    
        dataAll{1}(dataAll{1}>2500) = nan;
        dataAll{2}(dataAll{2}>2500) = nan;
        h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0.25, 'dot_dodge_amount', .3,...
         'box_col_match', 0);
        set(h1{2}, 'SizeData', 2);
        h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
        set(h2{2}, 'SizeData', 2);
        box off    
        stats.amplitude{ii} = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false);  
        [stats.amplitude{ii}.signrank.p,~,stats.amplitude{ii}.signrank.stats] = ranksum(compiledMiceRipples.amplitude{ii}(:,2),compiledMiceRipples.amplitude{ii}(:,3));
        title(num2str(stats.amplitude{ii}.signrank.p));        
        set(gca,'view',[90 -90])
        xlabel('Ripple amplitude')    
        xlim([0 3000])        
        
        dataAll = [];
        subplot((numAnalog+1),8,8*(ii-1)+3)
        title('Ripple duration')
        dataAll{1} = compiledMiceRipples.duration{ii}(:,2);
        dataAll{2} = compiledMiceRipples.duration{ii}(:,3);
        h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0.25, 'dot_dodge_amount', .3,...
         'box_col_match', 0);
        set(h1{2}, 'SizeData', 2);
        h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
        set(h2{2}, 'SizeData', 2);
        box off    
        stats.duration{ii} = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false);  
        [stats.duration{ii}.signrank.p,~,stats.duration{ii}.signrank.stats] = ranksum(compiledMiceRipples.duration{ii}(:,2),compiledMiceRipples.duration{ii}(:,3));
        title(num2str(stats.duration{ii}.signrank.p));         
        set(gca,'view',[90 -90])
        xlabel('Ripple duration')
        xlim([0 0.12])

        dataAll = [];
        subplot((numAnalog+1),8,8*(ii-1)+4)
        title('Ripple frequency')
        dataAll{1} = compiledMiceRipples.frequency{ii}(:,2);
        dataAll{2} = compiledMiceRipples.frequency{ii}(:,3);
        h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:),'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0.25, 'dot_dodge_amount', .3,...
         'box_col_match', 0);
        set(h1{2}, 'SizeData', 2);
        h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
        set(h2{2}, 'SizeData', 2);
        box off    
        stats.frequency{ii} = groupStats(dataAll,[],'plotType','violinPlot','inAxis',true,'color',col,'doPlot',false);        
        [stats.frequency{ii}.signrank.p,~,stats.frequency{ii}.signrank.stats] = ranksum(compiledMiceRipples.frequency{ii}(:,2),compiledMiceRipples.frequency{ii}(:,3));
        title(num2str(stats.frequency{ii}.signrank.p));            
        set(gca,'view',[90 -90])
        xlabel('Ripple frequency')
        xlim([130 200])
        
        dataAll = [];
        subplot((numAnalog+1),8,8*(ii-1)+5)
        title('SW magnitude')
        dataAll{1} = compiledMiceRipples.swMag{ii}(:,2);
        dataAll{2} = compiledMiceRipples.swMag{ii}(:,3);        
        h1 = raincloud_plot(dataAll{1}, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0.25, 'dot_dodge_amount', .3,...
         'box_col_match', 0);
        set(h1{2}, 'SizeData', 2);
        h2 = raincloud_plot(dataAll{2}, 'box_on', 1, 'color', col(2,:), 'alpha', 0.4,...
         'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', .05, 'box_col_match', 0);
        set(h2{2}, 'SizeData', 2);
        box off          
        stats.swMag{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','violinPlot','color',col,'doPlot',false);     
        [stats.swMag{ii}.signrank.p,~,stats.swMag{ii}.signrank.stats] = ranksum(compiledMiceRipples.swMag{ii}(:,2),compiledMiceRipples.swMag{ii}(:,3));
        title(num2str(stats.swMag{ii}.signrank.p));          
        set(gca,'view',[90 -90])
        xlabel('Sharp wave magnitude')
        xlim([-10 0])
        
        %% Only count those sessions where there were at least 4 ripples
        dataAll = [];
        subplot((numAnalog+1),8,8*(ii-1)+6)
        title('SW proportion')
        dataAll{1} = compiledMiceRipples.swProp{ii}(:,2);
        dataAll{2} = compiledMiceRipples.swProp{ii}(:,3); 
%         dataAll{1}(compiledMiceRipples.num{ii}(:,2)<0.015) = nan;
%         dataAll{2}(compiledMiceRipples.num{ii}(:,3)<0.015) = nan;
        stats.swProp{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','boxplot','color',col);  
        hold on
        plot([dataAll{1} dataAll{2}]')
        [stats.swProp{ii}.signrank.p,~,stats.swProp{ii}.signrank.stats] = signrank(dataAll{1},dataAll{2});
        title(num2str(stats.swProp{ii}.signrank.p));   
        xlabel('Sharp wave proportion') 
    
        subplot((numAnalog+1),8,8*(ii-1)+7)
        swMag1 = compiledMiceRipples.swMag{ii}(:,2);     
        swMag2 = compiledMiceRipples.swMag{ii}(:,3);     
        ripAmp1 = compiledMiceRipples.amplitude{ii}(:,2);    
        ripAmp2 = compiledMiceRipples.amplitude{ii}(:,3);  

        scatter(abs(swMag1),ripAmp1,45,col(1,:),'.')
        hold on
        scatter(abs(swMag2),ripAmp2, 45, col(2,:),'.') 
        [stats.swMagripAmp.Pre.R,stats.swMagripAmp.Pre.p]  = corr(abs(swMag1),ripAmp1);
        [stats.swMagripAmp.Post.R,stats.swMagripAmp.Post.p]  = corr(abs(swMag2),ripAmp2);
        ylabel('Ripple amplitude') 
        xlabel('SPW amplitude') 
        ylim([500 2000])
        
        subplot((numAnalog+1),8,8*(ii-1)+8)
        dataAll = [];
        dataAll{1} = compiledMiceRipples.numpost{ii}(:,1);
        dataAll{2} = compiledMiceRipples.numpost{ii}(:,2);    
        dataAll{3} = compiledMiceRipples.numpost{ii}(:,3);    
        dataAll{4} = compiledMiceRipples.numpost{ii}(:,4);    
        dataAll{5} = compiledMiceRipples.numpost{ii}(:,5);  
        stats.rippleratepost{ii} = groupStats(dataAll,[],'inAxis',true,'plotType','boxLinesSEM','repeatedMeasures',true,'color',colMat(1,:));      
        ylabel('Number of ripples in recovery')
    end    

end

if ~Control
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.png'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.fig'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.eps'),'epsc');
    save(strcat(parentDir,'Compiled\Ripples\compiledRipples',tag,'.mat'),'stats')
else
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.png'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.fig'));
    saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.eps'),'epsc');
    save(strcat(parentDir,'Compiled\Ripples\compiledRipplesCtrl',tag,'.mat'),'stats')
end
end