function sessRippleCortexDS

%Quantify number of ripples around each manipulation
%Look at distribution of ripples around downstates for all, pre, stim, post

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cond = 1; % 0 is mEC only, 1 is mEC bilateral
downUP = 0; % 0 is UP-DOWN transition, 1 if DOWN-UP transition
down = 1; %1 if looking at down states, 0 if looking at up states
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);

if cond == 0
    pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
    'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};
    analogCh = 2;
elseif cond == 1
    pathToSessionsAll = {'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3'};
    analogCh = [1 2 3];
end

distRipAll = zeros(length(pathToSessionsAll),101);
distDSDurAll  = zeros(length(pathToSessionsAll),50);

for ii = analogCh
    
    distRipPulses{ii} = zeros(length(pathToSessionsAll),51);
    
    distRipPre{ii} = zeros(length(pathToSessionsAll),35);
    distRipStim{ii} = zeros(length(pathToSessionsAll),35);
    distRipPost{ii} = zeros(length(pathToSessionsAll),35);
    
    ratioRipPre{ii} = zeros(length(pathToSessionsAll),1);
    ratioRipStim{ii} = zeros(length(pathToSessionsAll),1);
    ratioRipPost{ii} = zeros(length(pathToSessionsAll),1);
    
    distDurPre{ii} = [];
    distDurStim{ii} = [];
    distDurPost{ii} = [];
    
    distNumPre{ii} = zeros(length(pathToSessionsAll),1);
    distNumStim{ii} = zeros(length(pathToSessionsAll),1);
    distNumPost{ii} = zeros(length(pathToSessionsAll),1);    
end

for ii=1:size(pathToSessionsAll,2)
    
    cd(strcat(expPath,pathToSessionsAll{ii}));
    file = dir('*.SlowWaves.events.mat');
    load(file.name);    
    if downUP ==0
        st = SlowWaves.ints.DOWN(:,1);
    elseif downUP == 1
        st = SlowWaves.ints.DOWN(:,2);
    end
    
    if down == 1
        swINTS = SlowWaves.ints.DOWN;
    else
        swINTS = SlowWaves.ints.UP;
        st = SlowWaves.ints.UP(:,2);
    end
    
    file = dir('*.ripples.events.mat');
    load(file.name);    
    
    file = dir('*.SleepState.states.mat');
    load(file.name);   
        
    % Load pulses
    disp('Getting analog-in inputs...');
    [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
    
    %% First look at distribution across all DownStates and ripples
    [stccg, tRip] = CCG({ripples.peaks st},[],'binSize',0.01,'duration',1,'norm','rate');
    distRipAll(ii,:) = stccg(:,2,1);
    % average duration of down states 
    distDSDurAll(ii,:) = histcounts(swINTS(:,2)-swINTS(:,1),0:0.02:1, 'Normalization', 'probability');
    
    for i = analogCh
        keepIdx = [];
        
        if cond ==1
            if exist('pulses')
                if i<=2
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
            end   
        elseif cond==0
            if ii< 16
                if i<=2
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
            else
                pulTr = (pulses.stimPerID'==1 & (pulses.stimComb==2 | pulses.stimComb==3)); % for the CA3 mice, both 2 and 3 are the same
            end
        end
        events = pulses.intsPeriods(:,pulTr)';       
        
        keepIdx(:,1) = InIntervals(st,events-5);    
        keepIdx(:,2) = InIntervals(st,events);           
        keepIdx(:,3) = InIntervals(st,events+5);       
        keepIdx = logical(keepIdx);
        
        if downUP ==0
            events = pulses.intsPeriods(1,pulTr)'; 
            events_st = events(InIntervals(events,SleepState.ints.NREMstate));
        elseif downUP==1
            events = pulses.intsPeriods(2,pulTr)';  
            events_st = events(InIntervals(events,SleepState.ints.NREMstate));
        end
        
        [stccg, tPul] = CCG({ripples.peaks events_st},[],'binSize',0.02,'duration',1,'norm','rate');
        distRipPulses{i}(ii,:) = stccg(:,2,1);
    
        %% Also confirm that the number and duration of down states remains the same
        distDurPre{i} = [distDurPre{i};swINTS(keepIdx(:,1),2)-swINTS(keepIdx(:,1),1)];
        distDurStim{i} = [distDurStim{i};swINTS(keepIdx(:,2),2)-swINTS(keepIdx(:,2),1)];
        distDurPost{i} = [distDurPost{i};swINTS(keepIdx(:,3),2)-swINTS(keepIdx(:,3),1)];
        
        distNumPre{i}(ii,1) = sum(keepIdx(:,1))./(size(events,1)*5);
        distNumStim{i}(ii,1) = sum(keepIdx(:,2))./(size(events,1)*5);
        distNumPost{i}(ii,1) = sum(keepIdx(:,3))./(size(events,1)*5);           

        if sum(keepIdx(:,1))>10
            [stccg, t] = CCG({ripples.peaks st(keepIdx(:,1)) },[],'binSize',0.03,'duration',1,'norm','rate');            
            distRipPre{i}(ii,:) = stccg(:,2,1);
        else
            distRipPre{i}(ii,:) = nan;     
        end
        
        if sum(keepIdx(:,2))>10
            [stccg, t] = CCG({ripples.peaks st(keepIdx(:,2)) },[],'binSize',0.03,'duration',1,'norm','rate');            
            distRipStim{i}(ii,:) = stccg(:,2,1);
        else
            distRipStim{i}(ii,:) = nan;     
        end        

        if sum(keepIdx(:,3))>10
            [stccg, t] = CCG({ripples.peaks st(keepIdx(:,3)) },[],'binSize',0.03,'duration',1,'norm','rate');            
            distRipPost{i}(ii,:) = stccg(:,2,1);
        else
            distRipPost{i}(ii,:) = nan;     
        end
        
    end
    
end

figure
set(gcf,'Color','w')
set(gcf, 'Renderer','painters')
set(gcf,'Position',[1 41 1920 970])
subplot(5,6,1)
imagesc(tRip,1:1:size(distRipAll,1),zscore(distRipAll,[],2))
colormap(YlGnBu)
colorbar    
caxis([-2 2])

subplot(5,6,7)
meanpsth = nanmean(distRipAll,1);
stdpsth = nanstd(distRipAll,1)./sqrt(size(distRipAll,1));        
% meanpsth = nanmean(zscore(distRipAll,[],2),1);
% stdpsth = nanstd(zscore(distRipAll,[],2),1)./sqrt(size(distRipAll,1));               
hold on
fill([tRip; flipud(tRip)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(tRip,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
title('Ripple dist')

subplot(5,6,13)
%imagesc(tPul,1:1:size(distRipPulses{2},1),zscore(distRipPulses{2},[],2))
imagesc(tPul,1:1:size(distRipPulses{2},1),distRipPulses{2})
colormap(YlGnBu)
colorbar    
caxis([-2 2])

subplot(5,6,19)
% meanpsth = nanmean(zscore(distRipPulses{2},[],2),1);
% stdpsth = nanstd(zscore(distRipPulses{2},[],2),1)./sqrt(size(distRipPulses{2},1));       
meanpsth = nanmean(distRipPulses{2});
stdpsth = nanstd(distRipPulses{2})./sqrt(size(distRipPulses{2},1));               
hold on
fill([tPul; flipud(tPul)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(tPul,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
title('Ripple dist')

a = 0:0.02:1;
a = a(1:50);
subplot(5,6,25)
meanpsth = nanmean(distDSDurAll,1);
stdpsth = nanstd(distDSDurAll,1)./sqrt(size(distDSDurAll,1));  
hold on
fill([a'; flipud(a')],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(a,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
title('Down state duration')

for ii = analogCh
    subplot(5,6,ii+1)
    aa = ~isnan(distRipPre{ii}(:,1));
    RipPre = distRipPre{ii}(aa,:);
    %imagesc(t,1:1:size(distRipPre{ii},1),zscore(RipPre,[],2))
    imagesc(t,1:1:size(distRipPre{ii},1),RipPre)
    colormap(YlGnBu)
    colorbar    
    caxis([-2 4])
    
    subplot(5,6,ii+7)
    aa = ~isnan(distRipStim{ii}(:,1));    
    RipStim = distRipStim{ii}(aa,:);    
    %imagesc(t,1:1:size(distRipStim{ii},1),zscore(RipStim,[],2))
    imagesc(t,1:1:size(distRipStim{ii},1),RipStim)
    colormap(YlGnBu)
    colorbar  
    caxis([-2 4])
    
    subplot(5,6,ii+13)
    aa = ~isnan(distRipPost{ii}(:,1));    
    RipPost = distRipPost{ii}(aa,:);     
    %imagesc(t,1:1:size(distRipPost{ii},1),zscore(RipPost,[],2))
    imagesc(t,1:1:size(distRipPost{ii},1),RipPost)
    colormap(YlGnBu)
    colorbar  
    caxis([-2 4])
    
    subplot(5,6,ii+19)
%     meanpsth = nanmean(zscore(distRipPre{ii},[],2),1);
%     stdpsth = nanstd(zscore(distRipPre{ii},[],2),1)./sqrt(size(distRipPre{ii},1));               
    meanpsth = nanmean(distRipPre{ii},1);
    stdpsth = nanstd(distRipPre{ii},1)./sqrt(size(distRipPre{ii},1));            
    hold on
    fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0 0 0],'linestyle','none','FaceAlpha',0.5);                    
    hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0 0 0]);    

%     meanpsth = nanmean(zscore(distRipPost{ii},[],2),1);
%     stdpsth = nanstd(zscore(distRipPost{ii},[],2),1)./sqrt(size(distRipPost{ii},1));          
    meanpsth = nanmean(distRipPost{ii},1);
    stdpsth = nanstd(distRipPost{ii},1)./sqrt(size(distRipPost{ii},1)); 
    hold on
    fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
    hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);    

%     meanpsth = nanmean(zscore(distRipStim{ii},[],2),1);
%     stdpsth = nanstd(zscore(distRipStim{ii},[],2),1)./sqrt(size(distRipStim{ii},1));     
    meanpsth = nanmean(distRipStim{ii},1);
    stdpsth = nanstd(distRipStim{ii},1)./sqrt(size(distRipStim{ii},1)); 
    hold on
    fill([t; flipud(t)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0 0.5 1],'linestyle','none','FaceAlpha',0.5);                    
    hi = line(t,meanpsth,'LineWidth',1.5,'Color',[0 0.5 1]);    
    ylim([0 1.5])
%     distDurPre{ii}(distDurPre{ii}>1.2) = nan;
%     distDurStim{ii}(distDurStim{ii}>1.2) = nan;
%     distDurPost{ii}(distDurPost{ii}>1.2) = nan;
    
    subplot(5,6,(6*(ii-1))+5)  
    stats.Dur{ii} = groupStats([{distDurPre{ii}} {distDurStim{ii}} {distDurPost{ii}}],[1 2 3],'plotType','boxplot','inAxis',true);     

    subplot(5,6,(6*(ii-1))+6) % number of downstate events is similar
    stats.Num{ii} = groupStats([{distNumPre{ii}} {distNumStim{ii}} {distNumPost{ii}}],[1 2 3],'inAxis',true,'plotType','BoxLinesSEM','repeatedMeasures',true);
    
    subplot(5,6,25+ii) 
    preidx = t>-0.1 & t<=0;
    ratPrePre = nanmean(distRipPre{ii}(:,preidx),2);
    ratStimPre = nanmean(distRipStim{ii}(:,preidx),2);
    ratPostPre = nanmean(distRipPost{ii}(:,preidx),2);
    stats.ripPre{ii} = groupStats([{ratPrePre} {ratStimPre} {ratPostPre}],[1 2 3],'inAxis',true,'plotType','BoxLinesSEM','repeatedMeasures',true);
    
end

saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.png'));
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.fig'));
save(strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.mat'),'stats');

end