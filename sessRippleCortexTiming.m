function sessRippleCortexTiming

% How long after a DOWN-UP transition does a ripple happen?

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cond = 1; % 0 is mEC only, 1 is mEC bilateral
downUP = 1; % 0 is UP-DOWN transition, 1 if DOWN-UP transition
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

distRipAll = zeros(length(pathToSessionsAll),35);
distDSDurAll  = zeros(length(pathToSessionsAll),50);

for ii = analogCh
    distRipPre{ii} = [];
    distRipStim{ii} = [];
    distRipPost{ii} = [];
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
    
    file = dir('*.ripples.events.mat');
    load(file.name);    
    
    % Load pulses
    disp('Getting analog-in inputs...');
    [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);

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
        
        keepIdx(:,1) = InIntervals(ripples.peaks,events-5);    
        keepIdx(:,2) = InIntervals(ripples.peaks,events);           
        keepIdx(:,3) = InIntervals(ripples.peaks,events+5);       
        keepIdx = logical(keepIdx);
        
        for ss = 1:3            
            if sum(keepIdx(:,ss))>1
                riptoKeep = ripples.peaks(keepIdx(:,ss));
                for rr = 1:length(riptoKeep)
                    %Find DOWN UP transition immediately preceding ripple
                    tDiff = riptoKeep(rr)-st;
                    tLastDiff = find(tDiff>0,1,'last');
                    if isempty(tLastDiff) || tDiff(tLastDiff)>0.5
                        continue
                    end
                    if ss == 1
                        distRipPre{i} = [distRipPre{i} tDiff(tLastDiff)];
                    elseif ss == 2
                        distRipStim{i} = [distRipStim{i} tDiff(tLastDiff)];
                    elseif ss == 3
                        distRipPost{i} = [distRipPost{i} tDiff(tLastDiff)];
                    end 
                end   
            end
        end
    end   
end


figure
set(gcf,'Color','w')
set(gcf, 'Renderer','painters')
set(gcf,'Position',[1 41 1920 970])

for ii = analogCh
    dist{1} = distRipPre{ii};
    dist{2} = distRipStim{ii};
    dist{3} = distRipPost{ii};
    
    subplot(3,2,2*(ii-1)+1)
    nhist(dist,'proportion','samebins','binfactor',6,'color','sequential')
    subplot(3,2,2*(ii-1)+2)
    stats{ii} =  groupStats(dist,[1 2 3],'inAxis',true);
end
    

%     ratPrePre = nanmean(distRipPre{ii}(:,preidx),2);
%     ratStimPre = nanmean(distRipStim{ii}(:,preidx),2);
%     ratPostPre = nanmean(distRipPost{ii}(:,preidx),2);
%     stats.ripPre{ii} = groupStats([{ratPrePre} {ratStimPre} {ratPostPre}],[1 2 3],'inAxis',true,'plotType','BoxLinesSEM','repeatedMeasures',true);
% 
% saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.png'));
% saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.eps'),'epsc');
% saveas(gcf,strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.fig'));
% save(strcat(expPath,'Compiled\Ripples\DownState\RippleDS_',num2str(cond),'_',num2str(downUP),'_',num2str(down),'.mat'),'stats');

end