function compileMiceDentateSpikes

tag = 'mEC';% mEC, CA3, Bilateral mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ20\Final','IZ21\Final','IZ30\Final',};
    reg = {'CA1','mEC','CA1Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ20\Final','IZ21\Final','IZ27\Saline',...
    'IZ28\Saline','IZ33\Saline'}; %'IZ30\Final','IZ29\Saline','IZ32\Saline',
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ29\Final','IZ32\Final','IZ33\Final','IZ27\Final','IZ28\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ33\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Bilateral'};
end

numAnalog = 2;
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
cmap = cbrewer('qual','Pastel2',length(mice));

for ii = 1:3
    DS.number{ii} = [];
    DS.amplitude{ii,1} = [];
    DS.amplitude{ii,2} = [];
    DS.amplitude{ii,3} = [];
end

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m}));  
    allSess = dir('*_sess*');
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));        
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        %Load DS        
        if exist([sessionInfo.FileName '.DS2.events.mat'],'file') 
            load([sessionInfo.FileName '.DS2.events.mat']);
        else
            disp('First calculate .DS2 file! Skipping');
            continue
        end
        
        % Load pulses
        load([sessionInfo.FileName '.pulses.events.mat']);
            
        for i = 1:(numAnalog+1)
            keepIdx = [];
            if exist('pulses')
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
            end
            
            %Only select the pulses that happened in the home cage,
            %i.e., which were 5 seconds long
            homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
            homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
            pulTr = pulTr & homeCagePulseidx;
            events = pulses.intsPeriods(:,pulTr)';
            
            keepIdx(:,1) = InIntervals(DS2.peaks,events);    
            keepIdx(:,2) = InIntervals(DS2.peaks,(events-5));    
            keepIdx(:,3) = InIntervals(DS2.peaks,(events+5));  
            
            DS.number{i} = [DS.number{i}; sum(keepIdx(:,2))./(size(events,1)*5) sum(keepIdx(:,1))./(size(events,1)*5) sum(keepIdx(:,3))./(size(events,1)*5)];
            DS.amplitude{i,1} = [DS.amplitude{i,1}; DS2.amplitudes(logical(keepIdx(:,2)))];
            DS.amplitude{i,2} = [DS.amplitude{i,2}; DS2.amplitudes(logical(keepIdx(:,1)))];
            DS.amplitude{i,3} = [DS.amplitude{i,3}; DS2.amplitudes(logical(keepIdx(:,3)))];
        end

    end    
end

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    figure    
    colMat = [0.5 0.5 0.5;8/243 133/243 161/243;193/243 90/243 99/243; 133/243 128/243 177/243];        
    set(gcf,'Position',[50 50 1200 800])
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    subplot(1,2,1)
    dataAll  = [DS.number{2}(:,1);DS.number{2}(:,2);DS.number{3}(:,1);DS.number{3}(:,2)];
    dataID = [ones(size(DS.number{2},1),1)*1;ones(size(DS.number{2},1),1)*2;ones(size(DS.number{2},1),1)*3;ones(size(DS.number{2},1),1)*4];
    stats.DSrate = groupStats(dataAll,dataID,'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM','color',colMat);    
    title('Dentate Spike 2 Rate')
    
    subplot(1,2,2)
    dataAll  = [DS.amplitude{2,1};DS.amplitude{2,2};DS.amplitude{3,1};DS.amplitude{3,2}];
    dataID = [ones(length(DS.amplitude{2,1}),1)*1;ones(length(DS.amplitude{2,2}),1)*2;ones(length(DS.amplitude{3,1}),1)*3;ones(length(DS.amplitude{3,2}),1)*4];
    stats.DSrate = groupStats(dataAll,dataID,'inAxis',true,'color',colMat);    
    title('Dentate Spike 2 Amplitude')    
else
    figure  
    set(gcf,'Position',[50 50 1200 800])
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243;
    0.5 0.5 0.5];

    for ii = 1:3
        subplot(2,3,ii)
        dataAll  = [DS.number{ii}(:,1);DS.number{ii}(:,2);DS.number{ii}(:,3)];
        dataID = [ones(size(DS.number{ii},1),1)*1;ones(size(DS.number{ii},1),1)*2;ones(size(DS.number{ii},1),1)*3];
        stats.DSrate = groupStats(dataAll,dataID,'inAxis',true,'repeatedMeasures',true,'plotType','boxLinesSEM','color',colMat);    
        title('Dentate Spike 2 Rate')
        
        subplot(2,3,ii+3)
        dataAll  = [DS.amplitude{ii,1};DS.amplitude{ii,2};DS.amplitude{ii,3}];
        dataID = [ones(size(DS.amplitude{ii,1},1),1)*1;ones(size(DS.amplitude{ii,2},1),1)*2;ones(size(DS.amplitude{ii,3},1),1)*3];
        stats.DSamplitude = groupStats(dataAll,dataID,'inAxis',true,'color',colMat);    
        title('Dentate Spike 2 Amplitude')        
    end
    
end

saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\compiledDSpikes',tag,'.png'));
saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\compiledDSpikes',tag,'.fig'));
saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\compiledDSpikes',tag,'.eps'),'epsc');
save(strcat(parentDir,'Compiled\DentateSpikes\compiledDSpikes',tag,'.mat'),'stats')
end