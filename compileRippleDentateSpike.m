function compileRippleDentateSpike


mice =  {'IZ20\Final','IZ21\Final','IZ27\Saline',...
    'IZ28\Saline','IZ33\Saline'}; 
% {'IZ20\Final','IZ21\Final','IZ30\Final','IZ27\Saline',...
%     'IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline'}; 

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

DSpsthAll = [];
DSpsthPre = [];
DSpsthPost = [];

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
        load([sessionInfo.FileName '.ripples.events.mat']);
        
        
        pulTr = (pulses.stimComb==2);
        %Only select the pulses that happened in the home cage,
        %i.e., which were 5 seconds long
        homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
        homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
        pulTr = pulTr & homeCagePulseidx;
        events = pulses.intsPeriods(:,pulTr)';
            
        keepIdx(:,1) = InIntervals(ripples.peaks,events);    
        keepIdx(:,2) = InIntervals(ripples.peaks,(events-5));   
        
        %calculate psth of ds around ripples
        [stccg, t] = CCG({DS2.peaks ripples.peaks},[],'binSize',0.01,'duration',1);            
        psth = stccg(:,2,1)./sum(stccg(:,2,1));
        DSpsthAll = [DSpsthAll; psth'];
         
        if sum(keepIdx(:,1))>5
            [stccg, t] = CCG({DS2.peaks ripples.peaks(keepIdx(:,1)) },[],'binSize',0.01,'duration',1);            
            psth = stccg(:,2,1)./sum(stccg(:,2,1));
            DSpsthPost = [DSpsthPost; psth'];
        end
         
        if sum(keepIdx(:,2))>5
            [stccg, t] = CCG({DS2.peaks ripples.peaks(keepIdx(:,2))},[],'binSize',0.01,'duration',1);            
            psth = stccg(:,2,1)./sum(stccg(:,2,1));
            DSpsthPre = [DSpsthPre; psth'];
        end
         
         clear keepIdx
    end
end

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

subplot(2,3,1)
imagesc(t, 1:size(DSpsthAll,1), DSpsthAll)

subplot(2,3,4)
plot(t, nanmean(DSpsthAll))
xlim([-0.5 0.5])
line([0 0],[0 0.025])

subplot(2,3,2)
imagesc(t, 1:size(DSpsthPre,1), DSpsthPre)

subplot(2,3,5)
plot(t, nanmean(DSpsthPre))
xlim([-0.5 0.5])
line([0 0],[0 0.025])


subplot(2,3,3)
imagesc(t, 1:size(DSpsthPost,1), DSpsthPost)

subplot(2,3,6)
plot(t, nanmean(DSpsthPost))
xlim([-0.5 0.5])
line([0 0],[0 0.025])

saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\rippleDSpikes.png'));
saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\rippleDSpikes.fig'));
saveas(gca,strcat(parentDir,'Compiled\DentateSpikes\rippleDSpikes.eps'),'epsc');

end