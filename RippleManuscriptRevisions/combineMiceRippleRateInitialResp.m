function combineMiceRippleRateInitialResp

tag = 'mEC';% mEC, CA3, Bilateral mEC
if strcmp(tag,'mEC')==1    
    analogCh = 2;
else
    analogCh = 1;
end
 
if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ21\Final',...
         'IZ18\Final','IZ20\Final','IZ24\Final','IZ25\Final','IZ26\Final',...
         'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
end
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);
 
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

kk = 1; tt =1;

unitPSTH = [];
distRipPulses = [];
region = [];
putativeClass = [];

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m},'\'));  
    allSess = dir('*_sess*');
    
    for ii = 1:size(allSess,1)
%         if m==1 && ii==4
%             continue
%         end
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

        if exist([sessionInfo.FileName '.ripples.events.mat'],'file') 
           load([sessionInfo.FileName '.ripples.events.mat']);
        else
           disp('First calculate .ripple file! Skipping');
           return
        end
        %Load ripple LFP
        pyrCh = ripples.detectorinfo.detectionchannel;
%         lfp = bz_GetLFP(pyrCh);

        %% Load cell metrics
        if exist([sessionInfo.FileName '.cell_metrics.cellinfo.mat'],'file') 
            load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        else
            continue
        end          
            
        if exist([sessionInfo.FileName '.SleepState.states.mat'],'file') 
           load([sessionInfo.FileName '.SleepState.states.mat']);
        else
           disp('No sleep states associated with this session');
        end

        % Load pulses
        load([sessionInfo.FileName '.pulses.events.mat']);
        
        % Load spikes
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);        

        if m< 12
            pulTr = (pulses.stimComb==analogCh);
        else
            pulTr = (pulses.stimPerID'==1 & (pulses.stimComb==2 | pulses.stimComb==3)); % for the CA3 mice, both 2 and 3 are the same
        end
   
        % Only take pulses that happened during NREM sleep
        events = pulses.intsPeriods(1,pulTr)';         
        eventsUD = events(InIntervals(events,SleepState.ints.NREMstate));
        
        events = pulses.intsPeriods(2,pulTr)';  
        eventsDU  = events(InIntervals(events,SleepState.ints.NREMstate));   
       
        if ~isempty(eventsUD) && ~isempty(eventsDU)
%             [wavAvg, lfpAvg] = bz_eventWavelet(lfp, eventsUD, 'twin',[2 2],'plotWave',false,'plotLFP',false); 
%             freqs = wavAvg.freqs;
%             timest = wavAvg.timestamps;            
%             wavAvgeventUD(:,:,kk) = wavAvg.data;
%             lfpAvgeventUD(:,:,kk) = lfpAvg.data;
% 
%             [wavAvg, lfpAvg] = bz_eventWavelet(lfp, eventsDU, 'twin',[2 2],'plotWave',false,'plotLFP',false); 
%             wavAvgeventDU(:,:,kk) = wavAvg.data;
%             lfpAvgeventDU(:,:,kk) = lfpAvg.data;        
%             kk = kk+1;
            
            [stccg, tPul] = CCG({ripples.peaks eventsUD},[],'binSize',0.1,'duration',20,'norm','rate');
            distRipPulses = [distRipPulses; stccg(:,2,1)'];

             MUA = [];
             numMUA = 0;
             for jj = 1:size(spikes.UID,2)
                if strcmp(cell_metrics.brainRegion{jj},'CA1') == 1 && strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell') == 1
                    MUA = [MUA; spikes.times{jj}];
                    numMUA = numMUA+1;
                end
             end
            
            if ~isempty(MUA)
                [stccg, psthtimes] = CCG({MUA eventsUD},[],'binSize',0.05,'duration',20,'norm','rate');
                temp = stccg(:,2,1)./numMUA;
                unitPSTH = [unitPSTH; temp'];   
            end
%                 % Get psth
%                 [stccg, t] = CCG({spikes.times{jj} eventsUD},[],'binSize',0.1,'duration',20,'norm','rate');
%                 psthtimes = t;
%                 unitPSTH = [unitPSTH; stccg(:,2,1)'];
%                 % Assign numerical tag to putative class
%                 if strcmp(cell_metrics.putativeCellType{jj},'Narrow Interneuron') == 1
%                     putativeClass = [putativeClass; 1];
%                 elseif strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell') == 1
%                     putativeClass = [putativeClass; 2];
%                 else 
%                     putativeClass = [putativeClass; 3];                       
%                 end                    
% 
%                 % Assign numerical tag to region
%                 if strcmp(cell_metrics.brainRegion{jj},'CA1') == 1
%                     region = [region; 1];
%                 elseif strcmp(cell_metrics.brainRegion{jj},'DG') == 1
%                     region = [region; 2];
%                 elseif strcmp(cell_metrics.brainRegion{jj},'CA3') == 1
%                     region = [region; 2];
%                 else 
%                     region = [region; 3];                       
%                 end       
%             end
        end
    end
end

figure
set(gcf,'Renderer','painters')
% subplot(3,2,1)
% contourf(timest*1000,freqs,nanmean(wavAvgeventUD,3)',30,'LineColor','none');hold on;
% set(gca,'YScale','log');
% ylim([freqs(1) freqs(end)]);
% colormap jet;
% xlabel('time (ms)'); ylabel('frequency (Hz)');
% fmin =  freqs(round(length(freqs)/4)); 
% fmax =  fmin*2;
% lfp_avg = nanmean(lfpAvgeventUD,3);
% lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
% plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;
% 
% 
% subplot(3,2,2)
% contourf(timest*1000,freqs,nanmean(wavAvgeventDU,3)',30,'LineColor','none');hold on;
% set(gca,'YScale','log');
% ylim([freqs(1) freqs(end)]);
% colormap jet;
% xlabel('time (ms)'); ylabel('frequency (Hz)');
% fmin =  freqs(round(length(freqs)/4)); 
% fmax =  fmin*2;
% lfp_avg = nanmean(lfpAvgeventDU,3);
% lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
% plot(timest*1000,lfpSc,'w','LineWidth',2);hold on;

% subplot(3,2,3)
% imagesc(tPul,1:1:size(distRipPulses,1),zscore(distRipPulses,[],2))
% colormap(YlGnBu)
% colorbar    
% caxis([-2 2])

subplot(3,2,3)
imagesc(tPul,1:1:size(distRipPulses,1),distRipPulses)
colormap(YlGnBu)
colorbar    
caxis([-1 1])

subplot(3,2,5)
meanpsth = nanmean(distRipPulses,1);
stdpsth = nanstd(distRipPulses,1)./sqrt(size(distRipPulses,1));   
% meanpsth = nanmean(zscore(distRipPulses,[],2),1);
% stdpsth = nanstd(zscore(distRipPulses,[],2),1)./sqrt(size(distRipPulses,1));                  
hold on
fill([tPul; flipud(tPul)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(tPul,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
title('Ripple dist')
xlim([-2 7])

unitPSTH = unitPSTH./2;
cellIdx(1:size(unitPSTH,1)) = 1;
cellIdx = logical(cellIdx);

subplot(3,2,4)

imagesc(psthtimes,1:1:size(unitPSTH(cellIdx,:),1),unitPSTH(cellIdx,:))
colormap(YlGnBu)
colorbar    
caxis([-2 2])

subplot(3,2,6)
meanpsth = nanmean(unitPSTH(cellIdx,:),1);
stdpsth = nanstd(unitPSTH(cellIdx,:),1)./sqrt(size(unitPSTH(cellIdx,:),1));   
% meanpsth = nanmean(zscore(unitPSTH,[],2),1);
% stdpsth = nanstd(zscore(unitPSTH,[],2),1)./sqrt(size(unitPSTH,1));                  
hold on
fill([psthtimes; flipud(psthtimes)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
hi = line(psthtimes,meanpsth,'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
xlim([-2 7])
title('Ripple dist')

saveas(gcf,strcat(parentDir,'Compiled\Ripples\CortexCA1Comparison\PSTH2s',tag,'.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\CortexCA1Comparison\PSTH2s',tag,'.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Ripples\CortexCA1Comparison\PSTH2s',tag,'.fig'));
end
