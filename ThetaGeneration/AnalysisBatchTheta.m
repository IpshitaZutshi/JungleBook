
% summarizeScript
% Preprocessing and folder summary for expPipeline
% IZ-BuzsakiLab 2019
% 
% 1. LFP power profile
% 2. LFP peri-stimulus wavelet power
% 3. LFP theta CSD
% 4. Ripple detection
% 5. Ripple CSD
% 6. Theta-gamma modulation
% 7. Spike-waveform, autocorrelogram and spatial position 
% 8. Psth from analog-in inputs
% 9. Scoring behavioral task (alternation task)

%clear 

if ~exist('analogEv','var')
    analogEv = 64;
    numAnalog = 2;
    if ~isempty(analogEv)
        for ii = 1:numAnalog
            analogCh(ii) = (analogEv-1)+ii;
        end
    end
end

disp('Get components...');
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
% if isempty(dir('*.lfp'))
%     try 
%         bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
%     catch
%         disp('Problems with bz_LFPfromDat, resampling...');
%         ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
%             sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
%     end
% end


% badChannels = [64 65 66 67 68];
%  try SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels,'noPrompts',true); % try to sleep score
%  catch err
%      disp(err.message)
%  end
if ~exist('summ','dir')
    mkdir('summ'); % create folder
end
close all
nf = 1;
if ~exist('analysisList')
    %analysisList = [0 1 1 1 1 0 1 1 0];
    analysisList = [0 0 0 0 0 0 1 0 0];
end

% C = strsplit(pwd,'\');
% saveDir = strcat('F:\Dropbox\Project Data\mEC and CA1(old - cell death)\',C{3},'\',C{end},'\');%'Summ';
% if ~exist(saveDir,'dir')
%     mkdir(saveDir); % create folder
% end
stimloc = {'CA1','mEC','Both'};

%% 1. LFP Power spectrum
if analysisList(1)==1
   
    % Calculate power profile code
    try 
        disp('Theta modulation...');
        % Theta profile
        powerProfile_theta = bz_PowerSpectrumProfile_IZ([6 12],'channels',[0:(analogEv-1)],'showfig',true);
        % Gamma profile
        powerProfile_sg = bz_PowerSpectrumProfile_IZ([30 60],'channels',[0:(analogEv-1)],'showfig',true);
        % HFO profile
        powerProfile_hfo = bz_PowerSpectrumProfile_IZ([120 250],'channels',[0:(analogEv-1)],'showfig',true);
        % get channels of interest % max theta power above pyr layer
        [~, a] = min(abs(powerProfile_hfo.mean - powerProfile_theta.mean));
        region.CA1sp = powerProfile_hfo.channels(a);

        for ii = 1:size(sessionInfo.AnatGrps,2)
            if ismember(region.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                p_th = powerProfile_theta.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                p_hfo = powerProfile_hfo.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                p_gs = powerProfile_sg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                channels = sessionInfo.AnatGrps(ii).Channels;
                [~,a] = max(p_th(1:find(channels == region.CA1sp)));
                region.CA1so = channels(a);
                [~,a] = max(p_th(find(channels == region.CA1sp):end));
                region.CA1slm = channels(a - 1+ find(channels == region.CA1sp));
                
                figure(nf)
                hold on
                plot(1:length(channels), zscore(p_th));
                plot(1:length(channels), zscore(p_hfo));
                plot(1:length(channels), zscore(p_gs));
                xlim([1 length(channels)]);
                set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
                ax = axis;
                xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
                plot(find(channels == region.CA1sp)*ones(2,1),ax([3 4]),'-k');
                plot(find(channels == region.CA1so)*ones(2,1),ax([3 4]),'--k');
                plot(find(channels == region.CA1slm)*ones(2,1),ax([3 4]),'-.k');
                legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');
                saveas(figure(nf),'Summ\regionDef.png');
                nf = nf + 1;
            end
        end
        save([sessionInfo.FileName,'.region.mat'],'region');

        % power profile of pyr channel of all session
        lfpT = bz_GetLFP(region.CA1sp,'noPrompts',true);
        params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
        [S,t,f] = mtspecgramc(single(lfpT.data),[2 1],params);
        S = log10(S); % in Db
        S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending

        figure(nf);
        subplot(3,3,[1 2 3 4 5 6])
        imagesc(t,f,S_det',[-1.5 1.5]);
        set(gca,'YDir','normal','XTick',[]); ylabel('Freqs');
        subplot(3,3,[7 8 9]);
        plot(f,mean(S,1));
        ylabel('Power'); xlabel('Freqs');
        saveas(figure(nf),'Summ\Spectogram.png');
        nf = nf + 1;

    catch err
        disp(err.message)
        warning('Error in LFP profile calculation');
    end

end

%% 2. LFP power before and after stim

if analysisList(2)==1
   % try
        %Calculate power and frequency profile before and after stim
        % First calculate wavelet
        [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
        channels = [1:sessionInfo.nChannels]-1;
        if exist([sessionInfo.FileName '.PowerSpectrumProfile_6_12.channelinfo.mat'],'file')  && ...
                exist([sessionInfo.FileName '.PowerSpectrumProfile_30_60.channelinfo.mat'],'file')
            powerProfile_theta = load([sessionInfo.FileName '.PowerSpectrumProfile_6_12.channelinfo.mat']);
            powerProfile_sg = load([sessionInfo.FileName '.PowerSpectrumProfile_30_60.channelinfo.mat']);
        else
            disp('First calculate the theta and gamma profile!')
            return
        end
        figure(nf)    
        for i = 1:(numAnalog+1)
            try
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
                if isempty(pulTr)
                    continue;
                end
                events = pulses.intsPeriods(1,pulTr);
                events = events(((events + 5) <= max(powerProfile_theta.powerProfile.time(1,:))) & ((events - 5) > 0));
                %Load theta and gamma power
                thetaProfile_Pre = [];
                thetaProfile_Post = [];
                gammaProfile_Pre = [];
                gammaProfile_Post = [];

                for pp = 1:length(events)
                    events_tmp = abs(powerProfile_theta.powerProfile.time(1,:)-events(pp));
                    [~,idx_temp] = min(events_tmp);
                    thetaProfile_Pre(:,pp) = powerProfile_theta.powerProfile.power(:,(idx_temp-1));
                    thetaProfile_Post(:,pp) = powerProfile_theta.powerProfile.power(:,(idx_temp+1));

                    gammaProfile_Pre(:,pp) = powerProfile_sg.powerProfile.power(:,(idx_temp-1));
                    gammaProfile_Post(:,pp) = powerProfile_sg.powerProfile.power(:,(idx_temp+1));
                end

                subplot((numAnalog+1),2,2*i-1)
                cmap = jet(size(sessionInfo.AnatGrps,2));

                for ii = 1:(size(sessionInfo.AnatGrps,2)-1)
                    [Lia] = ismember(sessionInfo.AnatGrps(ii).Channels, channels);
                    nC = length(sessionInfo.AnatGrps(ii).Channels):-1:1;
                    nC = nC(Lia)';
                    hold on
                    mean_pre = mean(thetaProfile_Pre(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),2);
                    std_pre = std(thetaProfile_Pre(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),0,2);
                    %fill([(mean_pre-std_pre) flip(mean_pre+std_pre)],[nC flip(nC)],cmap(ii,:),'FaceAlpha',.2,'EdgeColor','none')
                    if (strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ11')==1) && ii ==3
                        plot(mean_pre(1:(end-5)),nC((6:end)),'color','k','LineWidth',1);
                    else
                        plot(mean_pre,nC(Lia),'color','k','LineWidth',1);
                    end

                    hold on
                    mean_post = mean(thetaProfile_Post(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),2);
                    std_post = std(thetaProfile_Post(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),0,2);
                    %fill([(mean_post-std_post) flip(mean_post+std_post)],[nC flip(nC)],cmap(ii,:),'FaceAlpha',.2,'EdgeColor','none')
                    if (strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ11')==1) && ii ==3
                        plot(mean_post(1:(end-5)),nC(6:end),'color','b','LineWidth',1);
                    else
                        plot(mean_post,nC(Lia),'color','b','LineWidth',1);
                    end
                 %   xlim([15 50])
                end
                xlabel('Amplitude')
                ylabel('Channel Number')
                title(['Theta Power' stimloc{i}]);

                subplot((numAnalog+1),2,2*i)
                for ii = 1:(size(sessionInfo.AnatGrps,2)-1)
                    [Lia] = ismember(sessionInfo.AnatGrps(ii).Channels, channels);
                    nC = length(sessionInfo.AnatGrps(ii).Channels):-1:1;
                    nC = nC(Lia)';
                    hold on
                    mean_pre = mean(gammaProfile_Pre(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),2);
                    std_pre = std(gammaProfile_Pre(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),0,2);
                    if (strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ11')==1) && ii ==3
                        plot(mean_pre(1:(end-5)),nC((6:end)),'color','k','LineWidth',1);
                    else
                        plot(mean_pre,nC(Lia),'color','k','LineWidth',1);
                    end

                    hold on
                    mean_post = mean(gammaProfile_Post(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),2);
                    std_post = std(gammaProfile_Post(sessionInfo.AnatGrps(ii).Channels(Lia)+1,:),0,2);
                    if (strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ11')==1) && ii ==3
                        plot(mean_post(1:(end-5)),nC((6:end)),'color','b','LineWidth',1);
                    else
                        plot(mean_post,nC(Lia),'color','b','LineWidth',1);
                    end
                 %   xlim([15 40])
                end
                xlabel('Amplitude')
                ylabel('Channel Number')            
                title(['Gamma Power' stimloc{i}]);
            catch err
                disp(err.message)
                warning('Error in LFP wavelet calculation');
            end
        end
        %saveas(figure(nf),strcat(saveDir,'powerProfile_Perievent.png'));
        saveas(figure(nf),'Summ\powerProfile_Perievent.png');
        saveas(figure(nf),'Summ\powerProfile_Perievent.eps','epsc');
        saveas(figure(nf),'Summ\powerProfile_Perievent.fig');
        nf = nf + 1;
end


%% 3. LFP CSD
if analysisList(3)==1
    
    if exist([sessionInfo.FileName '.region.mat'],'file') 
        load([sessionInfo.FileName '.region.mat']);
        pyrCh = 7;%region.CA1sp;
        for ii = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ii).Channels)
                thetaCh = find(sessionInfo.AnatGrps(ii).Channels==pyrCh);
            end
        end
    else
        thetaCh = ceil(size(sessionInfo.AnatGrps(1).Channels,2)/2);
    end

    for i = 1:(numAnalog)%+1)
        try
            disp('LFP CSD...');
            disp('Psth and CSD from analog-in inputs...');
            [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);

            %Calculate CSD PSTH
            if exist('pulses') && ~isempty(pulses.intsPeriods(1,:)')
                for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
                    if (strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ11')==1) && jj ==3
                        lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:14),'noPrompts', true);
                    else
                        lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels,'noPrompts', true);
                        if ~ismember(thetaCh, sessionInfo.AnatGrps(jj).Channels)
                            thetaChTemp = ceil(size(sessionInfo.AnatGrps(jj).Channels,2)/2);
                        else
                            thetaChTemp = thetaCh;
                        end
                    end
                    twin = 0.2;
                    if i<=numAnalog
                        pulTr = (pulses.stimComb==i);
                    else
                        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                    end
                    if isempty(pulTr)
                        continue;
                    end
                    [csd,lfpAvg] = bz_eventThetaCSD(lfp,pulses.intsPeriods(1,pulTr)','twin',[twin twin],'plotLFP',false,'plotCSD',false,'thetaCh',thetaChTemp);
                    taxis = linspace(-twin,twin,size(csd.data_pre,1));
                    cmax = max(max(csd.data_pre)); 
                    figure(nf);
                    %set(gcf,'Position',[100 100 2500 1200])
                    subplot((size(sessionInfo.AnatGrps,2)-1),3,3*(jj-1)+1);
                    contourf(taxis,1:size(csd.data_pre,2),csd.data_pre',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('Baseline #',num2str(jj),' # ',num2str(csd.num(1))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data_pre,2)
                        plot(taxis,(lfpAvg.data_pre(:,kk)/1000)+kk-1,'k')
                    end
                    
                    %cmax = max(max(csd.data_post)); 
                    subplot((size(sessionInfo.AnatGrps,2)-1),3,3*(jj-1)+2);
                    contourf(taxis,1:size(csd.data_post,2),csd.data_post',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('Stim #',num2str(jj),' # ',num2str(csd.num(2))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data_post,2)
                        plot(taxis,(lfpAvg.data_post(:,kk)/1000)+kk-1,'k')
                    end
                    
                    subplot((size(sessionInfo.AnatGrps,2)-1),3,3*(jj-1)+3);
                    contourf(taxis,1:size(csd.data_post,2),(csd.data_pre'-csd.data_post'),40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('#',num2str(jj),'Baseline-Stim'),'FontWeight','normal'); 
                    caxis([-100 100])
                    colorbar
                end
                %saveas(figure(nf),strcat(saveDir,'LFP_CSD_',stimloc{i},'.png'));
                saveas(figure(nf),strcat('Summ\LFP_CSD_Ch',num2str((analogEv-1)+i),'.png'));
                saveas(figure(nf),strcat('Summ\LFP_CSD_Ch',num2str((analogEv-1)+i),'.eps'),'epsc');
                saveas(figure(nf),strcat('Summ\LFP_CSD_Ch',num2str((analogEv-1)+i),'.fig'));
                nf = nf + 1;
            end
        catch err
            disp(err.message)
            warning('Error in LFP amplitude and CSD calculation');
        end
    end
end

%% 4. Ripple detection
if analysisList(4)==1
    try
        if exist([sessionInfo.FileName '.region.mat'],'file') 
            load([sessionInfo.FileName '.region.mat']);
            pyrCh = region.CA1sp;
        else
            disp('First calculate .region file to identify ripple channel! Skipping');
            return
        end
        if strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ23')==1
            noiseCh = 29;
        else
            noiseCh = 15;
        end
        [ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true);
        passband = [100 200];
        
        lfp = bz_GetLFP(pyrCh,'noPrompts', true);
        signal = bz_Filter(double(lfp.data),'filter','butter','passband',passband,'order', 3);
        timestamps = lfp.timestamps;
        [maps,data] = bz_RippleStats(signal,timestamps,ripples);
        
        disp('Getting analog-in inputs...');
        [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
        
        figure(nf)
        
        if exist('pulses') && ~isempty(pulses.intsPeriods(1,:)')
            for i = 1:(numAnalog+1)
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
                if isempty(pulTr)
                    continue;
                end
                events = pulses.intsPeriods(1,pulTr);
                events = events(((events + 5) <= max(ripples.peaks)) & ((events - 5) > 0));
                
                %Generate logicals for ripples in pre versus post
                ripple_pre(1:length(ripples.peaks)) = 0;
                ripple_prestim(1:length(ripples.peaks)) = 0;
                ripple_post(1:length(ripples.peaks)) = 0;
                ripple_poststim(1:length(ripples.peaks)) = 0;
                
                for pp = 1:length(ripples.peaks)
                    tempDiff = ripples.peaks(pp) - events;
                    
                    if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                       [~,idxmin] =  min(abs(tempDiff));
                       if tempDiff(idxmin) > 0
                           ripple_post(pp) = 1;
                       elseif tempDiff(idxmin) < 0
                           ripple_pre(pp) = 1;
                       end
                    elseif min(abs(tempDiff)) >=5 && min(abs(tempDiff)) < 10
                        [~,idxmin] =  min(abs(tempDiff));
                        if tempDiff(idxmin) > 0
                            ripple_poststim(pp) = 1;
                       elseif tempDiff(idxmin) < 0
                           ripple_prestim(pp) = 1;
                       end
                    else
                        continue
                    end
                end
                
                ripple_pre = logical(ripple_pre);
                ripple_post = logical(ripple_post);
                ripple_poststim = logical(ripple_poststim);
                ripple_prestim = logical(ripple_prestim);
                
                % Plot the results
                subplot((numAnalog+1),4,4*(i-1)+1)  
                bar([sum(ripple_prestim) sum(ripple_pre) sum(ripple_post) sum(ripple_poststim)]);
                ylabel('Number of ripples')
                
                subplot((numAnalog+1),4,4*(i-1)+2)
                title('Ripple power')
                if sum(ripple_pre>0)
                    histogram(data.peakAmplitude(ripple_pre),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                hold on
                if sum(ripple_post>0)
                    histogram(data.peakAmplitude(ripple_post),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                if sum(ripple_poststim>0)
                    histogram(data.peakAmplitude(ripple_poststim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                if sum(ripple_prestim>0)
                    histogram(data.peakAmplitude(ripple_prestim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                legend('Pre','Post','PostStim','PreStim','Location','best')
                ylabel('Proportion')
                xlabel('Peak ripple power')
                
                subplot((numAnalog+1),4,4*(i-1)+3)
                title('Ripple duration')
                if sum(ripple_pre>0)
                    histogram(data.duration(ripple_pre),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                hold on
                if sum(ripple_post>0)
                    histogram(data.duration(ripple_post),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                if sum(ripple_poststim>0)
                    histogram(data.duration(ripple_poststim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end         
                if sum(ripple_prestim>0)
                    histogram(data.duration(ripple_prestim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end                   
                legend('Pre','Post','PostStim','PreStim','Location','best')
                ylabel('Proportion')
                xlabel('Ripple duration')
                title([stimloc{i} 'Silencing']);
                
                subplot((numAnalog+1),4,4*(i-1)+4)
                title('Ripple frequency')
                if sum(ripple_pre>0)
                    histogram(data.peakFrequency(ripple_pre),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                hold on
                if sum(ripple_post>0)
                    histogram(data.peakFrequency(ripple_post),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end
                if sum(ripple_poststim>0)
                    histogram(data.peakFrequency(ripple_poststim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end   
                if sum(ripple_prestim>0)
                    histogram(data.peakFrequency(ripple_prestim),25,'EdgeAlpha',0,'FaceAlpha',0.5);
                end                 
                legend('Pre','Post','PostStim','PreStim','Location','best')
                ylabel('Proportion')
                xlabel('Ripple frequency')
                title([stimloc{i} 'Silencing']);
                
            end
            %saveas(figure(nf),strcat(saveDir,'Ripples_perievent.png'));
            saveas(figure(nf),'Summ\Ripples_perievent.png');
            saveas(figure(nf),'Summ\Ripples_perievent.eps');
            saveas(figure(nf),'Summ\Ripples_perievent.fig');
            nf = nf + 1;
        end
     
    catch err
        disp(err.message)
        warning('Error in ripple analysis');
    end
end


%% 5. Ripple CSD
if analysisList(5)==1
    try
        if exist([sessionInfo.FileName '.region.mat'],'file') 
            load([sessionInfo.FileName '.hippocampalLayers.channelInfo.mat']);
            pyrCh = hippocampalLayers.pyramidal;
        else
            disp('First calculate .region file to identify ripple channel! Skipping');
            return
        end
        if strcmp(sessionInfo.FileName(1:4),'IZ11')==1 || strcmp(sessionInfo.FileName(1:4),'IZ15')==1 || strcmp(sessionInfo.FileName(1:4),'IZ23')==1
            noiseCh = 29;
        else
            noiseCh = [];
        end
        [ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true);

        disp('Getting analog-in inputs...');
        [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
  
        if exist('pulses') && ~isempty(pulses.intsPeriods(1,:)')
            
            figure(nf)
            for i = 1:(numAnalog+1)
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
                if isempty(pulTr)
                    continue;
                end
                events = pulses.intsPeriods(1,pulTr);
                events = events(((events + 5) <= max(ripples.peaks)) & ((events - 5) > 0));
                
                %Generate logicals for ripples in pre versus post
                ripple_pre(1:length(ripples.peaks)) = 0;
                ripple_post(1:length(ripples.peaks)) = 0;
                
                for pp = 1:length(ripples.peaks)
                    tempDiff = ripples.peaks(pp) - events;
                    
                    if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                       [~,idxmin] =  min(abs(tempDiff));
                       if tempDiff(idxmin) > 0
                           ripple_post(pp) = 1;
                       elseif tempDiff(idxmin) < 0
                           ripple_pre(pp) = 1;
                       end
                    else
                        continue
                    end
                end
                
                ripple_pre = logical(ripple_pre);
                ripple_post = logical(ripple_post);
                for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
                    lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:numel(sessionInfo.AnatGrps(jj).Channels)-mod(numel(sessionInfo.AnatGrps(jj).Channels),8)),'noPrompts', true);
                    twin = 0.2;
                    [csd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks(ripple_pre),'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    figure(nf);
                    subplot(2,(size(sessionInfo.AnatGrps,2)-1),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('CTRL, #',num2str(sum(ripple_pre))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                    
                    [csd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks(ripple_post),'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    %cmax = max(max(csd.data)); 
                    figure(nf);
                    subplot(2,(size(sessionInfo.AnatGrps,2)-1),jj+(size(sessionInfo.AnatGrps,2)-1));
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('STIM, #',num2str(sum(ripple_post))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                %saveas(figure(nf),strcat(saveDir,'RipplesCSD_',stimloc{i},'.png'));
                saveas(figure(nf),strcat('Summ\RipplesCSD_Ch',num2str((analogEv-1)+i),'.png'));
                saveas(figure(nf),strcat('Summ\RipplesCSD_Ch',num2str((analogEv-1)+i),'.eps'));
                saveas(figure(nf),strcat('Summ\RipplesCSD_Ch',num2str((analogEv-1)+i),'.fig'));
                nf = nf + 1;
            end
        end
     
    catch err
        disp(err.message)
        warning('Error in ripple analysis');
    end
end

%% 6. Theta-gamma comodulation
if analysisList(6)==1
    try
        getPeriStimModIdx;
        nf = nf+1;
    catch
        disp(err.message)
        warning('Error on theta-gamma coherence! ');
    end
end
%% Load spikes
% try
%     spikes = bz_LoadPhy('noPrompts',true);
% catch
% end

%% 7. Spike-waveform, autocorrelogram and spatial position 

if analysisList(7)==1
   try 
    % plot spikes summary
        %if ~exist('spikes','var')
            spikes = bz_LoadPhy('noPrompts',true);
        %end
    
        disp('Plotting spikes summary...');
        load('chanMap.mat','xcoords','ycoords');
        xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
        ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
        for jj = 1:size(spikes.UID,2)
            figure(nf)
            dur = 0.08;
            set(gcf,'Position',[100 100 2500 1200])
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            ccg=CCG(spikes.times,[],'binSize',0.001,'duration',0.08);
            xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
            area(xt,ccg(:,jj,jj),'LineStyle','none');

            xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,size(spikes.filtWaveform{jj},1)); % waveform
            spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
            hold on
            plot(xt2,spk)

            plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
            plot((xcoords(spikes.maxWaveformCh(jj)+1)*30) + dur/2*1000+5+60, ycoords(spikes.maxWaveformCh(jj)+1)*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('Counts/ norm');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (100 ms /1.5ms)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        %saveas(figure(nf),strcat(saveDir,'spikes.png'));
        saveas(figure(nf),'Summ\spikes.png');
        saveas(figure(nf),'Summ\spikes.eps');
        saveas(figure(nf),'Summ\spikes.fig');
        nf = nf + 1;
        
        win = [-0.3 0.3];
        disp('Plotting CCG...');
        figure(nf);
        set(gcf,'Position',[100 100 2500 1200])
        [allCcg, t] = CCG(spikes.times,[],'binSize',0.005,'duration',0.6);
        indCell = [1:size(allCcg,2)];
        for jj = 1:size(spikes.UID,2)
            subplot(7,ceil(size(spikes.UID,2)/7),jj);
            cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
            imagesc(t,1:max(indCell)-1,cc)
            set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
            hold on
            zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
            zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
            plot(t, zmean,'k','LineWidth',2);
            xlim([win(1) win(2)]); ylim([0 max(indCell)-1]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('Cell');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(figure(nf),'Summ\CrossCorr.png'); 
        saveas(figure(nf),'Summ\CrossCorr.eps'); 
        nf = nf + 1;
    catch err
        disp(err.message)
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end

%% 8. Spike PSTH from analog-in inputs
if analysisList(8)==1
    for i = 1:(numAnalog)%+1)
        try 
            
            disp('Psth and CSD from analog-in inputs...');
            [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);
            % copy NPY file to Kilosort folder, for phy psth option
%             writeNPY(pulses.intsPeriods(1,:)','pulTime.npy'); % 
%             kilosort_path = dir('*Kilosort*');
%             copyfile('pulTime.npy', strcat(kilosort_path.name,'\','pulTime.npy')); % copy pulTime to kilosort folder

            if exist('pulses') && ~isempty(pulses.intsPeriods(1,:)')
                % CSD
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
                if isempty(pulTr)
                    continue;
                end
                for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
                    lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:numel(sessionInfo.AnatGrps(jj).Channels)-mod(numel(sessionInfo.AnatGrps(jj).Channels),8)),'noPrompts', true);
                    twin = 0.5;
                    [csd,lfpAvg] = bz_eventCSD(lfp,pulses.intsPeriods(1,pulTr)','twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    figure(nf);
                    set(gcf,'Position',[100 100 2500 1200])
                    subplot(1,(size(sessionInfo.AnatGrps,2)-1),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('STIMULATION, Shank #',num2str(jj)),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k','LineWidth',2)
                    end
                end
                %saveas(figure(nf),strcat(saveDir,'stimulationCSD_',stimloc{i},'.png'));
                saveas(figure(nf),strcat('Summ\stimulationCSDVer2_Ch',num2str((analogEv-1)+i),'.png'));
                saveas(figure(nf),strcat('Summ\stimulationCSDVer2_Ch',num2str((analogEv-1)+i),'.eps'));
                saveas(figure(nf),strcat('Summ\stimulationCSDVer2_Ch',num2str((analogEv-1)+i),'.fig'));
                nf = nf + 1;

%                 % PSTH
                if ~exist('spikes','var')
                    spikes = bz_LoadPhy('noPrompts',true);
                end

                st = pulses.intsPeriods(1,pulTr)';
                win = [-1 1];
                disp('Plotting spikes raster and psth...');
                figure(nf);
                set(gcf,'Position',[100 100 2500 1200])
                for jj = 1:size(spikes.UID,2)
                    rast_x = []; rast_y = [];
                    for kk = 1:length(st)
                        temp_rast = spikes.times{jj} - st(kk);
                        temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                        rast_x = [rast_x temp_rast'];
                        rast_y = [rast_y kk*ones(size(temp_rast))'];
                    end
                    [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.1,'duration',12);
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    plot(rast_x, rast_y,'.','MarkerSize',2)
                    hold on
                    plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',1);
                    xlim([win(1) win(2)]); ylim([0 kk]);
                    line([0 0],[0 kk])
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('Trial');
                    elseif jj == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                %saveas(figure(nf),strcat(saveDir,'SpikePSTH_',stimloc{i},'.png'));
                saveas(figure(nf),strcat('Summ\SpikePSTH_Ch',num2str((analogEv-1)+i),'.png'));
                saveas(figure(nf),strcat('Summ\SpikePSTH_Ch',num2str((analogEv-1)+i),'.eps')); 
                saveas(figure(nf),strcat('Summ\SpikePSTH_Ch',num2str((analogEv-1)+i),'.fig')); 
                nf = nf + 1;
             end  
        catch err
            disp(err.message)
            warning('Error on Psth and CSD from analog-in inputs! ');
        end
        clear temp_rast rast_x rast_y st
    end
end
%% 9. Behavior, alternation
if analysisList(9)==1
    try 
        getSessionTracking('convFact',0.1149);
        getSessionArmChoice('task','alternation','force',true);
        getSessionLinearize('forceReload',true);
        getSessionPulses('force',true);
%         
%         position{1} = [behavior.timestamps(behavior.masks.recording==1)...
%             behavior.position.x(behavior.masks.recording==1)...
%             behavior.position.y(behavior.masks.recording==1)];
%         spikes = bz_LoadPhy;
%         fm = bz_firingMapAvg(position,spikes,'nBins',120,'minTime',0.05);
%         
%         figure
%         for ii = 1:size(fm.rateMaps,1)
%             imagesc(fm.rateMaps{ii}{1})
%             pause;
%         end
    catch
        warning('It has not been possible to run the behaviour code...');
    end
end

% 
% %% 5. THETA AND GAMMA PROFILE/ MODULATION
% if analisysList(7)
%     try 
%         disp('Theta modulation...');
%         % Theta profile
%         powerProfile_theta = bz_PowerSpectrumProfile([6 12],'channels',[0:63],'showfig',true);
%         saveas(figure(nf),'Summ\powerProfile_6_12.png');
%         nf = nf + 1;
%         % Gamma profile
%         powerProfile_sg = bz_PowerSpectrumProfile([30 60],'channels',[0:63],'showfig',true);
%         saveas(figure(nf),'Summ\powerProfile_30_60.png');
%         nf = nf + 1;
%         % HFO profile
%         powerProfile_hfo = bz_PowerSpectrumProfile([120 250],'channels',[0:63],'showfig',true);
%         saveas(figure(nf),'Summ\powerProfile_120_250.png');
%         nf = nf + 1;
%         % get channels of interest % max theta power above pyr layer
%         [~, a] = max(powerProfile_hfo.mean);
%         region.CA1sp = powerProfile_hfo.channels(a);
% 
%         for ii = 1:size(sessionInfo.AnatGrps,2)
%             if ismember(region.CA1sp, sessionInfo.AnatGrps(ii).Channels)
%                 p_th = powerProfile_theta.mean(sessionInfo.AnatGrps(ii).Channels + 1);
%                 p_hfo = powerProfile_hfo.mean(sessionInfo.AnatGrps(ii).Channels + 1);
%                 p_gs = powerProfile_sg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
%                 channels = sessionInfo.AnatGrps(ii).Channels;
%                 [~,a] = max(p_th(1:find(channels == region.CA1sp)));
%                 region.CA1so = channels(a);
%                 [~,a] = max(p_th(find(channels == region.CA1sp):end));
%                 region.CA1slm = channels(a - 1+ find(channels == region.CA1sp));
%                 
%                 figure
%                 hold on
%                 plot(1:length(channels), zscore(p_th));
%                 plot(1:length(channels), zscore(p_hfo));
%                 plot(1:length(channels), zscore(p_gs));
%                 xlim([1 length(channels)]);
%                 set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
%                 ax = axis;
%                 xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
%                 plot(find(channels == region.CA1sp)*ones(2,1),ax([3 4]),'-k');
%                 plot(find(channels == region.CA1so)*ones(2,1),ax([3 4]),'--k');
%                 plot(find(channels == region.CA1slm)*ones(2,1),ax([3 4]),'-.k');
%                 legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');
%                 saveas(figure(nf),'Summ\regionDef.png');
%                 nf = nf + 1;
%             end
%         end
%         save([sessionInfo.FileName,'.region.mat'],'region');
% 
%         % power profile of pyr channel of all session
%         lfpT = bz_GetLFP(region.CA1so,'noPrompts',true);
%         params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
%         [S,t,f] = mtspecgramc(single(lfpT.data),[2 1],params);
%         S = log10(S); % in Db
%         S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending
% 
%         figure(nf);
%         subplot(3,4,[1 2 3 5 6 7])
%         imagesc(t,f,S_det',[-1.5 1.5]);
%         set(gca,'XTick',[]); ylabel('Freqs');
%         subplot(3,4,[4 8]);
%         plot(mean(S,1),f);
%         set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
%         subplot(3,4,[9 10 11]);
%         % to do
%         xlabel('s'); ylabel('mov (cm/s)');
%         saveas(figure(nf),'Summ\spectrogramAllSession.png');
%         nf = nf + 1;
%         % phase modulation by spikes
%         spikes = bz_LoadPhy('noPrompts',true);
%         for ii = 1:size(spikes.filtWaveform,2)
%             try spikes.region{ii} = sessionInfo.region(spikes.maxWaveformCh(ii));
%             catch
%                 spikes.region{ii} = 'NA';
%             end
%         end
%         save([sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes');
%         
%         % thetaMod modulation
%         PLD = bz_PhaseModulation(spikes,lfpT,[4 12],'plotting',false,'powerThresh',1);  
%         disp('Theta modulation...');
%         figure(nf)
%         set(gcf,'Position',[100 100 2500 1200]);
%         for jj = 1:size(spikes.UID,2)
%             subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
%             area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
%             hold on
%             ax = axis;
%             x = 0:.001:4*pi;
%             y = cos(x);
%             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%             xlim([0 4*pi]);
%             title(num2str(jj),'FontWeight','normal','FontSize',10);
% 
%             if jj == 1
%                 ylabel('prob');
%             elseif jj == size(spikes.UID,2)
%                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                 xlabel('phase (rad)');
%             else
%                 set(gca,'YTick',[],'XTick',[]);
%             end
%         end
%         saveas(figure(nf),'Summ\thetaMod.png');
%         nf = nf + 1;
%         
%         % gamma modulation
%         PLD = bz_PhaseModulation(spikes,lfpT,[30 60],'plotting',false,'powerThresh',2);     
%         disp('Gamma modulation...');
%         figure(nf)
%         set(gcf,'Position',[100 100 2500 1200]);
%         for jj = 1:size(spikes.UID,2)
%             subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
%             area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
%             hold on
%             ax = axis;
%             x = 0:.001:4*pi;
%             y = cos(x);
%             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%             xlim([0 4*pi]);
%             title(num2str(jj),'FontWeight','normal','FontSize',10);
% 
%             if jj == 1
%                 ylabel('prob');
%             elseif jj == size(spikes.UID,2)
%                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                 xlabel('phase (rad)');
%             else
%                 set(gca,'YTick',[],'XTick',[]);
%             end
%         end
%         saveas(figure(nf),'Summ\gammaMod.png');
%         nf = nf + 1;  
%     catch
%         warning('It has not been possible to run theta and gamma mod code...');
%     end
% end
% 
