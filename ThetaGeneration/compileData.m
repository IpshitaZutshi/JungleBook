% Compile data across all sessions

function compiledData = compileData(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isdir);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'classify',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
classify = p.Results.classify;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');


if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

compiledData.FR = [];
compiledData.TroughtoPeak = [];
compiledData.putativeClass = [];
compiledData.classifier = [];
compiledData.psthtimes = [];
compiledData.psthstccg = [];
compiledData.psthAnalogEv = [];
compiledData.region = [];

%% Start collecting data
for ii = 1:size(allSess,1)

        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
        % Load spikes
       if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
            spikes = bz_LoadPhy('noPrompts',true);
       else
           continue
       end

        % Load pulses

        sess = bz_getSessionInfo(pwd,'noPrompts',true);
        if exist([sess.FileName '.pulses.events.mat'],'file') 
            load([sess.FileName '.pulses.events.mat']);
        else
            continue
        end

        % calculate cell metrics for classification
        [cell_metrics] = bz_CellMetricsSimpleIZ(spikes,sessionInfo.rates.wideband);
        compiledData.FR = [compiledData.FR; cell_metrics.FR]; 
        compiledData.TroughtoPeak = [compiledData.TroughtoPeak; cell_metrics.filtWaveform.TroughtoPeak];
        compiledData.putativeClass = [compiledData.putativeClass; cell_metrics.putativeClass];

        for i = 1:(numAnalog+1)
            if exist('pulses')
                if i<=numAnalog
                    pulTr = (pulses.stimComb==i);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
                end
            end

            st = pulses.intsPeriods(1,pulTr)';
            for jj = 1:size(spikes.UID,2)
                % Get psth
                [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.1,'duration',12);
                compiledData.psthtimes = [compiledData.psthtimes; t'];
                compiledData.psthstccg = [compiledData.psthstccg; stccg(:,2,1)'];
                compiledData.psthAnalogEv = [compiledData.psthAnalogEv; i];
                compiledData.classifier = [compiledData.classifier;cell_metrics.putativeClass(jj)]; 

                if ~classify
                    compiledData.region = [compiledData.region; 1];
                else
                    maxChan = spikes.maxWaveformCh(jj);
                    idxChan = [];
                    for aa = 1:(size(sessionInfo.AnatGrps,2)-1)
                        if isempty(idxChan)
                            idxChan = find(sessionInfo.AnatGrps(aa).Channels==maxChan);
                        end
                    end        
                    if ~isempty(idxChan)
                        if idxChan <= 34
                            compiledData.region = [compiledData.region; 1];
                        else
                            compiledData.region = [compiledData.region; 2];
                        end
                    else
                        compiledData.region = [compiledData.region; 1];
                        continue;
                    end
                end
            end
 %       end
        end  
end


figure
scatter(compiledData.FR, compiledData.TroughtoPeak,'ro')
hold on
scatter(compiledData.FR(compiledData.putativeClass==2), compiledData.TroughtoPeak(compiledData.putativeClass==2),'bo')
scatter(compiledData.FR(compiledData.putativeClass==1), compiledData.TroughtoPeak(compiledData.putativeClass==1),'ko')
xlabel('Firing Rate')
ylabel('Trough to Peak Ratio')

CellClass = {'Pyr','IN'};
col = [0.9856 0.7372 0.2537;0.1986 0.7214 0.6310];%;0.4733 0.3374 0.2798];
%col = [0.2305 0.2510 0.6173;0.1986 0.7214 0.6310];
colormap default

for cc = 1:size(unique(compiledData.region),1)
    figure(cc)
    set(gcf,'Position',[100 100 1000 600])
    for class = 1:2
        psth_ca1 = compiledData.psthstccg(compiledData.psthAnalogEv==1 & compiledData.region==cc & compiledData.classifier==class,:);
        psth_mec = compiledData.psthstccg(compiledData.psthAnalogEv==2 & compiledData.region==cc & compiledData.classifier==class,:);
        psth_ca1_mec = compiledData.psthstccg(compiledData.psthAnalogEv==3 & compiledData.region==cc & compiledData.classifier==class,:);

        kk = 1;
        for ii = 1:size(psth_ca1,1)
            FR = max([nanmean(psth_ca1(ii,1:60)), nanmean(psth_mec(ii,1:60)),nanmean(psth_ca1_mec(ii,1:60))]);
            if FR <0.2
                continue
            else
                FR = max(psth_ca1(ii,:));
                psth_norm_ca1(kk,:) = psth_ca1(ii,:)/FR;

                FR = max(psth_mec(ii,:));
                psth_norm_mec(kk,:) = psth_mec(ii,:)/FR;

                FR = max(psth_ca1_mec(ii,:));
                psth_norm_ca1_mec(kk,:) = psth_ca1_mec(ii,:)/FR;

                FR_ratio_ca1(kk) = nanmean(psth_ca1(ii,61:110))/nanmean(psth_ca1(ii,1:60));
                FR_ratio_mec(kk) = nanmean(psth_mec(ii,61:110))/nanmean(psth_mec(ii,1:60));
                FR_ratio_both(kk) = nanmean(psth_ca1_mec(ii,61:110))/nanmean(psth_ca1_mec(ii,1:60));
                kk = kk+1;
            end
        end
        [~, idxsort_ca1] = sort(FR_ratio_ca1);
        [~, idxsort_mec] = sort(FR_ratio_mec);
        [~, idxsort_both] = sort(FR_ratio_both);
        psth_norm_ca1 = psth_norm_ca1(idxsort_ca1,:);
        psth_norm_mec = psth_norm_mec(idxsort_mec,:);
        psth_norm_ca1_mec = psth_norm_ca1_mec(idxsort_both,:);
        Rates.FR_ratio_ca1 = FR_ratio_ca1;
        Rates.FR_ratio_mec = FR_ratio_mec;
        %Rates.FR_ratio_both = FR_ratio_both;

        subplot(4,3,3*(class-1)+1)
        imagesc(psth_norm_ca1)
        title(strcat('CA1 Silencing,', CellClass{class}))
        subplot(4,3,3*(class-1)+2)
        imagesc(psth_norm_mec)
        title(strcat('mEC Silencing,', CellClass{class}))
        subplot(4,3,3*(class-1)+3)
        imagesc(psth_norm_ca1_mec)
        title(strcat('CA1 & mEC Silencing,', CellClass{class}))
              
        if class == 1
            subplot(4,3,7)
            a = nanmean(psth_norm_ca1,1);
            % %plot(a,'k')
            plot(a/mean(a(1:60)),'Color',[0.9856 0.7372 0.2537],'LineWidth',1.5)
            %plot(a/mean(a(1:60)),'Color',[0.2305 0.2510 0.6173],'LineWidth',2)
            hold on
            a = nanmean(psth_norm_mec,1);
            %plot(a,'r')
            plot(a/mean(a(1:60)),'Color',[0.1986 0.7214 0.6310],'LineWidth',1.5)
            a = nanmean(psth_norm_ca1_mec,1);
            % %plot(a,'b')
            plot(a/mean(a(1:60)),'Color',[0.4733 0.3374 0.2798],'LineWidth',1.5)
            %plot(a/mean(a(1:60)),'Color',[0.2140 0.2140 0.2140],'LineWidth',2)
            xlim([0 120])
            ylim([0 1.5]) 
            
            subplot(4,3,10)
            colormap(col)
            if cc == 1
                binnums = 10;
            else
                binnums = 3;
            end
                
            nhist(Rates,'binfactor',binnums,'maxx',5,'proportion','samebins','color','colormap')
            colormap default
%             scatter(FR_ratio_ca1,zeros(1,length(FR_ratio_ca1)),[],[0.9856 0.7372 0.2537],'filled')
%             hold on
%             scatter(FR_ratio_mec,zeros(1,length(FR_ratio_mec))+1,[],[0.1986 0.7214 0.6310],'filled')
%             scatter(FR_ratio_both,zeros(1,length(FR_ratio_both))+2,[],[0.4733 0.3374 0.2798],'filled')
%             ylim([-1 3])
        else
            subplot(4,3,8)
            a = nanmean(psth_norm_ca1,1);
            % %plot(a,'k')
            plot(a/mean(a(1:60)),'Color',[0.9856 0.7372 0.2537],'LineWidth',1.5)
            %plot(a/mean(a(1:60)),'Color',[0.2305 0.2510 0.6173],'LineWidth',2)
            hold on
            a = nanmean(psth_norm_mec,1);
            %plot(a,'r')
            plot(a/mean(a(1:60)),'Color',[0.1986 0.7214 0.6310],'LineWidth',1.5)
            a = nanmean(psth_norm_ca1_mec,1);
            % %plot(a,'b')
            plot(a/mean(a(1:60)),'Color',[0.4733 0.3374 0.2798],'LineWidth',1.5)
            %plot(a/mean(a(1:60)),'Color',[0.2140 0.2140 0.2140],'LineWidth',2)
            xlim([0 120])
            ylim([0 1.5])  
            
            subplot(4,3,11)
            colormap(col)            
            nhist(Rates,'binfactor',3,'maxx',6,'proportion','samebins','color','colormap')
            colormap default
%             scatter(FR_ratio_ca1,zeros(1,length(FR_ratio_ca1)),[],[0.9856 0.7372 0.2537],'filled')
%             hold on
%             scatter(FR_ratio_mec,zeros(1,length(FR_ratio_mec))+1,[],[0.1986 0.7214 0.6310],'filled')
%             scatter(FR_ratio_both,zeros(1,length(FR_ratio_both))+2,[],[0.4733 0.3374 0.2798],'filled')
%             ylim([-1 3])
        end
        set(gcf,'renderer','Painters')
        clear psth_norm_ca1 psth_norm_mec psth_norm_ca1_mec FR_ratio_ca1 FR_ratio_mec FR_ratio_both
    end
end
end

