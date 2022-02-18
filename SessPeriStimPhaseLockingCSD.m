function ModData = SessPeriStimPhaseLockingCSD(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'nbins',180,@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
nbins = p.Results.nbins;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
passband = [5 12];
pyrChPlus = 0;

if exist(strcat('Summ\PhaseLocking.mat'),'file') && ~force 
    disp('Phase Locking already computed! Loading file.');
    load(strcat('Summ\PhaseLocking.mat'));
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                phaseData.phasedistros{rr,cc}{zz} = [];
                phaseData.p{rr,cc}{zz} = [];
                phaseData.m{rr,cc}{zz} = [];
                phaseData.r{rr,cc}{zz} = [];
            end
        end
    end    

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        file = dir(('*.region.mat'));
        load(file.name);  
        file = dir(('*.csd.mat'));
        load(file.name);  
        pyrCh = 9;%region.CA1sp; IZ8-19, IZ9-56, IZ10-9
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                %pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus);
            end
        end    
        
        csd.data = csd.data(:,(Chstart+pyrChPlus-1));
        csd.channels = csd.channels((Chstart+pyrChPlus-1));
        lfp = csd;
        
        if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
             spikes = bz_LoadPhy('noPrompts',true);
             spikes.region = 'CA1';
        else
            continue;
        end

        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

            for zz = 1:6
                %Extract relevant intervals for cross-frequency coupling - 4 cross
                %modulograms
                switch zz
                    case 1  %First, no stim trials, return        
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        events = [startTS'; endTS'];
                    case 2  %No stim, stem
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
                        events = [startTS';endTS'];
                    case 3 %No stim, delay
                        startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                        endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
                        events = [startTS';endTS'];  
                    case 4  % Stim, return
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        events = [startTS';endTS'];                    
                    case 5   % Stim, stem
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
                        events = [startTS';endTS'];                      
                    case 6    %stim, delay
                        startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                        endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
                        events = [startTS';endTS'];
                end
                
                if (zz == 3 || zz == 6) && sessionArmChoice.(efields{jj}).delay.dur < 1 
                    phaseLocking.phasedistros = nan;
                    phaseLocking.phasestats.p = nan;
                    phaseLocking.phasestats.r = nan;
                    phaseLocking.phasestats.m = nan;
                else
                    phaseLocking = bz_PhaseModulation_IZ(spikes,lfp,'passband',passband,'intervals',events','numbins',nbins,'powerThresh',1,'plotting',false);
                end
                phaseData.phasedistros{region,target}{zz} = catpad(3,phaseData.phasedistros{region,target}{zz},phaseLocking.phasedistros);    
                phaseData.p{region,target}{zz} = catpad(3,phaseData.p{region,target}{zz},phaseLocking.phasestats.p); 
                phaseData.r{region,target}{zz} = catpad(3,phaseData.r{region,target}{zz},phaseLocking.phasestats.r); 
                phaseData.m{region,target}{zz} = catpad(3,phaseData.m{region,target}{zz},phaseLocking.phasestats.m); 

            end   
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'PhaseLocking.mat'], 'phaseData');
    end
end

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

% Angle bins
binned = linspace(0,2*pi,nbins+1)';binned(end) = [];
binSize = binned(2)-binned(1);
binned = binned + binSize/2;

for ii = 1:length(reg)
    figure(ii)
    set(gcf,'Position',[100 100 1700 600])    
    figure(ii+3)
    set(gcf,'Position',[100 100 1700 600])    
    for jj = 1:length(target)
         for kk = 1:length(zone)      
            data(kk).m = [];
            data(kk).p = [];
            data(kk).r = [];
            data(kk).phasedistros= [];

            for pp = 2:size(phaseData.p{ii,jj}{kk},3) % First array is nans
                data(kk).m  = [data(kk).m  phaseData.m{ii,jj}{kk}(:,:,pp)]; 
                data(kk).p = [data(kk).p phaseData.p{ii,jj}{kk}(:,:,pp)]; 
                data(kk).r = [data(kk).r phaseData.r{ii,jj}{kk}(:,:,pp)]; 
                data(kk).phasedistros = [data(kk).phasedistros phaseData.phasedistros{ii,jj}{kk}(:,:,pp)];
            end   
            if size(data(kk).phasedistros,1)==180
                for bb = 1:20
                    phasedist(bb,:) = nanmean(data(kk).phasedistros((9*(bb-1)+1):(9*bb),:),1);
                end
                figure(ii)
                hax = subplot(2,6,6*(jj-1)+kk);
                bins = [0:18:360];
                nanRows = mean(~isnan(phasedist))>0;
                imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),zscore(phasedist(:,nanRows),0,1)')
                %imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),phasedist(:,nanRows)')
                hold on
                caxis([-3 3])            
                %imagesc((binned*180/pi)+360,1:size(data(kk).phasedistros,2),zscore(data(kk).phasedistros,0,1)')
                xlim([0 360])
                set(hax,'XTick',[0 90 180 270 360]) 
                plot([0:360],cos((pi/180)*[0:360])*0.05*size(phasedist(:,nanRows),2)+0.95*size(phasedist(:,nanRows),2),'color',[1 1 1],'LineWidth',1.5)
                set(gca,'YDir','normal')
                xlabel('Theta phase');
                ylabel('Neuron ID');            
                title(strcat(target(jj),'.',zone(kk)));            

                figure(ii+3)
                subplot(1,2,jj)
                switch kk
                    case 1
                        col = [85/243 85/243 85/243];
                        lstyle = '-';
                    case 2
                        col = [85/243 85/243 85/243];
                        lstyle = '--';
                    case 3
                        col = [85/243 85/243 85/243];
                        lstyle = ':';
                    case 4
                        col = [8/243 133/243 161/243];
                        lstyle = '-';                    
                    case 5
                        col = [8/243 133/243 161/243];
                        lstyle = '--';                    
                    case 6
                        col = [8/243 133/243 161/243];
                        lstyle = ':';      
                end
                avgDist= nanmean(phasedist(:,nanRows),2);
                stderr = nanstd(phasedist(:,nanRows),0,2)/sqrt(size(phasedist(:,nanRows),2));
                fill([bins(2:end)'; fliplr(bins(2:end))'],[avgDist-stderr;flipud(avgDist+stderr)],col,'linestyle','none','FaceAlpha',0.2);
                hold on            
                plot(bins(2:end),avgDist,'Color',col,'LineWidth',1.5,'LineStyle',lstyle)
                hold on
                xlabel('Theta phase');
                clear phasedist            
            end
        end
             
    %          meanRes = [];pVal = [];meanAng = [];
    %          meanAng = [data(1).m;data(2).m;data(3).m;data(4).m;data(5).m;data(6).m];
    %          pVal = [data(1).p;data(2).p;data(3).p;data(4).p;data(5).p;data(6).p];
    %          meanRes = [data(1).r;data(2).r;data(3).r;data(4).r;data(5).r;data(6).r];
    end
    saveas(figure(ii),strcat(expPath,'\Summ\PhaseLocking',reg{ii},'.png'));
    saveas(figure(ii),strcat(expPath,'\Summ\PhaseLocking',reg{ii},'.eps'),'epsc');
    saveas(figure(ii),strcat(expPath,'\Summ\PhaseLocking',reg{ii},'.fig'));

    saveas(figure(ii+3),strcat(expPath,'\Summ\AvgPhaseLocking',reg{ii},'.png'));
    saveas(figure(ii+3),strcat(expPath,'\Summ\AvgPhaseLocking',reg{ii},'.eps'),'epsc');
    saveas(figure(ii+3),strcat(expPath,'\Summ\AvgPhaseLocking',reg{ii},'.fig'));

end

end