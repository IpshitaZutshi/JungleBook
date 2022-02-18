function SessBehaviorLFPSpeed(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
pyrChPlus = 0;

if exist(strcat('Summ\LFPSpeed.mat'),'file') && ~force 
    disp('LFP Speed relationship already computed! Loading file.');
    load(strcat('Summ\LFPSpeed.mat'));
else

    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                LFPSpeed.specgram{rr,cc}{zz} = [];
                LFPSpeed.ts{rr,cc}{zz} = [];
                LFPSpeed.vel{rr,cc}{zz} = [];
            end
        end
    end

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.Behavior.mat'));
        load(file(1).name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        file = dir(('*.region.mat'));
        load(file.name);    
        pyrCh = region.CA1sp;        
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus);
            end
        end  
        lfp = bz_GetLFP(pyrChNew,'noPrompts', true);
        
        [~,~,~,vx,vy,~,~] = KalmanVel(behavior.position.x,behavior.position.y,behavior.timestamps,2);
        vel = sqrt(vx.^2+vy.^2);

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
                    specslope.specgram = nan;
                    specslope.timestamps = nan;
                    specslope.vel = nan;
                else
                    specslope = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[1 150],'spectype','fft','ints',events');
                    for pp = 1:length(specslope.timestamps)
                        events_tmp = abs(behavior.timestamps-specslope.timestamps(pp));
                        [~,idx] = min(events_tmp);
                        specslope.vel(pp,1) = vel(idx);
                    end                        
                end
                LFPSpeed.specgram{region,target}{zz} = [LFPSpeed.specgram{region,target}{zz}; specslope.specgram];
                LFPSpeed.ts{region,target}{zz} = [LFPSpeed.ts{region,target}{zz}; specslope.timestamps];
                LFPSpeed.vel{region,target}{zz} = [LFPSpeed.vel{region,target}{zz}; specslope.vel];
                
                clear specslope

            end        
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\LFPSpeed.mat'], 'LFPSpeed');
    end
end

freqs = logspace(log10(1),log10(150),200);
reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

speedRange = 0:1:25;

for ii = 2:length(reg)
     figure(ii)
     set(gcf,'Position',[50 50 1500 500])
     for jj = 1:length(target)
         for kk = 1:length(zone)      
             subplot(2,length(zone),length(zone)*(jj-1)+kk)  
             for ss = 1:(length(speedRange)-1)
                 idxCurr = LFPSpeed.vel{ii,jj}{1,kk}>=speedRange(ss) & LFPSpeed.vel{ii,jj}{1,kk} < speedRange(ss+1);
                 specmean(:,ss) = nanmean((LFPSpeed.specgram{ii,jj}{1,kk}(idxCurr,:)).^10)';
             end
             %specmean = log10(specmean);
             %contourf(speedRange(1:(end-1)),freqs,specmean,40,'LineColor','none')
             %specmean = zscore(specmean,0,'all');
             imagesc(speedRange(1:(end-1)),freqs,specmean)
             set(gca,'YDir','normal')
             set(gca,'YScale','log')
             colormap jet
             colorbar
             %caxis([0 2*10^7])
             ylim([5 150])
             xlabel('Running speed (cm/s)');
             ylabel('Frequency (Hz)');
             title(strcat(target(jj),'.',zone(kk)));
         end
     end
        
     saveas(figure(ii),strcat(expPath,'\Summ\LFPSpeed',reg{ii},'.png'));
     saveas(figure(ii),strcat(expPath,'\Summ\LFPSpeed',reg{ii},'.fig'));
     saveas(figure(ii),strcat(expPath,'\Summ\LFPSpeed',reg{ii},'.eps'));
end

end