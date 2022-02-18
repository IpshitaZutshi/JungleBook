function SessBehaviorPowerSpectrum(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
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

if exist(strcat('Summ\PowerSpectrumSlopePyrChPlus',num2str(pyrChPlus),'.mat'),'file') && ~force 
    disp('Power spectrum already computed! Loading file.');
    load(strcat('Summ\PowerSpectrumSlopePyrChPlus',num2str(pyrChPlus),'.mat'));
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                PSS.specData{rr,cc}{zz} = [];
                PSS.specTS{rr,cc}{zz} = [];
                PSS.specIntercept{rr,cc}{zz} = [];
                PSS.specGram{rr,cc}{zz} = [];
                PSS.specResid{rr,cc}{zz} = [];
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
%         file = dir(('*.hippocampalLayers.channelinfo.mat'));
%         load(file.name);   
        file = dir(('*.region.mat'));
        load(file.name);   
%         region.CA1sp = 8;
%         save(file.name,'region');
        
        %pyrCh = hippocampalLayers.pyramidal;
        pyrCh = region.CA1so;        
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus);
            end
        end  
        lfp = bz_GetLFP(pyrChNew,'noPrompts', true);

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
                    specslope.data = nan;
                    specslope.timestamps = nan;
                    specslope.intercept = nan;
                    specslope.specgram = nan;
                    specslope.resid = nan;
                    
                else
                    specslope = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[1 300],'spectype','fft','ints',events');
                    %specslope = bz_PowerSpectrumSlope(lfp,5,0.1,'frange',[1 300],'spectype','wavelet','ints',events');
                end
                PSS.specData{region,target}{zz} = catpad(3,PSS.specData{region,target}{zz},specslope.data);    
                PSS.specTS{region,target}{zz} = catpad(3,PSS.specTS{region,target}{zz},specslope.timestamps); 
                PSS.specIntercept{region,target}{zz} = catpad(3,PSS.specIntercept{region,target}{zz},specslope.intercept); 
                PSS.specGram{region,target}{zz} = catpad(3,PSS.specGram{region,target}{zz},specslope.specgram); 
                PSS.specResid{region,target}{zz} = catpad(3,PSS.specResid{region,target}{zz},specslope.resid);

            end        
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ContraPowerSpectrumSlopePyrChPlus' num2str(pyrChPlus) '.mat'], 'PSS');
    end
end

freqs = logspace(log10(1),log10(300),200);

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};
nf  = 1;

for ii = 1:length(reg)
     for jj = 1:length(target)
         for kk = 1:length(zone)      
            dat{kk}.slope = [];
            dat{kk}.intercept = [];
            data.spec = [];
            data.resid = [];
            for pp = 2:size(PSS.specData{ii,jj}{kk},3) % First array is nans
                dat{kk}.slope = [dat{kk}.slope PSS.specData{ii,jj}{kk}(:,:,pp)']; 
                dat{kk}.intercept = [dat{kk}.intercept PSS.specIntercept{ii,jj}{kk}(:,:,pp)']; 
                data.spec = [data.spec (PSS.specGram{ii,jj}{kk}(:,:,pp).^10)']; 
                data.resid = [data.resid PSS.specResid{ii,jj}{kk}(:,:,pp)']; 
            end                
            if kk < 4
                loc = kk;
                col = [85/243 85/243 85/243];
            else
                loc = kk-3;
                col = [8/243 133/243 161/243];
            end
            if ~isempty(dat{kk}.slope)
%                 figure(nf)
%                 set(gcf,'Position',[100 100 1700 600])
%                 subplot(2,6,6*(jj-1)+loc)  
%                 hold on
%               %  if kk ~= 3 && kk ~=6
%                     h = cdfplot(dat{kk}.slope);
%                     h.Color = col;
%                     h.LineWidth = 1;
%                     xlabel('Specgram slope')
%                % end
% 
%                 subplot(2,6,6*(jj-1)+loc+3)   
%                 hold on
%               %  if kk ~= 3 && kk ~=6
%                     h1 = cdfplot(dat{kk}.intercept);
%                     h1.Color = col;
%                     h1.LineWidth = 1;
%                     xlabel('Specgram intercept')       
%               %  end
% 
%                 figure(nf+1)
%                 set(gcf,'Position',[100 100 1700 600])
%                 subplot(2,6,6*(jj-1)+kk)  
%                 hold on
%                % if kk ~= 3 && kk ~=6
%                     imagesc(1:size(data.spec,2),freqs,data.spec)
%                     ylabel('f (Hz)')
%                     
%                     set(gca,'YDir','normal')
%                     set(gca, 'YScale', 'log')
%               %  end

                figure(nf)
                set(gcf,'Position',[100 100 1700 600])
                subplot(2,9,9*(jj-1)+loc)                 
                meanPower = nanmean(data.spec,2);
                stdPower = nanstd(data.spec,[],2)./sqrt(size(data.spec,2));
                dev1 = meanPower - stdPower;
                dev2 = meanPower + stdPower;
                hold on
                fill([freqs flip(freqs)],[dev1' flip(dev2)'],col,'FaceAlpha',.2,'EdgeColor','none');
                plot(freqs,meanPower,'color',col,'LineWidth',1.5);
                xlim([0 300])
                xlabel('Frequency')
                ylabel('Spec Power')    
                ylim([0 9*10^6])

                subplot(2,9,9*(jj-1)+loc+3)   
                hold on
                fill([freqs flip(freqs)],[dev1' flip(dev2)'],col,'FaceAlpha',.2,'EdgeColor','none');
                plot(freqs,meanPower,'color',col,'LineWidth',1.5);
                xlim([0 50])
                xlabel('Frequency')
                ylabel('Spec Power')    
                ylim([0 9*10^6])
                
                subplot(2,9,9*(jj-1)+loc+6)   
                hold on
                fill([freqs flip(freqs)],[dev1' flip(dev2)'],col,'FaceAlpha',.2,'EdgeColor','none');
                plot(freqs,meanPower,'color',col,'LineWidth',1.5);
                xlim([150 300])
                xlabel('Frequency')
                ylabel('Spec Power')    
                ylim([0 5*10^5])               

            end
         end
         
    end
%     saveas(figure(nf),strcat(expPath,'\Summ\PSS_SlopeIntercept',reg{ii},'.png'));
%     saveas(figure(nf),strcat(expPath,'\Summ\PSS_SlopeIntercept',reg{ii},'.eps'),'epsc');
%     saveas(figure(nf),strcat(expPath,'\Summ\PSS_SlopeIntercept',reg{ii},'.fig'));
%     
%     saveas(figure(nf+1),strcat(expPath,'\Summ\PSS_SpecGram',reg{ii},'.png'));
%     saveas(figure(nf+1),strcat(expPath,'\Summ\PSS_SpecGram',reg{ii},'.eps'),'epsc');
%     saveas(figure(nf+1),strcat(expPath,'\Summ\PSS_SpecGram',reg{ii},'.fig'));
    
    saveas(figure(nf),strcat(expPath,'\Summ\PSS_PowervsFreq',reg{ii},'.png'));
    saveas(figure(nf),strcat(expPath,'\Summ\PSS_PowervsFreq',reg{ii},'.eps'),'epsc');
    saveas(figure(nf),strcat(expPath,'\Summ\PSS_PowervsFreq',reg{ii},'.fig'));
    nf = nf+1;
end
close all
end