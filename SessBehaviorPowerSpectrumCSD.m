function SessBehaviorPowerSpectrumCSD(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'isCA3',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
plotfig = p.Results.plotfig;
isCA3 = p.Results.isCA3;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
chName = {'pyramidal','radiatum','slm','ml'};


if exist(strcat('Summ\PowerSpectrumCSD.mat'),'file') && ~force 
    disp('Power spectrum already computed! Loading file.');
    load(strcat('Summ\PowerSpectrumCSD.mat'));
else
    for ch = 1:length(chName)
        for rr = 1:3
            for cc = 1:2
                for zz = 1:6
                    PSS.(chName{ch}).specData{rr,cc}{zz} = [];
                    PSS.(chName{ch}).specTS{rr,cc}{zz} = [];
                    PSS.(chName{ch}).specIntercept{rr,cc}{zz} = [];
                    PSS.(chName{ch}).specGram{rr,cc}{zz} = [];
                    PSS.(chName{ch}).specResid{rr,cc}{zz} = [];
                end
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
        file = dir(('*.hippocampalLayersCSD.channelinfo.mat'));
        load(file.name);   
        
        chan = [];
        for ch = 1:length(chName)
            if ~isnan(hippocampalLayers.(chName{ch}))
                chan(ch) = hippocampalLayers.(chName{ch});            
                chanIdx(ch) = find(hippocampalLayers.channelOrder==chan(ch));
            end
        end
        
        lfp = bz_GetLFP(hippocampalLayers.channelOrder,'noPrompts', true);
        csd = bz_CSDIZ(lfp);
       
        lfp = csd;
        lfp.channels = chan;
        lfp.data = csd.data(:,(chanIdx-1));
        
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
                for ch = 1:length(chan)
                    PSS.(chName{ch}).specData{region,target}{zz} = catpad(3,PSS.(chName{ch}).specData{region,target}{zz},specslope.data(:,ch));    
                    PSS.(chName{ch}).specTS{region,target}{zz} = catpad(3,PSS.(chName{ch}).specTS{region,target}{zz},specslope.timestamps); 
                    PSS.(chName{ch}).specIntercept{region,target}{zz} = catpad(3,PSS.(chName{ch}).specIntercept{region,target}{zz},specslope.intercept(:,ch)); 
                    PSS.(chName{ch}).specGram{region,target}{zz} = catpad(3,PSS.(chName{ch}).specGram{region,target}{zz},specslope.specgram(:,:,ch)); 
                    PSS.(chName{ch}).specResid{region,target}{zz} = catpad(3,PSS.(chName{ch}).specResid{region,target}{zz},specslope.resid(:,:,ch));
                end

            end        
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'PowerSpectrumCSD.mat'], 'PSS');
    end
end


if plotfig
    freqs = logspace(log10(1),log10(300),200);
    colMat = [85/243 85/243 85/243;...
    224/243 163/243 46/243;... 
    8/243 133/243 161/243;...
    56/243 61/243 150/243];  

    reg = {'CA3','mEC','Both'};
    zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
    target = {'STEM', 'RETURN'};
    nf  = 1;
    figure
    set(gcf,'Renderer','painters');  
    set(gcf,'Position',[100 100 1700 600])
    if ~isCA3
        for ch = 1:length(fieldnames(PSS))
            for kk =1:2
                specbase{kk}= [];
            end
            for ii = 1:length(reg)
                 for jj = 1%:length(target)
                     for kk = [1 2]   
                        spec = [];
                        for pp = 2:size(PSS.(chName{ch}).specData{ii,jj}{kk},3) % First array is nans
                            specbase{kk} = [specbase{kk} (PSS.(chName{ch}).specGram{ii,jj}{kk}(:,:,pp).^10)']; 
                            spec = [spec (PSS.(chName{ch}).specGram{ii,jj}{kk+3}(:,:,pp).^10)'];
                        end                
                        if ~isempty(spec)                                                
                            subplot(4,2,2*(ch-1)+kk)                 
                            meanPower = nanmean(spec,2);
                            stdPower = nanstd(spec,[],2)./sqrt(size(spec,2));
                            dev1 = meanPower - stdPower;
                            dev2 = meanPower + stdPower;
                            hold on
                            fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(ii+1,:),'FaceAlpha',.2,'EdgeColor','none');
                            plot(freqs,meanPower,'color',colMat(ii+1,:),'LineWidth',1.5);
                            xlim([0 150])
                            xlabel('Frequency')
                            ylabel('Spec Power')    
                            hold on
                            %ylim([0 4*10^5])            
                        end
                     end
                 end
            end
            for kk = 1:2
                if ~isempty(specbase{kk}) 
                    subplot(4,2,2*(ch-1)+kk)                 
                    meanPower = nanmean(specbase{kk},2);
                    stdPower = nanstd(specbase{kk},[],2)./sqrt(size(specbase{kk},2));
                    dev1 = meanPower - stdPower;
                    dev2 = meanPower + stdPower;
                    hold on
                    fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(1,:),'FaceAlpha',.2,'EdgeColor','none');
                    plot(freqs,meanPower,'color',colMat(1,:),'LineWidth',1.5);
                    xlim([0 150])
                    xlabel('Frequency')
                    ylabel('Spec Power')    
                    hold on
                end
            end
        end
    else
        for ch = 1:4
             for jj = 1%:length(target)
                 for kk = [1 2]   
                    specbase =[]; 
                    spec1 = [];
                    spec2 = [];
                    spec3 = [];
                    for pp = 2:size(PSS.(chName{ch}).specData{2,jj}{kk},3) % First array is nans
                        specbase = [specbase (PSS.(chName{ch}).specGram{2,jj}{kk}(:,:,pp).^10)']; 
                        spec1 = [spec1 (PSS.(chName{ch}).specGram{2,jj}{kk+3}(:,:,pp).^10)'];
                        spec2 = [spec2 (PSS.(chName{ch}).specGram{3,jj}{kk}(:,:,pp).^10)'];
                        spec3 = [spec3 (PSS.(chName{ch}).specGram{3,jj}{kk+3}(:,:,pp).^10)'];
                    end                
                    if ~isempty(spec1)                                                
                        subplot(4,2,2*(ch-1)+kk)                  
                        meanPower = nanmean(specbase,2);
                        stdPower = nanstd(specbase,[],2)./sqrt(size(specbase,2));
                        dev1 = meanPower - stdPower;
                        dev2 = meanPower + stdPower;
                        hold on
                        fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(1,:),'FaceAlpha',.2,'EdgeColor','none');
                        plot(freqs,meanPower,'color',colMat(1,:),'LineWidth',1.5);
                        xlim([0 150])
                        xlabel('Frequency')
                        ylabel('Spec Power')    
                        hold on
                        
                        meanPower = nanmean(spec1,2);
                        stdPower = nanstd(spec1,[],2)./sqrt(size(spec1,2));
                        dev1 = meanPower - stdPower;
                        dev2 = meanPower + stdPower;
                        hold on
                        fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(3,:),'FaceAlpha',.2,'EdgeColor','none');
                        plot(freqs,meanPower,'color',colMat(3,:),'LineWidth',1.5);
                        xlim([0 150])
                        xlabel('Frequency')
                        ylabel('Spec Power')    
                        hold on
                        
                        meanPower = nanmean(spec2,2);
                        stdPower = nanstd(spec2,[],2)./sqrt(size(spec2,2));
                        dev1 = meanPower - stdPower;
                        dev2 = meanPower + stdPower;
                        hold on
                        fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(2,:),'FaceAlpha',.2,'EdgeColor','none');
                        plot(freqs,meanPower,'color',colMat(2,:),'LineWidth',1.5);
                        xlim([0 150])
                        xlabel('Frequency')
                        ylabel('Spec Power')    
                        hold on
                        
                        meanPower = nanmean(spec3,2);
                        stdPower = nanstd(spec3,[],2)./sqrt(size(spec3,2));
                        dev1 = meanPower - stdPower;
                        dev2 = meanPower + stdPower;
                        hold on
                        fill([freqs flip(freqs)],[dev1' flip(dev2)'],colMat(1,:),'FaceAlpha',.2,'EdgeColor','none');
                        plot(freqs,meanPower,'color',colMat(4,:),'LineWidth',1.5);
                        xlim([0 150])
                        xlabel('Frequency')
                        ylabel('Spec Power')    
                        hold on
                        %ylim([0 4*10^5])            
                    end
                 end
             end
        end               
    end
    saveas(gcf,strcat(expPath,'\Summ\CSDPSS_PowervsFreq.png'));
    saveas(gcf,strcat(expPath,'\Summ\CSDPSS_PowervsFreq.eps'),'epsc');
    saveas(gcf,strcat(expPath,'\Summ\CSDPSS_PowervsFreq.fig'));
end
end