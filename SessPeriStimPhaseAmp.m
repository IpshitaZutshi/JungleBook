function phaseAmp = SessPeriStimPhaseAmp(varargin)

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
phaserange =[5 15];
amprange = [25 300];

pyrChPlus = 5;

if exist(strcat('Summ\phaseAmpPyrChPlus',num2str(pyrChPlus),'.mat'),'file') && ~force 
    disp('phaseAmp already computed! Loading file.');
    load(strcat('Summ\phaseAmpPyrChPlus',num2str(pyrChPlus),'.mat'));
else
    for rr = 1:3
        for cc = 1:2
            phaseAmp{rr,cc}{1} = [];
            phaseAmp{rr,cc}{2} = [];
            phaseAmp{rr,cc}{3} = [];
            phaseAmp{rr,cc}{4} = [];
            phaseAmp{rr,cc}{5} = [];        
            phaseAmp{rr,cc}{6} = [];        
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
        pyrCh = region.CA1sp;
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus);
                lfp = bz_GetLFP(sessionInfo.AnatGrps(ch).Channels,'noPrompts', true);
                csd  = bz_CSD(lfp,'plotCSD',false,'plotLFP',false);
                clear lfp
            end
        end    
        
        csd.data = csd.data(:,csd.channels == pyrChNew);
        csd.channels = csd.channels(csd.channels == pyrChNew);
        lfp = csd;
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  

            %Extract relevant intervals for cross-frequency coupling - 4 cross
            %modulograms
            %First, no stim trials, return
            startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
            endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
            events = [startTS'; endTS'];
            comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false);  
            phaseAmp{region,target}{1} =  cat(3,phaseAmp{region,target}{1},comod.map);

            %No stim, stem
            startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
            endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
            events = [startTS';endTS'];
            comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false);  
            phaseAmp{region,target}{2} = cat(3,phaseAmp{region,target}{2},comod.map);

            %No stim, delay
            startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
            endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
            events = [startTS';endTS'];
            if sessionArmChoice.(efields{jj}).delay.dur >= 1
                comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false); 
            else 
                comod = nan;
            end
            phaseAmp{region,target}{3} = cat(3,phaseAmp{region,target}{3},comod.map);

            % Stim, return
            startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
            endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
            events = [startTS';endTS'];
            comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false); 
            phaseAmp{region,target}{4} = cat(3,phaseAmp{region,target}{4},comod.map);

            % Stim, stem
            startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
            endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
            events = [startTS';endTS'];     
            comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false); 
            phaseAmp{region,target}{5} = cat(3,phaseAmp{region,target}{5},comod.map);       

            %stim, delay
            startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
            endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
            events = [startTS';endTS'];
            if sessionArmChoice.(efields{jj}).delay.dur >= 1
                comod = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,'intervals',events','makePlot',false); 
            else
                comod = nan;
            end
            phaseAmp{region,target}{6} = cat(3,phaseAmp{region,target}{6},comod.map);

            clear rewardTS startDelay events
        end

    end
    
    if saveMat
        save([expPath '\Summ\' 'phaseAmpPyrChPlus' num2str(pyrChPlus) '.mat'], 'phaseAmp','comod');
    end
end

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ii = 1:size(phaseAmp,1)
    figure(ii)
    for jj = 1:size(phaseAmp,2)
        for kk = 1:length(zone)
            subplot(2,length(zone),length(zone)*(jj-1)+kk)  
            meancomod = nanmean(phaseAmp{ii,jj}{kk},3);
            imagesc(comod.phasecenters,log2(comod.amp_freq),meancomod');
            colormap jet
            hold on
            imagesc(comod.phasecenters+2*pi,log2(comod.amp_freq),meancomod')
            colormap jet
            plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(comod.amp_freq(end/2)),'k')
            xlabel(['Phase (',num2str(comod.phase_range),' Hz)']);
            ylabel('Amplitude Frequency (Hz)')
            LogScale('y',2)
            xlim([comod.phasecenters(1) comod.phasecenters(end)+2*pi]);
            colorbar
            axis xy
            title(strcat(target(jj),'.',zone(kk)));
        end
    end
    saveas(figure(ii),strcat(expPath,'\Summ\phaseAmpPyrChPlus',num2str(pyrChPlus),reg{ii},'.png'));
    saveas(figure(ii),strcat(expPath,'\Summ\phaseAmpPyrChPlus',num2str(pyrChPlus),reg{ii},'.eps'));
    saveas(figure(ii),strcat(expPath,'\Summ\phaseAmpPyrChPlus',num2str(pyrChPlus),reg{ii},'.fig'));
end

end