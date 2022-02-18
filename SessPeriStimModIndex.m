function ModData = SessPeriStimModIndex(varargin)

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
phaserange = 2:0.2:20;
amprange = 25:1:300;
pyrChPlus = 0;

if exist(strcat('Summ\ModIdx_PyrChPlus',num2str(pyrChPlus),'.mat'),'file') && ~force 
    disp('ModIdx already computed! Loading file.');
    load(strcat('Summ\ModIdx_PyrChPlus',num2str(pyrChPlus),'.mat'));
else
    for rr = 1:3
        for cc = 1:2
            ModData{rr,cc}{1} = [];
            ModData{rr,cc}{2} = [];
            ModData{rr,cc}{3} = [];
            ModData{rr,cc}{4} = [];
            ModData{rr,cc}{5} = [];        
            ModData{rr,cc}{6} = [];        
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

            %Extract relevant intervals for cross-frequency coupling - 4 cross
            %modulograms
            %First, no stim trials, return
            startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
            endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
            events = [startTS'; endTS'];
            comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            ModData{region,target}{1} =  cat(3,ModData{region,target}{1},comod);

            %No stim, stem
            startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
            endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
            events = [startTS';endTS'];
            comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            ModData{region,target}{2} = cat(3,ModData{region,target}{2},comod);

            %No stim, delay
            startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
            endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
            events = [startTS';endTS'];
            if sessionArmChoice.(efields{jj}).delay.dur >= 1
                comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            else 
                comod = nan;
            end
            ModData{region,target}{3} = cat(3,ModData{region,target}{3},comod);

            % Stim, return
            startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
            endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
            events = [startTS';endTS'];
            comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            ModData{region,target}{4} = cat(3,ModData{region,target}{4},comod);

            % Stim, stem
            startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
            endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
            events = [startTS';endTS'];     
            comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            ModData{region,target}{5} = cat(3,ModData{region,target}{5},comod);       

            %stim, delay
            startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
            endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
            events = [startTS';endTS'];
            if sessionArmChoice.(efields{jj}).delay.dur >= 1
                comod = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot',false,'phaserange',phaserange,'amprange',amprange);  
            else
                comod = nan;
            end
            ModData{region,target}{6} = cat(3,ModData{region,target}{6},comod);

            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ModIdx_PyrChPlus' num2str(pyrChPlus) '.mat'], 'ModData');
    end
end

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ii = 1:size(ModData,1)
    figure(ii)
    set(gcf,'Position',[100 100 1700 600])
    for jj = 1:size(ModData,2)
        for kk = 1:length(zone)
            subplot(2,length(zone),length(zone)*(jj-1)+kk)  
            meancomod = nanmean(ModData{ii,jj}{kk},3);
            imagesc(phaserange(2:end),amprange(2:end),meancomod)
            %imagesc(phaserange(2:end),log2(amprange(2:end)),meancomod)
            set(gca,'YDir','normal')
            colormap jet
            if kk <=3
                cmax = max(max(nanmean(ModData{ii,jj}{kk},3))); 
            else
                cmax = max(max(nanmean(ModData{ii,jj}{(kk-3)},3))); 
            end
            caxis([0 cmax])
            colorbar
            %LogScale('y',2)
            xlabel('Frequency phase');
            ylabel('Frequency amplitude');
            title(strcat(target(jj),'.',zone(kk)));
        end
    end
    saveas(figure(ii),strcat(expPath,'\Summ\ThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.png'));
    saveas(figure(ii),strcat(expPath,'\Summ\ThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.eps'));
    saveas(figure(ii),strcat(expPath,'\Summ\ThetaGammaCouplingPyrChPlus',num2str(pyrChPlus),reg{ii},'.fig'));
end

end