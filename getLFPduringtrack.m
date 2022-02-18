function getLFPduringtrack(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'numtrials',5,@isnumeric);
addParameter(p,'refChannel',[],@isnumeric);
addParameter(p,'pyrChPlus',[0 5 10 15],@isnumeric);
addParameter(p,'makePlots',true,@islogical);
parse(p,varargin{:});


expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
numtrials = p.Results.numtrials;
refChannel = p.Results.refChannel;
pyrChPlus = p.Results.pyrChPlus;
makePlots = p.Results.makePlots;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = pwd; % select folder
end

%% Load necessary files
cd(expPath);

if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
load([sessionInfo.FileName '.region.mat']);
pyrCh = region.CA1sp;
load([sessionInfo.FileName '.SessionPulses.Events.mat']);
load([sessionInfo.FileName '.SessionArmChoice.Events.mat']);
load([sessionInfo.FileName '.Behavior.mat']);

 if makePlots
     figure
     set(gcf,'renderer','painters')
     set(gcf,'Position',[100 100 900 1300])
 end
for kk = 1:length(pyrChPlus)
    if isempty(refChannel)
%         for ch = 1:size(sessionInfo.AnatGrps,2)
%             if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
%                 Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
%                 pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus(kk));   
%             end
%         end
        lfp = bz_GetLFP(pyrChPlus(kk),'noPrompts', true);   
        %lfp = bz_GetLFP(pyrChNew,'noPrompts', true);    
    else
        for ch = 1:size(sessionInfo.AnatGrps,2)
            if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
                Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
                pyrChNew = sessionInfo.AnatGrps(ch).Channels(Chstart+pyrChPlus(kk));   
            end
        end    
        lfp = bz_GetLFP([pyrChNew refChannel],'noPrompts', true);
        lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        lfp.data = lfp.data(:,1);
    end
    %% Calculate spectogram
    params.Fs = lfp.samplingRate; params.fpass = [2 100]; params.tapers = [3 5]; params.pad = 1;
    [S,t,f] = mtspecgramc(single(lfp.data),[2 1],params);
    S1 = log10(S); % in Db
    S_det= bsxfun(@minus,S1,polyval(polyfit(f,mean(S1,1),2),f)); % detrending

    %% Extract timestamps
    spec_arr = [];
    efields = fieldnames(sessionPulses);
    reg = {'CA1','mEC','Both'};
    zone = {'Stem','Return'};

    for jj = 1:length(efields)
        spec_stim = [];
        spec_bs = [];
        ts_stim = [];
        ts_bs = [];
        stem_stim = [];
        stem_bs = [];
        counts = 0;
        countb = 0;
         stim = sessionPulses.(efields{jj}).stim(1:(end-1));
         startTS = sessionArmChoice.(efields{jj}).timestamps;
         delayTS = sessionArmChoice.(efields{jj}).delay.timestamps(2,:);
    %      [~,idxstart] = min(abs(t-startTS(1)));
    %      [~,idxend] = min(abs(t-startTS(end)));

         for ss = 1:(size(startTS,1)-1)

             [~,idxstart] = min(abs(t-startTS(ss)));
             [~,idxend] = min(abs(t-startTS(ss+1)));
             [~,idxstem] = min(abs(t-delayTS(ss)));

             if stim(ss) == 1
    %             [~,idx] = min(abs(t-sessionPulses.(efields{jj}).timestamps(ss)));
    %             idxstim(ss,1) = idx-idxstart;
    %             idxstim(ss,2) = idxend-idx;
                if counts<numtrials 
                    spec_stim = [spec_stim S_det(idxstart:idxend,:)'];
                    ts_stim = [ts_stim; t(idxend)-t(idxstart)];
                    stem_stim = [stem_stim; t(idxstem)-t(idxstart)];
                end
                counts = counts+1;
             else
    %              idxstim(ss,1) = 1;
    %              idxstim(ss,2) = 1;
                if countb<numtrials     
                    spec_bs = [spec_bs S_det(idxstart:idxend,:)'];
                    ts_bs = [ts_bs; t(idxend)-t(idxstart)];
                    stem_bs = [stem_bs; t(idxstem)-t(idxstart)];
                end
                countb = countb+1;
             end
    %          if isempty(spec_arr)
    %              spec_arr = S_det(idxstart:idxend,:)';
    %          else
    %              spec_arr = catpad(3,spec_arr,S_det(idxstart:idxend,:)');
    %          end
         end

         %Align stimulus timebin in spectogram
    %      minstart = min(idxstim(idxstim(:,1)>1,1));
    %      minend = min(idxstim(idxstim(:,2)>1,2));
    %      
         if makePlots
             subplot(length(pyrChPlus),2*length(efields),(length(efields)*2*(kk-1))+(2*(jj-1)+1))
             imagesc(1:size(spec_stim,2),f,spec_stim,[-1 1]);
             set(gca,'YDir','normal'); ylabel('Freqs');
             hold on
             for lin = 1:(length(ts_stim)-1)
                line([sum(ts_stim(1:lin)) sum(ts_stim(1:lin))],[2 100],'Color','k')
                if lin ==1 
                    line([stem_stim(lin) stem_stim(lin)],[2 100],'Color','b')   
                else
                    line([sum(ts_stim(1:(lin-1)))+stem_stim(lin) sum(ts_stim(1:(lin-1)))+stem_stim(lin)],[2 100],'Color','b')  
                end
             end
             
             title(strcat('stim ',reg(sessionPulses.(efields{jj}).region),' ',zone(sessionPulses.(efields{jj}).target),' ',num2str(pyrChPlus(kk))))

             subplot(length(pyrChPlus),2*length(efields),(length(efields)*2*(kk-1))+(2*(jj-1)+2))
             imagesc(1:size(spec_bs,2),f,spec_bs,[-1 1]);
             set(gca,'YDir','normal'); ylabel('Freqs');
             title(strcat('base ',reg(sessionPulses.(efields{jj}).region),' ',zone(sessionPulses.(efields{jj}).target),' ',num2str(pyrChPlus(kk))))
             hold on
             for lin = 1:length(ts_bs)
                line([sum(ts_bs(1:lin)) sum(ts_bs(1:lin))],[2 100],'Color','k')
                if lin ==1 
                    line([stem_bs(lin) stem_bs(lin)],[2 100],'Color','b')   
                else
                    line([sum(ts_bs(1:(lin-1)))+stem_bs(lin) sum(ts_bs(1:(lin-1)))+stem_bs(lin)],[2 100],'Color','b')  
                end
             end
         end

    end
end
if makePlots
    saveas(gcf,strcat('Summ\BehaviorSpectogramSess.png'));
    saveas(gcf,strcat('Summ\BehaviorSpectogramSess.eps'),'epsc');
    saveas(gcf,strcat('Summ\BehaviorSpectogramSess.fig'));
end
end