function compileMiceCSDPowerCoherence(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
tag = 'CA3';
chName =  {'oriens','pyramidal','radiatum','slm','ml'};
   
if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ21\Final'};%,'IZ15\Final','IZ30\Final'};%
    reg = {'mEC','CA1','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ18\Final','IZ20\Final','IZ24\Final','IZ25\Final','IZ26\Final',...
        'IZ21\Final','IZ31\Final','IZ27\Final','IZ28\Final','IZ33\Final'};  
    reg = {'mEC','CA1','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ33\Final'};%,'IZ32\Final','IZ29\Final'};
    reg = {'mEC','CA3','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ33\Saline','IZ34\Saline'};%%'IZ29\Saline','IZ32\Saline',
    reg = {'mEC','CA3','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'ipsimEC','contramEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
targetReg = {'STEM'};%, 'RETURN'};
freqBand = {'theta','sg','mg','hfo'};
bandRange = {[6 11],[25 45],[60 120],[120 200]};

stimID = [];
targetID = [];

for ff = 1:length(freqBand)
    corrCSD.(freqBand{ff}) = [];
end

for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    % Start collecting data
    for ii = 1:size(allSess,1)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); 
        
        file = dir(('*.Behavior.mat'));
        load(file(1).name);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);
        file = dir(('*.coherenceRad.channelinfo.mat'));
        load(file.name);        
        file = dir(('*.hippocampalLayersCSD.channelinfo.mat'));
        load(file.name);  
        chan = hippocampalLayers.all;
        for ch = 1:length(chName)
            if ~isnan(hippocampalLayers.(chName{ch}))
                chan(ch) = hippocampalLayers.(chName{ch});            
                chanIdx(ch) = find(hippocampalLayers.channelOrder==chan(ch));
            else
                chanIdx(ch) = nan;
            end
        end
        
        efields = fieldnames(sessionPulses);       
        tCoh = coherence.timebins;
        tBehav = behavior.timestamps;
        linMap = behavior.position.lin;
        linMapDouble = linMap+170;
        idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
        linMapNew = linMap;
        linMapNew(idxMap) =  linMapDouble(idxMap);           
        linMapNew = linMapNew-75;        
        
        for tt = 1:length(tCoh)
            [tmin,idx] = min(abs(tBehav-tCoh(tt)));
            if tmin <1
                tPos(tt) = floor(linMapNew(idx));
                tTrial(tt) = behavior.masks.trials(idx);
                tRec(tt) = behavior.masks.recording(idx);
            else
                tPos(tt) = nan;
                tTrial(tt) = nan;
                tRec(tt) = nan;                
            end           
        end
        
        for jj = 1:length(efields)
            
                if strcmp(tag,'CA3') ==1 || strcmp(tag,'CA3Saline')==1
                    if jj == 1
                        regionBase = 1; 
                        regionStim = 2;
                    elseif jj == 2
                        regionBase = 3; 
                        regionStim = 4;
                    end
                elseif strcmp(tag,'mEC') == 1
                    if sessionPulses.(efields{jj}).region==2
                        regionBase = 1; 
                        regionStim = 2; 
                    else 
                        continue;
                    end
                else
                    regionBase = 1; 
                    if sessionPulses.(efields{jj}).region==1
                       regionStim = 3; 
                    elseif sessionPulses.(efields{jj}).region==2
                       regionStim = 2;
                    elseif sessionPulses.(efields{jj}).region==3
                       regionStim = 4;
                    end
                end
                    
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return  
            idxStim = find(sessionPulses.(efields{jj}).stim == 1);
            idxNoStim = find(sessionPulses.(efields{jj}).stim == 0);
            
            idxSelStim = tRec==jj & ismember(tTrial,idxStim);
            idxSelNoStim = tRec==jj & ismember(tTrial,idxNoStim);
            
            for ff = 1:length(freqBand)
                for bb = 1:ceil(max(linMapNew))
                   idxB = tPos == bb & idxSelNoStim;
                   idxS = tPos == bb & idxSelStim;               
                   if sum(idxB)>0
                       freqIdx = coherence.frequency>bandRange{ff}(1) & coherence.frequency<bandRange{ff}(2); 
                       if isnan(chanIdx(5))
                           corrMatBase(1,bb,1:4) = nanmean(coherence.coherogram(freqIdx,idxB,chanIdx(1:4)),[1 2]);
                           corrMatBase(:,:,5) = nan;
                       else
                           corrMatBase(1,bb,1:5) = nanmean(coherence.coherogram(freqIdx,idxB,chanIdx(1:5)),[1 2]);
                       end
                   else
                       corrMatBase(1,bb,1:5) = nan;
                   end
                   if sum(idxS)>0
                       freqIdx = coherence.frequency>bandRange{ff}(1) & coherence.frequency<bandRange{ff}(2); 
                       if isnan(chanIdx(5))
                           corrMatStim(1,bb,1:4) = nanmean(coherence.coherogram(freqIdx,idxS,chanIdx(1:4)),[1 2]);
                           corrMatStim(:,:,5) = nan;
                       else
                           corrMatStim(1,bb,1:5) = nanmean(coherence.coherogram(freqIdx,idxS,chanIdx(1:5)),[1 2]);
                       end
                   else
                       corrMatStim(1,bb,1:5) = nan;
                   end
               end
               corrCSD.(freqBand{ff}) = [corrCSD.(freqBand{ff});corrMatBase;corrMatStim];
               if ff == 1
                    stimID = [stimID;regionBase;regionStim];
                    targetID = [targetID;target;target];
               end
            end          
        end
        clear tPos tTrial tRec
    end
end

colMat = [52/243 52/243 52/243;...% gray - or
    224/243 163/243 46/243;...%yellow - pyr
    8/243 133/243 166/243;...%teal- rad
    187/243 86/243 149/243;...%pink - slm
    56/243 61/243 150/243];%dark blue

for jj = 1:length(targetReg)
    figure(jj)
    set(gcf,'Renderer','painters')
    set(gcf,'Position', [1921 41 1920 963]);
    for ff = 1:length(freqBand)
        for mm = 1:length(unique(stimID))
            subplot(length(freqBand),length(unique(stimID)),length(unique(stimID))*(ff-1)+mm)
            for ii = 1:size(corrCSD.(freqBand{ff}),3)
                meanCorr = nanmean(corrCSD.(freqBand{ff})(stimID==mm & targetID==jj,:,ii),1);
                stdCorr = nanstd(corrCSD.(freqBand{ff})(stimID==mm & targetID==jj,:,ii),[],1);
                dev1 = meanCorr(5:174)-stdCorr(5:174);
                dev2 = meanCorr(5:174)+stdCorr(5:174);
                hold on
                nC =5:1:174;
                fill([nC flip(nC)],[dev1 flip(dev2)],colMat(ii,:),'FaceAlpha',.2,'EdgeColor','none');
                plot(nC,meanCorr(5:174),'color',colMat(ii,:),'LineWidth',1);                
            end 
            if mm == 1
                title(strcat('Baseline ',freqBand{ff}))
            else
                title(strcat(reg{mm-1},' ',freqBand{ff}))
            end
            xlim([0 175])
            ylabel('Coherence')
            xlabel('Maze location')
            ylim([0.4 1])
        end
    end
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\CoherenceRad',tag,'.png'));
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\CoherenceRad',tag,'.eps'),'epsc');
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\CoherenceRad',tag,'.fig'));    
end

