function compileMiceCSDPowerCorr(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
tag = 'CA3';
chName =  {'pyramidal','radiatum','slm','ml'};
   
if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ21\Final'};%,'IZ15\Final','IZ30\Final'};%
    reg = {'mEC','CA1','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ18\Final','IZ20\Final','IZ21\Final','IZ31\Final','IZ27\Final','IZ28\Final','IZ33\Final'};  
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
targetReg = {'STEM', 'RETURN'};
freqBand = {'theta','sg','mg','hfo'};

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
        efields = fieldnames(sessionPulses);
        
        file = dir(('*.PowerSpectrumProfileCSD_6_12channelinfo.mat'));
        Pow.theta = load(file.name);
        file = dir(('*.PowerSpectrumProfileCSD_25_45channelinfo.mat'));
        Pow.sg = load(file.name);
        file = dir(('*.PowerSpectrumProfileCSD_60_120channelinfo.mat'));
        Pow.mg = load(file.name);
        file = dir(('*.PowerSpectrumProfileCSD_150_250channelinfo.mat'));
        Pow.hfo = load(file.name); 
        
        %Find linearized position for each timestamp
        tCSD = Pow.theta.powerProfile.time(1,:);
        tBehav = behavior.timestamps;
        linMap = behavior.position.lin;
        linMapDouble = linMap+170;
        idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
        linMapNew = linMap;
        linMapNew(idxMap) =  linMapDouble(idxMap);           
        linMapNew = linMapNew-75;
        
        for tt = 1:length(tCSD)
            [tmin,idx] = min(abs(tBehav-tCSD(tt)));
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
                    else
                       regionStim = sessionPulses.(efields{jj}).region;
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
                       pwr = Pow.(freqBand{ff}).powerProfile.power(:,idxB)';
                       corrPwr = corr(pwr,'Rows','pairwise','Type','Spearman');
                       if size(corrPwr,2) == 3
                          corrMatBase(1,bb,1:2) = reshape(corrPwr(1,2:size(corrPwr,2)),[1,1,(size(corrPwr,2)-1)]);
                          corrMatBase(1,bb,3) = nan;
                       else
                          corrMatBase(1,bb,:) = reshape(corrPwr(1,2:size(corrPwr,2)),[1,1,(size(corrPwr,2)-1)]);                           
                       end
                   else
                       corrMatBase(1,bb,1:3) = nan;
                   end
                   if sum(idxS)>0
                       pwr = Pow.(freqBand{ff}).powerProfile.power(:,idxS)';
                       corrPwr = corr(pwr,'Rows','pairwise','Type','Spearman');
                       if size(corrPwr,2) == 3
                          corrMatStim(1,bb,1:2) = reshape(corrPwr(1,2:size(corrPwr,2)),[1,1,(size(corrPwr,2)-1)]);
                          corrMatStim(1,bb,3) = nan;
                       else
                          corrMatStim(1,bb,:) = reshape(corrPwr(1,2:size(corrPwr,2)),[1,1,(size(corrPwr,2)-1)]);                           
                       end
                   else
                       corrMatStim(1,bb,1:3) = nan;
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

colMat = [52/243 52/243 52/243;...
    8/243 133/243 166/243;...
    187/243 86/243 149/243];

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
            ylabel('Power correlation')
            xlabel('Maze location')
            ylim([-1 1])
        end
    end
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\PowerCorr',tag,'.png'));
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\PowerCorr',tag,'.eps'),'epsc');
    saveas(gcf,strcat(parentDir,'Compiled\CSDPowerCorr\PowerCorr',tag,'.fig'));    
end

