function SessBehaviorPeerPrediction(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'savePeerMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[7 12],@isnumeric);
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
savePeerMat = p.Results.savePeerMat;
force = p.Results.force;
plotfig = p.Results.plotfig;
passband = p.Results.passband;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');
pyrChPlus = 0;

if exist(strcat('Summ\ThetaCompression.mat'),'file') && ~force 
    disp('Theta compression already computed! Loading file.');
    load(strcat('Summ\ThetaCompression.mat'));
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:2
                peerPrediction.dev{rr,cc}{zz} = [];
                peerPrediction.devControl{rr,cc}{zz} = [];
            
            end
            thetaComp.region{rr,cc} = [];
            thetaComp.putativeCellType{rr,cc} = [];          
            thetaComp.sessCount{rr,cc} = [];
            sessCount{rr,cc} = 1;
        end
    end    

    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        file = dir(('*.Behavior.mat'));
        load(file(1).name);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        file = dir(('*.hippocampalLayers.channelinfo.mat'));
        load(file.name);    
        pyrCh = hippocampalLayers.oriens; 
   
        if ~isfield(session.channelTags,'RippleNoise')
            lfp = bz_GetLFP(pyrCh,'noPrompts', true);
        else
            refChannel = session.channelTags.RippleNoise.channels-1;
            lfp = bz_GetLFP([pyrCh refChannel],'noPrompts', true);
            lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
        end
        
        %filterlfp
        [b,a] = butter(3,[passband(1)/(lfp.samplingRate/2) passband(2)/(lfp.samplingRate/2)],'bandpass'); % order 3
        filt = FiltFiltM(b,a,double(lfp.data(:,1)));
        hilb = hilbert(filt);
        lfpphase = mod(angle(hilb),2*pi);
            
        if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
             spikes = bz_LoadPhy('noPrompts',true);
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
            arm = sessionArmChoice.(efields{jj}).visitedArm(2:end);
            choice = sessionArmChoice.(efields{jj}).choice(2:end);
                                
            for zz = 1:2
                
               spikeTimes = [];
               extraPredictors = [];
               %Extract relevant intervals to build the spike matrix
               fprintf(' ** Examining zone%3.i\n',zz);
                
               switch zz %% for some reason armchoice arms are FUTURE arm, not current - correct
                    case 1  %First, no stim trials   
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
                        events = [startTS endTS];
                    case 2  %stim
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
                        events = [startTS endTS];                       
               end
               
               % Build matrix, concatenate trials
               for trials = 1:size(events,1)
                   %position
                   behavBools = InIntervals(behavior.timestamps,events(trials,:));
                   behavPos = behavior.position.lin(behavBools);
                   
                   behavPosCorr = correctPos(behavPos);
                   behavTS = behavior.timestamps(behavBools);
                   
                   lfpBools = InIntervals(lfp.timestamps,events(trials,:));
                   lfpPhaseTrial = lfpphase(lfpBools);
                   lfpTS = lfp.timestamps(lfpBools);
                   
                   trialTime = round(events(trials,1)*1000):1:round(events(trials,2)*1000);
                   [~,lfpidx] = unique(round(lfpTS*1000));
                   lfpTS = lfpTS(lfpidx);
                   lfpPhaseTrial = lfpPhaseTrial(lfpidx);
                   lfpTime = interp1(round(lfpTS*1000),lfpPhaseTrial,trialTime);
                   behavTime = interp1(round(behavTS*1000),behavPosCorr,trialTime);
                   extraMat = [lfpTime;behavTime];
                   %extraMat = reshape(extraMat,[2,1,size(extraMat,2)]);
                   extraPredictors = catpad(2,extraPredictors,extraMat);
                   k = 1;
                   spkTimes = [];
                   %spikes
                   for unit = 1:length(spikes.times)
                       if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1
                           continue
                       else
                           bools = InIntervals(spikes.times{unit},events(trials,:));
                           tSp =spikes.times{unit}(bools);
                           spkTimes(k,1,1:length(trialTime)) = 0;
                           tSp = round(tSp*1000);
                           spkTimes(k,1,ismember(trialTime,tSp)) = 1;
                           k = k+1;
                       end
                   end
                   spikeTimes = catpad(3,spikeTimes,spkTimes);
                  % spikeTimes = catpad(2,spikeTimes,spkTimes);
               end  
               [dev devControl] = bz_peerPrediction(spikeTimes(:,:,2:end),0:150,extraPredictors);
               %[dev devControl] = bz_peerPrediction(spikeTimes,0:150,extraPredictors);
               
               save(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj},'\PeerPredictionMat',num2str(zz),'.mat'),'dev','devControl','-v7.3')
               
               dev = nanmean(dev,3);
               devControl = nanmean(devControl,3);
               
               peerPrediction.dev{region,target}{zz} = catpad(3,peerPrediction.dev{region,target}{zz},dev);
               peerPrediction.devControl{region,target}{zz} = catpad(3,peerPrediction.devControl{region,target}{zz},devControl);
               minDev = min(dev,[],1)-mean(dev,1);
               minDevControl = min(devControl,[],1)-mean(devControl,1);
               peerPrediction.assemblyStrength{region,target}{zz}  = [peerPrediction.assemblyStrength{region,target}{zz} (minDev./minDevControl)];
            end
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'ThetaCompression.mat'], 'thetaComp','-v7.3');
    end
end

if plotfig
    figure
    for ii = 2:3       
        for jj = 1
            for kk = 1:2
                idxtoKeep =  thetaComp.placefield_center{ii,jj}{kk}==1 & (thetaComp.placefield_difference{ii,jj}{kk}>-15 & thetaComp.placefield_difference{ii,jj}{kk}<15) & ...
                    (thetaComp.ccgs_place_offset{ii,jj}{kk}>-600 & thetaComp.ccgs_place_offset{ii,jj}{kk}<600);%...
                    %& thetaComp.region{ii,jj} == 1 & thetaComp.putativeCellType{ii,jj} == 1;
                
                subplot(4,6,2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (cm)')
                xlabel('Place field dist (ms)')
                    
                subplot(4,6,6+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.placefield_difference{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (cm)')
                xlabel('Theta timescale (ms)')

                subplot(4,6,12+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_time_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (ms)')
                xlabel('Theta timescale (ms)')    
                
                subplot(4,6,18+2*(ii-1)+kk)   
                scatter(thetaComp.ccgs_phase_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),25,'.')
                [R, pVal] = corr(thetaComp.ccgs_phase_offset{ii,jj}{kk}(idxtoKeep),thetaComp.ccgs_place_offset{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','Type','Spearman');
                title(strcat(' R:',num2str(R),' p:',num2str(pVal)))
                ylabel('Place field dist (ms)')
                xlabel('Theta phase (rad)')                                     
            end
        end
    end
end
end

function behavPosCorr = correctPos(linMap)

    linMapDouble = linMap+170;
    idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
    linMapNew = linMap;
    linMapNew(idxMap) =  linMapDouble(idxMap);           
    linMapNew = linMapNew-75;
    %remove sections when the mouse went backwards
    a = diff(linMapNew);
    idxA= find(a>150);
    idxB= find(a<-150);
    idxtoKeep = ones(length(linMapNew),1);
    %Remove the sections between idxA and idxB
    for pp = 1:length(idxA)
        temp = idxB((idxB-idxA(pp))>0);
        idxEnd = min(temp);
        idxtoKeep(idxA(pp):idxEnd) = 0;
    end
    idxtoKeep = logical(idxtoKeep);
    behavPosCorr = linMapNew(idxtoKeep);
end