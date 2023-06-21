lin = behavior.position.lin;
linMap = lin;
linMapDouble = linMap+170;
idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
linMapNew = linMap;
linMapNew(idxMap) =  linMapDouble(idxMap);           
linMapNew = linMapNew-75;
linMapNew = round(linMapNew*(100/175));

efields = fields(sessionArmChoice);
fieldNum = 2;
cellNum = 3;

figure
set(gcf,'Renderer','painters')
% subplot(1,2,1)
% scatter(behavior.position.x,behavior.position.y,2,linMap,"filled");
% subplot(1,2,2)
scatter(behavior.position.x,behavior.position.y,2,linMapNew,"filled");

ts = sessionPulses.(efields{fieldNum}).timestamps(sessionPulses.(efields{fieldNum}).stim==1);

idx = [];
for ii = 1:length(ts)
    [~,idx(ii)] = min(abs(behavior.timestamps-ts(ii)));
end

hold on
plot(behavior.position.x(idx),behavior.position.y(idx),'ko')
    
ts = sessionArmChoice.(efields{fieldNum}).timestamps(sessionPulses.(efields{fieldNum}).stim==1);

idx = [];
for ii = 1:length(ts)
    [~,idx(ii)] = min(abs(behavior.timestamps-ts(ii)));
end

hold on
plot(behavior.position.x(idx),behavior.position.y(idx),'ro')


 for unit = cellNum
     posIdx = [];
        trialNum{1} = find(sessionArmChoice.(efields{fieldNum}).visitedArm==1 & sessionArmChoice.(efields{fieldNum}).choice==1 & sessionPulses.(efields{fieldNum}).stim'==0);
        trialNum{2} = find(sessionArmChoice.(efields{fieldNum}).visitedArm==0 & sessionArmChoice.(efields{fieldNum}).choice==1 & sessionPulses.(efields{fieldNum}).stim'==0);
        trialNum{3} = find(sessionArmChoice.(efields{fieldNum}).visitedArm==1 & sessionArmChoice.(efields{fieldNum}).choice==1 & sessionPulses.(efields{fieldNum}).stim'==1);
        trialNum{4} = find(sessionArmChoice.(efields{fieldNum}).visitedArm==0 & sessionArmChoice.(efields{fieldNum}).choice==1 & sessionPulses.(efields{fieldNum}).stim'==1);
        trialNum{2} = trialNum{2}(1:(end-1));
        for ii = 1:4
            trialsInt = [sessionArmChoice.(efields{fieldNum}).timestamps(trialNum{ii}) sessionArmChoice.(efields{fieldNum}).timestamps(trialNum{ii}+1)];
            [idx] = InIntervals(spikes.times{unit},trialsInt); 
            tsBehav = spikes.times{unit}(idx);
%             [idx] = InIntervals(spkT,trialsInt); 
%             tsBehav = spkT(idx);

            if isempty(tsBehav)
                posIdx{ii}(tt) = [];
            else
                for tt = 1:length(tsBehav)
                    [~,closestIndex] = min(abs(behavior.timestamps-tsBehav(tt)));
                    posIdx{ii}(tt) = closestIndex;
                end
            end
        end
 end
 figure
 set(gcf,'Renderer','painters')
 subplot(2,1,1)
 plot(behavior.position.x,behavior.position.y,'color',[0.5 0.5 0.5])
 hold on
 scatter(behavior.position.x(posIdx{1}),behavior.position.y(posIdx{1}),5,'r',"filled")
 scatter(behavior.position.x(posIdx{2}),behavior.position.y(posIdx{2}),5,'b',"filled")
 scatter(behavior.position.x(posIdx{3}),behavior.position.y(posIdx{3}),5,'g',"filled")
 scatter(behavior.position.x(posIdx{4}),behavior.position.y(posIdx{4}),5,'m',"filled")
 
 subplot(2,1,2)
 plot(behavior.timestamps(behavior.masks.recording == fieldNum),linMapNew(behavior.masks.recording == fieldNum),'color',[0.5 0.5 0.5])
 hold on
 scatter(behavior.timestamps(posIdx{1}),linMapNew(posIdx{1}),5,'r',"filled")
 scatter(behavior.timestamps(posIdx{2}),linMapNew(posIdx{2}),5,'b',"filled")
 scatter(behavior.timestamps(posIdx{3}),linMapNew(posIdx{3}),5,'g',"filled")
 scatter(behavior.timestamps(posIdx{4}),linMapNew(posIdx{4}),5,'m',"filled")
 
 
 for kk = 1:44
     ts = behavior.timestamps(behavior.masks.trials == kk & behavior.masks.recording == fieldNum);
     lin = linMapNew(behavior.masks.trials == kk & behavior.masks.recording == fieldNum);
     positions{kk} = [ts lin];
 end
 
%  spik.times{1} = spikes.times{cellNum};
%  spik.UID = 1;
%  firingMaps = bz_firingMapAvg_IZ(positions,spik,'nBins',100,'plotFig',false,'saveMat',false);
%  
 YlGnBu=cbrewer('seq', 'YlGnBu', 11);
 figure
 set(gcf,'Renderer','painters')
 for kk = 1:4
    subplot(2,4,kk)
    rm= [];
    for ll = 1:length(trialNum{kk})
        rm = [rm;firingMaps.trialMaps{cellNum}{trialNum{kk}(ll)}];
    end
    imagesc(rm)
    colormap(YlGnBu)
    caxis([0 30])
    
    subplot(2,4,kk+4)
    plot(nanmean(rm,1))
     ylim([0 20])
 end
