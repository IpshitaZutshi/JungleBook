gain =[400/220, 400/320, 400/400];

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    vy = tracking.position.vy;
    idx = vy>0.2 & idx;
    positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
    %Velocities
    if max(tracking.position.y(idx))<68
        yPos = tracking.position.y(idx);
        vel = tracking.position.v(idx);
        [~,idxY1] = min(abs(yPos-57));
        positions.vel(pf,:) = [vel(idxY1) nan nan];
    elseif max(tracking.position.y(idx))<90
        yPos = tracking.position.y(idx);
        vel = tracking.position.v(idx);
        [~,idxY1] = min(abs(yPos-57));
        [~,idxY2] = min(abs(yPos-83));            
        positions.vel(pf,:) = [vel(idxY1) vel(idxY2) nan];       
    else
        yPos = tracking.position.y(idx);
        vel = tracking.position.v(idx);
        [~,idxY1] = min(abs(yPos-57));
        [~,idxY2] = min(abs(yPos-83));  
        [~,idxY3] = min(abs(yPos-108)); 
        positions.vel(pf,:) = [vel(idxY1) vel(idxY2) vel(idxY3)];
    end
end
cellNum = 48;
ctrl = 0;
b = linspace(1,112,55);
% First lin trials
idx = find(behavTrials.linTrial);
figure
set(gcf,'Renderer','painters')
set(gcf,'Color','white')
subplot(4,2,1)
for ii = 1:length(idx)
    plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
    hold on
    posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
    plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
    xlim([0 115])
end
%set(gca,'YDir','reverse')
subplot(4,2,2)
plot(b,firingMaps.forward.rateMaps{cellNum}{1},'Color',[0.5 0.5 0.5],'LineWidth',1.5);
col = [103/243 189/243 170/243; 8/243 133/243 161/243; 56/243 61/243 150/243];

a = linspace(1000,22000,55);

for kk = 1:3
    idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & behavTrials.toneGain(1:(end-1)) ==(kk-1));
    subplot(4,2,3+(2*(kk-1)))
    for ii = 1:length(idx)
        plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
        hold on
        posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
        plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
        xlim([0 115])
    end
    %set(gca,'YDir','reverse')
    if ctrl
        subplot(4,2,2)
        hold on
        plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);
        xlim([1 112])
        
        subplot(4,2,4)
        hold on
        plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);    
        xlim([1000 22000])
    else
        subplot(4,2,2)
        hold on        
        plot(b,firingMaps.forward.rateMaps{cellNum}{4+kk},'Color',col(kk,:),'LineWidth',1.5);
        xlim([1 112])
        subplot(4,2,4)
        hold on
        plot(a,firingMaps.tone.rateMaps{cellNum}{4+kk},'Color',col(kk,:),'LineWidth',1.5);      
        xlim([1000 22000])
    end
end


%cell ids = 17, 48, 109, 186, 187,202 - sess 13
%cellNum = 57;
figure
set(gcf,'Color','w')
plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1.2)
hold on
%plot(tracking.timestamps,tracking.position.v)
scatter(tracking.timestamps(spikeData.posIdx{cellNum}),tracking.position.y(spikeData.posIdx{cellNum}),8,'r','filled')
hold on
scatter(behavTrials.timestamps((behavTrials.linTrial==1),2),ones(1,sum(behavTrials.linTrial==1))*120,25,[0.5 0.5 0.5],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==0),2),ones(1,sum(behavTrials.lickLoc==0))*130,25,[103/243 189/243 170/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==1),2),ones(1,sum(behavTrials.lickLoc==1))*130,25,[8/243 133/243 161/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==2),2),ones(1,sum(behavTrials.lickLoc==2))*130,25,[56/243 61/243 150/243],'filled')
scatter(behavTrials.timestamps((behavTrials.correct==1),2),ones(1,sum(behavTrials.correct==1))*140,25,[70/243 148/243 73/243],'filled')
scatter(behavTrials.timestamps((behavTrials.correct==0 & behavTrials.linTrial==0),2),ones(1,sum(behavTrials.correct==0 & behavTrials.linTrial==0))*140,25,[187/243 86/243 149/243],'filled')
xlim([8250 10650])
ylim([0 150])
xlabel('Time(s)')
ylabel('Position on track (cm)')
box off

% for pf = 1:(size(behavTrials.timestamps,1)-1)    
%     if behavTrials.linTrial(pf)==1
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color',[0.5 0.5 0.5])
%     elseif behavTrials.correct(pf)==0
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color','m')
%     elseif behavTrials.correct(pf)==1
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color','g')
%     end
% end

