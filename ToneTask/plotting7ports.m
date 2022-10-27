gain =[420/55, 420/130, 420/210, 420/290, 420/370, 420/420];

for cellNum = 59%1:249    
    for pf = 1:(size(behavTrials.timestamps,1)-1)    
        [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
        vy = tracking.position.vy;
        idx = vy>0.2 & idx;
        positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
    %     %Velocities
    %     if max(tracking.position.y(idx))<68
    %         yPos = tracking.position.y(idx);
    %         vel = tracking.position.v(idx);
    %         [~,idxY1] = min(abs(yPos-57));
    %         positions.vel(pf,:) = [vel(idxY1) nan nan];
    %     elseif max(tracking.position.y(idx))<90
    %         yPos = tracking.position.y(idx);
    %         vel = tracking.position.v(idx);
    %         [~,idxY1] = min(abs(yPos-57));
    %         [~,idxY2] = min(abs(yPos-83));            
    %         positions.vel(pf,:) = [vel(idxY1) vel(idxY2) nan];       
    %     else
    %         yPos = tracking.position.y(idx);
    %         vel = tracking.position.v(idx);
    %         [~,idxY1] = min(abs(yPos-57));
    %         [~,idxY2] = min(abs(yPos-83));  
    %         [~,idxY3] = min(abs(yPos-108)); 
    %         positions.vel(pf,:) = [vel(idxY1) vel(idxY2) vel(idxY3)];
    %     end
    end
    ctrl = 0;
    b = linspace(1,118,55);
    % First lin trials
    idx = find(behavTrials.linTrial);
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[2085 79 1057 863])
    set(gcf,'Color','white')
    subplot(7,2,1)
    for ii = 1:length(idx)-1
        plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
        hold on
        posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
        plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
        xlim([0 115])
    end
    %set(gca,'YDir','reverse')
    subplot(7,2,2)
    plot(b,firingMaps.forward.rateMaps{cellNum}{1},'Color',[0 0 0],'LineWidth',1.5);
    col = [176/243 223/243 229/243; 149/243 200/243 216/243; 137/243 207/243 240/243;
        70/243 130/243 180/243; 16/243 52/243 166/243;0/243 0/243 128/243];
    
    a = linspace(1000,25000,55);

    for kk = 1:6
        idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & behavTrials.toneGain(1:(end-1)) ==(kk-1));
        subplot(7,2,3+(2*(kk-1)))
        for ii = 1:length(idx)
            plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
            hold on
            posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
            plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
            xlim([0 115])
        end
        %set(gca,'YDir','reverse')
        if ctrl
            subplot(7,2,2)
            hold on
            plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);
            xlim([1 112])

            subplot(4,2,4)
            hold on
            plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);    
            xlim([1000 22000])
        else
            subplot(7,2,2)
            hold on        
            plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);
            xlim([1 118])
            subplot(7,2,4)
            hold on
            plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);      
            xlim([1000 25000])
        end
    end
    saveas(gcf,['RateMap',filesep ,'cell_' num2str(cellNum) '.png'],'png');
    saveas(gcf,['RateMap',filesep ,'cell_' num2str(cellNum) '.fig'],'fig');
  %  close all;

    %cell ids = 17, 48, 109, 186, 187,202 - sess 13
    %cellNum = 57;
    figure
    set(gcf,'Color','w')
    set(gcf,'Position',[2050 181 1585 762])
    plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1.2)
    hold on
    %plot(tracking.timestamps,tracking.position.v)
    scatter(tracking.timestamps(spikeData.posIdx{cellNum}),tracking.position.y(spikeData.posIdx{cellNum}),8,'r','filled')
    hold on
    scatter(behavTrials.timestamps((behavTrials.linTrial==1),2),ones(1,sum(behavTrials.linTrial==1))*120,25,[0.5 0.5 0.5],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==0),2),ones(1,sum(behavTrials.lickLoc==0))*130,25,[176/243 223/243 229/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==1),2),ones(1,sum(behavTrials.lickLoc==1))*130,25,[149/243 200/243 216/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==2),2),ones(1,sum(behavTrials.lickLoc==2))*130,25,[137/243 207/243 240/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==3),2),ones(1,sum(behavTrials.lickLoc==3))*130,25,[70/243 130/243 180/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==4),2),ones(1,sum(behavTrials.lickLoc==4))*130,25,[16/243 52/243 166/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.lickLoc==5),2),ones(1,sum(behavTrials.lickLoc==5))*130,25,[0/243 0/243 128/243],'filled')   
    scatter(behavTrials.timestamps((behavTrials.correct==1),2),ones(1,sum(behavTrials.correct==1))*140,25,[70/243 148/243 73/243],'filled')
    scatter(behavTrials.timestamps((behavTrials.correct==0 & behavTrials.linTrial==0),2),ones(1,sum(behavTrials.correct==0 & behavTrials.linTrial==0))*140,25,[187/243 86/243 149/243],'filled')
    %xlim([8250 10650])
    ylim([0 150])
    xlabel('Time(s)')
    ylabel('Position on track (cm)')
    box off

    saveas(gcf,['RateMap',filesep ,'Avgcell_' num2str(cellNum) '.png'],'png');
    saveas(gcf,['RateMap',filesep ,'Avgcell_' num2str(cellNum) '.fig'],'fig');
   % close all;
end
% for pf = 1:(size(behavTrials.timestamps,1)-1)    
%     if behavTrials.linTrial(pf)==1
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color',[0.5 0.5 0.5])
%     elseif behavTrials.correct(pf)==0
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color','m')
%     elseif behavTrials.correct(pf)==1
%         line([behavTrials.timestamps(pf,2) behavTrials.timestamps(pf,2)],[0 120],'Color','g')
%     end
% end

