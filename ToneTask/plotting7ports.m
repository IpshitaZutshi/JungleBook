function plotting7ports(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveLoc',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveLoc = p.Results.saveLoc;

if isempty(saveLoc)
    saveLoc = strcat(basepath,'\Maps');
    if ~isfolder('Maps')
        mkdir('Maps')
    end    
end
%% Deal with inputs
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikeData.cellinfo.mat']))
    file = dir([basepath filesep '*.spikeData.cellinfo.mat']);
    load(file.name);
end

if ~isempty(dir([basepath filesep '*.rateMapsAvg.cellinfo.mat']))
    file = dir([basepath filesep '*.rateMapsAvg.cellinfo.mat']);
    load(file.name);
end


[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
b = linspace(0,125,50);
a = linspace(2000,22000,50);

for cellNum = 1:length(spikeData.pos)   
    for pf = 1:(size(behavTrials.timestamps,1)-1)    
        [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
      %  vy = tracking.position.vy;
      %  idx2 = vy>0.2 & idx;
        positions.forward{pf} = [tracking.position.x(idx) tracking.position.y(idx) tracking.position.v(idx) find(idx==1)];   
    end
    
    % First lin trials
    linIdx = find(behavTrials.linTrial==1);
    jumpLin = find(diff(linIdx)>1);  
    if isempty(jumpLin)
        idx = linIdx;
    else
        idx = linIdx(1:jumpLin);
    end
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[2085 79 1057 863])
    set(gcf,'Color','white')
    subplot(8,5,1)
    for ii = 1:length(idx)-1
        plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
        hold on
        posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
        plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
        xlim([0 122])
    end    
    if isempty(jumpLin)
        idx = [];
    else
        idx = linIdx(jumpLin+1:end);
    end
    subplot(8,5,6) % Second half of the lin trials
    for ii = 1:length(idx)-1
        plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
        hold on
        posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
        plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
        xlim([0 122])
    end
    
    %set(gca,'YDir','reverse')
    subplot(8,5,5)
    plot(b,firingMaps.forward.rateMaps{cellNum}{1},'Color',[0 0 0],'LineWidth',1.5);
    hold on
    subplot(8,5,5)
    plot(b,firingMaps.forward.rateMaps{cellNum}{28},'Color',[0.5 0.5 0.5],'LineWidth',1.5);
    
    col = [176/243 223/243 229/243; 149/243 200/243 216/243; 137/243 207/243 240/243;
        70/243 130/243 180/243; 16/243 52/243 166/243;0/243 0/243 128/243];   
    
    %Ratemaps separated by stim and no stim
    for ss = 1:2
        for kk = 1:6
            idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & ...
                behavTrials.toneGain(1:(end-1)) ==(kk-1) & behavTrials.stim(1:(end-1))==(ss-1));
            subplot(8,5,11+(5*(kk-1))+(2*(ss-1)))
            for ii = 1:length(idx)
                plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
                hold on
                posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
                plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
                xlim([0 115])
            end
            %set(gca,'YDir','reverse')
            subplot(8,5,10+(10*(ss-1)))
            hold on        
            plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk+(12*(ss-1))},'Color',col(kk,:),'LineWidth',1.5);
            xlim([1 118])
            subplot(8,5,15+(10*(ss-1)))
            hold on
            plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk+(12*(ss-1))},'Color',col(kk,:),'LineWidth',1.5);      
            xlim([1000 25000])
        end   

        for kk = 1:6
            idx = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==0 &...
                behavTrials.toneGain(1:(end-1)) ==(kk-1) & behavTrials.stim(1:(end-1))==(ss-1));
            subplot(8,5,12+(5*(kk-1))+(2*(ss-1)))
            for ii = 1:length(idx)
                plot(positions.forward{idx(ii)}(:,2),positions.forward{idx(ii)}(:,1),'Color',[0.5 0.5 0.5])
                hold on
                posID = ismember(positions.forward{idx(ii)}(:,4),spikeData.posIdx{cellNum});
                plot(positions.forward{idx(ii)}(posID,2),positions.forward{idx(ii)}(posID,1),'r.')
                xlim([0 115])
            end
    %         %set(gca,'YDir','reverse')
    %         subplot(8,3,6)
    %         hold on        
    %         plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);
    %         xlim([1 118])
    %         subplot(8,3,9)
    %         hold on
    %         plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1.5);      
    %         xlim([1000 25000])
        end
    end
 
    saveas(gcf,[saveLoc,filesep ,'cell_', num2str(cellNum),'.png'],'png');
    saveas(gcf,[saveLoc,filesep ,'cell_', num2str(cellNum),'.fig'],'fig');
    saveas(gcf,[saveLoc,filesep ,'cell_', num2str(cellNum),'.eps'],'epsc');
    close all;
  
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
    scatter(behavTrials.timestamps((behavTrials.stim==1),2),ones(1,sum(behavTrials.stim==1))*150,25,'m','filled')
    scatter(behavTrials.timestamps((behavTrials.correct==0 & behavTrials.linTrial==0),2),ones(1,sum(behavTrials.correct==0 & behavTrials.linTrial==0))*140,25,[187/243 86/243 149/243],'filled')
    %xlim([8250 10650])
    ylim([0 150])
    xlabel('Time(s)')
    ylabel('Position on track (cm)')
    box off 
   saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.png'],'png');
   saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.fig'],'fig');
   saveas(gcf,[saveLoc,filesep ,'Avgcell_', num2str(cellNum),'.eps'],'epsc');
   close all;
end

