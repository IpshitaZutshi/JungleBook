
function [armChoice] = getArmChoice(digitalIn)
% Compute T maze performance and gets timestamps
%
% USAGE
%
%   [armChoice, behavior] = getArmChoice(varargin)
%
% INPUTS
% 
% digitalIn                     digitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left Alternation,  4.Righ Alternation
%                                 5. Home Delay,        6. Is alternation forzed?
%                               If not called, look for it.
% OUTPUT
%       - armChoice.behaviour output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.arm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.ints          Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.delay.timestamps    Delay timestamps, in seconds
% armChoice.armAlternation      Performance vector, 1 is alternation, 0 is
%                                   repetition. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
%
%   MV-BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
% p = inputParser;
% addParameter(p,'fs',30,@isnumeric);
% addParameter(p,'artThr',10,@isnumeric);
% addParameter(p,'convFact',0.1149,@isnumeric);
% % maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% % 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 
% parse(p,varargin{:})
% fs = p.Results.fs;
% artThr = p.Results.artThr;
% convFact = p.Results.convFact;

if ~isempty(dir('*.ArmChoice.Events.mat'))
    disp('Arm choice already computed! Loading file.');
    file = dir('*.ArmChoice.Events.mat');
    load(file.name);
    return
end

%% Get basler TTL
% digital_in legend: 1. Basler, 2. maze LEd, 3. Left Alternation, 4.Righ
% Alternation, 5. Home Delay, 6. Is alternation forzed?
disp('Loading digital in...');
digitalIn = bz_getDigitalIn;
if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
    armChoice.timestamps = [digitalIn.timestampsOn{3}; digitalIn.timestampsOn{4}]; 
    % 0 is left, 1 is right
    armChoice.arm = [zeros(size(digitalIn.timestampsOn{3})); ones(size(digitalIn.timestampsOn{4}))];
    [armChoice.timestamps, idx] = sort(armChoice.timestamps);
    armChoice.arm = armChoice.arm(idx);
    armChoice.delay.ints = digitalIn.ints{5};
    armChoice.delay.timestamps = digitalIn.timestampsOn{5};
    armChoice.delay.dur = mean(armChoice.delay.ints(2,:) - armChoice.delay.ints(1,:));
    armChoice.armAlternation = [NaN; abs(diff(armChoice.arm))]; % 1 is alternation, 0 is no alternation
    armChoice.performance = nansum(armChoice.armAlternation)/(length(armChoice.armAlternation) - 1);
    
    if size(digitalIn.timestampsOn,2) >=6
        if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
            armChoice.forzed = 1;
            desc = 'forzed alternation';
        else
            armChoice.forzed = 0;
            desc = 'spontaneous alternation';
        end
    else
        if armChoice.performance == 1 
            armChoice.forzed = 1;
            desc = 'forzed alternation';
        else
            armChoice.forzed = 0;
            desc = 'spontaneous alternation';
        end
    end
    
    h = figure;
    hold on
    plot(armChoice.timestamps, armChoice.arm,'color',[.7 .7 .7]);
    scatter(armChoice.timestamps(isnan(armChoice.armAlternation)),...
        armChoice.arm(isnan(armChoice.armAlternation)),100,[.8 .8 .8],'filled');
    scatter(armChoice.timestamps(find(armChoice.armAlternation == 1)),...
        armChoice.arm(find(armChoice.armAlternation == 1)),100,[.6 .9 .7],'filled');
    scatter(armChoice.timestamps(find(armChoice.armAlternation == 0)),...
        armChoice.arm(find(armChoice.armAlternation == 0)),100,[.9 .6 .7],'filled');
    for ii = 1:size(armChoice.delay.timestamps,1)
        fill([armChoice.delay.ints(:,ii); flip(armChoice.delay.ints(:,ii))],[1 1 1.2 1.2]',...
        [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
    end
    xlabel('seconds'); ylim([-.2 1.2]);
    text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
        desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
        num2str(length(armChoice.arm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
        's'));
    set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});

    mkdir('Behavior');
    saveas(h,'Behavior\armChoice.png');

    C = strsplit(pwd,'\');
    save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
else
    warning('DigitalIn format does not match. Was T maze performed? ');
end
end

%     h1 = figure;
%     subplot(3,1,[1 2])
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     scatter(mazeVirtual(:,1), mazeVirtual(:,2),5,linearizeMazeVirtual);
%     ylabel('cm'); 
%     subplot(3,1,[3])
%     hold on
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMaze)),...
%         vlinMaze,5,vlinMaze);
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMazeCont)),...
%         vlinMazeCont,5,vlinMazeCont);
%     xlabel('cm');  ylabel('cm');
%     colormap jet
%     saveas(h1,'Behavior\virtualMaze.png');

%     centroidR = [30, 22]; 
%     centroidL = [30, 55];
%     x = behavior.positions.x';
%     y = behavior.positions.y';
%     t = behavior.timestamps;
%     [rangle,rradius] = cart2pol(x-centroidR(1), y - centroidR(2));
%     [langle,lradius] = cart2pol(x-centroidL(1), y - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + .15);
%     langle = wrapTo2Pi(-langle + 5.3658 + .15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial; % Distance traveled in cm, left trial after righ trial
%     
%     h1 = figure;
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     xlabel('cm'); ylabel('cm'); 
%     plot(centroidR(1),centroidR(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     plot(centroidL(1),centroidL(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     saveas(h1,'Behavior\centroids.png');
%     
%     % 
%     h2 = figure;
%     subplot(3,1,[1 2])
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     prev = 0;
%     linearized = []; arm = [];
%     for ii = 1:size(armChoice.delay.timestamps,1)
%         winTrial = [prev armChoice.delay.timestamps(ii)];
%         prev = armChoice.delay.timestamps(ii);
%         xspam = behavior.timestamps >= winTrial(1) & behavior.timestamps <= winTrial(2);
%         subplot(3,1,[1 2])
%         hold on
%         p = plot(x(xspam),y(xspam),'lineWidth',2);
%         ylabel('cm'); 
%         
%         if armChoice.arm(ii) == 1
%             linearized = [linearized rangle(xspam)];
%             arm = [arm ones(size(rangle(xspam)))];
%         elseif armChoice.arm(ii) == 0
%             linearized = [linearized langle(xspam)];
%             arm = [arm zeros(size(rangle(xspam)))]; 
%         end
%         subplot(3,1,3)
%         hold on
%         plot(t(xspam),linearized(xspam),'lineWidth',2);
%         xlabel('samples');
%         drawnow;
%         pause;    
%         delete(p);
%     end
%     % close(h2);
%       
% %     for ii = 1:size(behavior.timestamps)
% %     end    
% while 1
%     [xg,yg]=ginput(1);
%     fprintf('x: %3.2f, y: %3.2f \n',xg,yg);
%     
%     [rangle,rradius] = cart2pol(xg-centroidR(1), yg - centroidR(2));
%     [langle,lradius] = cart2pol(xg-centroidL(1), yg - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + 0.15);
%     langle = wrapTo2Pi(-langle + 5.3658 + 0.15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial ;
%     fprintf('linearized R: %3.2f \n',rangle);
%     fprintf('linearized L: %3.2f \n\n',langle);
% end
