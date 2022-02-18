function plotHistologyRegion

%% List of animals to use
animalID = {'IZ12','IZ13','IZ15', 'IZ18','IZ19','IZ20','IZ21','IZ23',...
    'IZ24','IZ25','IZ26','IZ27','IZ28','IZ29','IZ30','IZ31','IZ32', 'IZ33'};
% Excluded animals with poor/incomplete profile = 'IZ17','IZ26','IZ11',,'IZ34'

%% Location of histology images
% (Update depending on where the folder is saved)
histfolder = 'C:\Users\Ipshita\Desktop\TrackData\';
% sr_channel_dist = [2, 7, 7,1, 6, 7, 10, 10, 2, 7, 7, 8, 10, 3, 2, 8, 2, 11]
%     slm_channel_dist = [4,17,17,5,15,17,19,21,3,16,18,17,18,4,4,17,4,23]
radii_list = [];
slope_list = [];
sr_plot_list = [];
slm_plot_list = [];

%% Plot all the histology images first 
% figure
% set(gcf,'renderer','painters');
% set(gcf,'Position',[100 100 1500 800])   

% for ii = 1:length(animalID)
%     subplot(4,4*3,[3*(ii-1)+1 3*(ii-1)+2]) 
%     currImg = imread(strcat(histfolder,animalID{ii},'.png'));
%     imshow(currImg(:,:,3))
%     title(animalID{ii})
% end

for ii = 1:length(animalID)
    
    figure(1)
    set(gcf,'Position',[10 10 800 1200])
    currImg = imread(strcat(histfolder,animalID{ii},'.png'));
    imshow(currImg)
    hold on
    title(animalID{ii})
    
    if ii == 1
        disp('Draw a line over the scale bar to determine the scale. Double click when done.')
        scalebar = 1; %in mm
        roi = drawline;
        linelength = pdist(roi.Position,'Euclidean');
        scale = linelength/scalebar;
    end
    disp('Next, draw along the probe track')
    roiProbe{ii} = drawline('Color','w','LineWidth', 1);
    disp('Next, mark the point where the probe and CA1 intersect')
    roiPyr{ii} = drawpoint;
    disp('Next, mark your best estimate of the mid point of sr')
    roisr{ii} = drawpoint;
    disp('Next, mark your best estimate of the fissure')
    roifiss{ii} = drawpoint;    
    
    %% Extract the profile file and determine number of channels between s.p, s.r., andn s.l.m.
    % Load the pyr, rad and slm channels
    load(strcat(histfolder,animalID{ii},'.mat'));
    channelOrder = hippocampalLayers.channelOrder;
    idxPyr = find(channelOrder==hippocampalLayers.pyramidal);
    idxrad = find(channelOrder==hippocampalLayers.radiatum);
    idxslm = find(channelOrder==hippocampalLayers.slm);
    sr_channel_dist = idxrad-idxPyr;
    slm_channel_dist = idxslm-idxPyr;    
    
    % Convert physical distance to the scale on the images
    if strcmp(animalID{ii},'IZ11')==1 || strcmp(animalID{ii},'IZ15')==1 || strcmp(animalID{ii},'IZ23')==1 || ...
            strcmp(animalID{ii},'IZ30')==1 || strcmp(animalID{ii},'IZ29')==1 || strcmp(animalID{ii},'IZ32')==1
        distbetweenCh = 0.1; %in mm
    else
        distbetweenCh = 0.02; %in mm
    end
        
    %% Now perform the computation
    sr_radius = distbetweenCh * sr_channel_dist * scale(1,1); 
    slm_radius = distbetweenCh * slm_channel_dist * scale(1,1); 
    radii_list(1,ii) = sr_radius;
    radii_list(2,ii) = slm_radius;
    % find tan 
    slope = (roiProbe{1,ii}.Position(2,2) - roiProbe{1,ii}.Position(1,2)) / (roiProbe{1,ii}.Position(2,1) - roiProbe{1,ii}.Position(1,1)); % %
    
    slope_list(ii) = slope;
    theta = atan(abs(slope));
    sr_cos = sr_radius * cos(theta); %xdist
    sr_sin = sr_radius * sin(theta); %ydist
    slm_cos =  slm_radius * cos(theta);
    slm_sin = slm_radius * sin(theta);
    
    %% Find the projection of sp, sr and fiss onto the probe track
    vector = roiProbe{1,ii}.Position;
    projPointPyr = proj(vector, roiPyr{1,ii}.Position);
    projPointsr = proj(vector, roisr{1,ii}.Position);
    projPointfiss = proj(vector, roifiss{1,ii}.Position);
    
    %% Find the intersection point of the s.r. and s.l.m on the drawn probe track
    if slope < 0 % x2<x1
        sr_plot_list(1,ii) = projPointPyr(1) - sr_cos; 
        sr_plot_list(2,ii) = projPointPyr(2) + sr_sin;
        slm_plot_list(1,ii) =  projPointPyr(1) - slm_cos; 
        slm_plot_list(2,ii) = projPointPyr(2) + slm_sin;
    else 
        sr_plot_list(1,ii) = projPointPyr(1) + sr_cos; 
        sr_plot_list(2,ii) = projPointPyr(2) + sr_sin;
        slm_plot_list(1,ii) =  projPointPyr(1) + slm_cos; 
        slm_plot_list(2,ii) = projPointPyr(2) + slm_sin;
    end
    
    %% Distance between determined points and estimated points.
    distsr(ii) = pdist([sr_plot_list(:,ii)';projPointsr'])*scale*1000;
    if sr_plot_list(2,ii)>=projPointsr(2) % Point is above
        distsr(ii) = distsr(ii)*-1;
    end

    distslm(ii) = pdist([slm_plot_list(:,ii)';projPointfiss'])*scale*1000;
    if slm_plot_list(2,ii)>=projPointfiss(2) % Point is above
        distslm(ii) = distslm(ii)*-1;
    end
        
    
%     else slope > 0
%         sr_plot_list(1,ii) = roiPyr{1,ii}.Position(1,1) - sr_cos 
%         sr_plot_list(2,ii) = roiPyr{1,ii}.Position(1,2) + sr_sin
%         slm_plot_list(1,ii) = roiPyr{1,ii}.Position(1,1) - slm_cos 
%         slm_plot_list(2,ii) = roiPyr{1,ii}.Position(1,2) + slm_sin
%     end
     
    % Now replot everything
    figure(2)
    set(gcf,'Renderer','painters')
    imshow(currImg(:,:,:));
    hold on;
    lh = line(vector(:,1), vector(:,2),'LineWidth',1.5);
    lh.Color=[1,1,1,0.2];
    axis on
    plot(sr_plot_list(1,ii), sr_plot_list(2,ii), 'co','MarkerFaceColor','c','MarkerSize',7) 
    plot(slm_plot_list(1,ii), slm_plot_list(2,ii), 'co','MarkerFaceColor','c','MarkerSize',7)
    
    plot(projPointPyr(1),projPointPyr(2),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.9763,0.9831,0.0538],'MarkerSize',5)
    plot(projPointsr(1),projPointsr(2),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.9763,0.9831,0.0538],'MarkerSize',5)
    plot(projPointfiss(1),projPointfiss(2),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.9763,0.9831,0.0538],'MarkerSize',5)

    saveas(figure(2),strcat(histfolder,'\Images\',animalID{ii},'_new.png'))
    saveas(figure(2),strcat(histfolder,'\Images\',animalID{ii},'_new.fig'))
    saveas(figure(2),strcat(histfolder,'\Images\',animalID{ii},'_new.eps'),'epsc')   
    close all
end
figure
set(gcf,'Renderer','painters')
stats = groupStats([distsr distslm],[ones(1,length(distsr)) ones(1,length(distslm))+1],'plotData',true,'repeatedMeasures',true);
[stats.signrank{1}.p,~,stats.signrank{1}.h] = signrank(distsr);
[stats.signrank{2}.p,~,stats.signrank{2}.h] = signrank(distslm);
title(strcat(num2str(stats.signrank{1}.p),'-----',num2str(stats.signrank{2}.p)))
Distance.distsr = distsr;
Distance.dislm = distslm;
save(strcat(histfolder,'\Images\SummaryStatsFinal.mat'),'Distance','stats')
saveas(gcf,strcat(histfolder,'\Images\Summary.png'))
saveas(gcf,strcat(histfolder,'\Images\Summary.fig'))
saveas(gcf,strcat(histfolder,'\Images\Summary.eps'),'epsc')  

end

function [ProjPoint] = proj(vector, q)
p0 = vector(1,:);
p1 = vector(2,:);
a = [-q(1)*(p1(1)-p0(1)) - q(2)*(p1(2)-p0(2)); ...
    -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))]; 
b = [p1(1) - p0(1), p1(2) - p0(2);...
    p0(2) - p1(2), p1(1) - p0(1)];
ProjPoint = -(b\a);
end
    

figure(1)
tcl = tiledlayout(5,3);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
set(gcf,'Position',[10 10 900 1500])
set(gcf,'Renderer','painters')
animalID = {'IZ12','IZ13','IZ15', 'IZ18','IZ19','IZ20','IZ21','IZ23',...
     'IZ24','IZ27','IZ28','IZ30','IZ31','IZ32', 'IZ33'};%'IZ25','IZ26','IZ29',
histfolder = 'C:\Users\Ipshita\Desktop\TrackData\Images\';
for ii = 1:length(animalID)   
%     subplot(6,3,ii)
    f1 = openfig(strcat(strcat(histfolder,animalID{ii},'_new.fig')))
    ax1=gca;
    ax1.Parent=tcl;
    ax1.Layout.Tile=ii;
%     currImg = imread(strcat(histfolder,animalID{ii},'_new.png'));
%     currImg = imcrop(currImg,[135.51,47.51,1223.98,643.98]);
%     imshow(currImg)
%     axis off
    hold on
    title(animalID{ii})
    nexttile
end


% open(strcat(strcat(histfolder,animalID{ii},'_new.fig')))
% axis off
% saveas(gcf,strcat(histfolder,animalID{ii},'_new.png'))
% saveas(gcf,strcat(histfolder,animalID{ii},'_new.fig'))
% saveas(gcf,strcat(histfolder,animalID{ii},'_new.eps'),'epsc')   
% close all
    %% Extract theta profile and plot it next to the image
    % subplot(4,4*3,[3*(ii-1)+3]) 
    
    % Mark s.p., s.r. and s.l.m. on the profile 
% end
%% Now determine the scale. IMPORTANT - do not change the size of the figure after this, until after the program finishes!
% disp('Draw a line over the scale bar to determine the scale. Double click when done.')
% scalebar = 1; %in mm
% subplot(2,2*3,[1 2])
% roi = drawline;
% linelength = pdist(roi.Position,'Euclidean');
% scale = linelength/scalebar;

