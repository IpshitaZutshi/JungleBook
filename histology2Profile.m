function histology2Profile

%% List of animals to use
animalID = {'IZ11','IZ12','IZ13','IZ15','IZ18','IZ19','IZ20','IZ21','IZ23',...
    'IZ24','IZ25','IZ27','IZ28','IZ29','IZ30','IZ31'};
% Excluded animals with poor/incomplete profile = 'IZ17','IZ26'

%% Location of histology images
% (Update depending on where the folder is saved)
histfolder = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\Compiled\Histology\Track\';


%% Plot all the histology images first 
figure
set(gcf,'renderer','painters');
set(gcf,'Position',[100 100 1500 800])   

for ii = 1:length(animalID)
    subplot(4,4*3,[3*(ii-1)+1 3*(ii-1)+2]) 
    currImg = imread(strcat(histfolder,animalID{ii},'.png'));
    imshow(currImg)
    title(animalID{ii})
end

%% Now determine the scale. IMPORTANT - do not change the size of the figure after this, until after the program finishes!
disp('Draw a line over the scale bar to determine the scale. Double click when done.')
scalebar = 1; %in mm
subplot(4,4*3,[1 2])
roi = drawline;
linelength = pdist(roi.Position,'Euclidean');
scale = linelength/scalebar;

%% Next draw the probe track (can be approximate) for each image. 
% Also mark the point where the track and CA1 pyramidal cell layer intersect.
for ii = 1:length(animalID)
    subplot(4,4*3,[3*(ii-1)+1 3*(ii-1)+2]) 
    disp('Next, draw along the probe track')
    roiProbe{ii} = drawline('Color','w');
    disp('Next, mark the point where the probe and CA1 intersect')
    roiPyr{ii} = drawpoint;
    
    %% Extract the profile file and determine number of channels between s.p, s.r., andn s.l.m.
    % load(strcat(histfolder,animalID{ii},'.mat'))
    
    % Convert physical distance to the scale on the images
    if strcmp(animalID{ii},'IZ11')==1 || strcmp(animalID{ii},'IZ15')==1 || strcmp(animalID{ii},'IZ23')==1 || ...
            strcmp(animalID{ii},'IZ30')==1 || strcmp(animalID{ii},'IZ29')==1
        distbetweenCh = 0.1; %in mm
    else
        distbetweenCh = 0.02; %in mm
    end
    
    %% Find the intersection point of the s.r. and s.l.m on the drawn probe track
    
    
    %% Draw those points on the images
    
    
    %% Extract theta profile and plot it next to the image
    % subplot(4,4*3,[3*(ii-1)+3]) 
    
    % Mark s.p., s.r. and s.l.m. on the profile 
    
end

end
