function [placeFieldStats] = definePlaceFieldBoundaries(varargin)

% Parse inputs 
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMapsAvg',{},@isstruct);
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.50,@isnumeric);
addParameter(p,'minPeak',2,@isnumeric);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
firingMaps = p.Results.firingMapsAvg;
% Get session info
basename = bz_BasenameFromBasepath(basepath);
% load([basepath filesep basename '.sessionInfo.mat']);
% Default firingMapsAvg
if isempty(firingMaps)
    firingMaps = load([basepath filesep basename '.firingMapsAvg.cellinfo.mat']);
    firingMaps = firingMaps.firingMaps;
end
sizeMaze = length(firingMaps.rateMaps{1}{1});
threshold = p.Results.threshold;
minSize = p.Results.minSize * sizeMaze;
maxSize = p.Results.maxSize * sizeMaze;
minPeak = p.Results.minPeak;
saveMat = p.Results.saveMat;

%% Find place fields
for unit = 1:length(firingMaps.rateMaps)     
    for c = 1:length(firingMaps.rateMaps{1})
        % Default values
        mapStats{unit,1}{c}.x = NaN;
        mapStats{unit,1}{c}.field = [];
        mapStats{unit,1}{c}.size = [];
        mapStats{unit,1}{c}.fieldX = [NaN NaN];
        mapStats{unit,1}{c}.fieldPeak = [];
        mapStats{unit,1}{c}.peak = [];
        mapStats{unit,1}{c}.mean = [];
     
        
         % Determine the field as the convex area around the peak where the value or rate is > threshold*peak
        % There can be two or more fields
        z = firingMaps.rateMaps{unit}{c};
        x = 1:length(firingMaps.rateMaps{1}{1});
        
        % If there is no firing rate, go to next unit
        if max(z) == 0
          mapStats{unit,1}{c}.field = logical(zeros(size(z)));
          continue;
        end
        
        [peakValues, peakLocations] = findpeaks(z, 'minpeakheight',minPeak);
        mapStats{unit,1}{c}.field = logical(zeros(size(z)));
        for j = 1:length(peakLocations)
            FieldPeak = peakLocations(j);
            % FieldPeak must be 2 Hz or more
            if peakValues(j) < minPeak, continue, end
            LookForward = [FieldPeak+1:sizeMaze,1:FieldPeak-1];
            LookBack = fliplr(LookForward);
            PercentPeakRate_Forward = z(LookForward)./peakValues(j);
            PercentPeakRate_Back = z(LookBack)./peakValues(j);
            tempInd1 = find(PercentPeakRate_Forward < threshold);
            if isempty(tempInd1), continue, end
            FieldEnd = LookForward(tempInd1(1)); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
            tempInd2 = find(PercentPeakRate_Back < threshold);
            FieldStart = LookBack(tempInd2(1)); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
            % Field must be at least 10 bins and less than 40 bins in the length (more than 20 cm and less than 80cm)
            if FieldEnd>FieldStart && FieldEnd-FieldStart < minSize, continue, end
            if FieldEnd<FieldStart && FieldEnd+(sizeMaze-FieldStart) < minSize, continue, end
            if FieldEnd>FieldStart && FieldEnd-FieldStart > maxSize, continue, end
            if FieldEnd<FieldStart && FieldEnd+(sizeMaze-FieldStart) > maxSize, continue, end
            if FieldEnd<FieldStart
                FieldEnd = sizeMaze+FieldEnd;
            end
            mapStats{unit,1}{c}.field(FieldStart:FieldEnd) = 1;
            mapStats{unit,1}{c}.fieldPeak = FieldPeak;
            mapStats{unit,1}{c}.fieldX(j,:) = [FieldStart FieldEnd];
            mapStats{unit,1}{c}.size(j) = FieldEnd-FieldStart;
            mapStats{unit,1}{c}.peak(j) = peakValues(j);
            mapStats{unit,1}{c}.mean(j) = mean(z(mapStats{unit,1}{c}.field));
        end
    end
end

%% =================
%   WRITE OUTPUT    
% =================

placeFieldStats = {};

% inherit required fields from spikes cellinfo struct
placeFieldStats.UID = firingMaps.UID;
placeFieldStats.sessionName = firingMaps.sessionName;
try
placeFieldStats.region = firingMaps.region; 
catch
   %warning('spikes.region is missing') 
end

placeFieldStats.params.sizeMaze = sizeMaze;
placeFieldStats.params.threshold = threshold;
placeFieldStats.params.minSize = minSize;
placeFieldStats.params.maxSize = maxSize;
placeFieldStats.params.minPeak = minPeak;
placeFieldStats.params.saveMat = saveMat;

placeFieldStats.mapStats = mapStats;

if saveMat
   save([basepath,filesep,placeFieldStats.sessionName '.placeFields.cellinfo.mat'],'placeFieldStats'); 
end


%% ==========
%   PLOT    
% ==========
for unit = 1:length(firingMaps.rateMaps)
    figure;
    for c = 1:length(firingMaps.rateMaps{1})
        subplot(2,2,c)
        plot(firingMaps.rateMaps{unit}{c},'k')
        if sum(firingMaps.rateMaps{unit}{c})>0
            hold on
            for ii = 1:size(mapStats{unit}{c}.field,2)
                plot(find(mapStats{unit}{c}.field(:,ii)),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field(:,ii)==1),'linewidth',2)
            end
        end
        if c==1 || c==3, ylabel('FR(Hz)'); end
        if c>2, xlabel('Track (cm)'); end
        if c==1, title(['                                                                  Cell ' num2str(unit)]); end
        %ylim([0,12])
    end
    if ~isfolder([basepath,filesep,'newPC1s'])
        mkdir(basepath,'newPC1s')
    end
    saveas(gcf,[basepath,filesep,'newPC1s',filesep ,'cell_' num2str(unit) '.png'],'png');
    close all;
end

end
